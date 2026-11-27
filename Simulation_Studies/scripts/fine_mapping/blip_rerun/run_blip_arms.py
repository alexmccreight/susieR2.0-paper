#!/usr/bin/env python3
"""
BLiP re-run, step 2: SuSiE+BLiP arms on the exported replicates
================================================================
Reads the folders written by export_blip_inputs.R
(<input_dir>/K<k>/rep_<id>/{alpha.npy, corr_sub.npy, meta.json}) and writes
blip_arms.json next to them: one group list per arm and per q.

Every arm follows the BLiP authors' UK Biobank wrapper
(github.com/amspector100/ukbb_blip, blip_wrapper.py) and calls
pyblip.blip.BLiP(error='fdr', weight_fn='inverse_size', deterministic=True):

  A_capped            SNPs with PIP >= --pip-threshold are retained; candidate
                      groups are pyblip sequential groups over the retained
                      alphas (size <= --max-size, PEP < --max-pep) plus SuSiE's
                      own per-component sets formed on the full alphas (a
                      component is kept when its top --max-size alphas exceed
                      1 - q; PEP = 1 - that sum) remapped into the retained
                      index; then BLiP at level q with --max-pep.
  B_uncapped          A with the size cap lifted (--uncapped-max-size).
  C_hier_capped       A / B plus hierarchical-clustering groups (single,
  C_hier_uncapped     average and complete linkage on 1 - |r| over the retained
                      SNPs, as pyblip.susie_groups does when given X), with
                      PEP = prod_l (1 - sum_{j in G} alpha_lj).
  D_purity_capped     A / B with candidate groups whose min |r| < 0.5 removed
  D_purity_uncapped   before BLiP (pyblip.create_groups.filter_by_purity).
  D_purity08_uncapped B with the same filter at min |r| >= 0.8 (the analogue
                      of SuSiE's purity-0.8 sets).

Indices in the output are 1-based, matching meta.json.

Examples
  HPC, one job per replicate chunk:
    python run_blip_arms.py --input-dir /home/apm2217/output/blip_rerun \
        --k 1 --start-replicate 1 --end-replicate 25
  Local, every K, 16 worker processes:
    python run_blip_arms.py --input-dir <export root> --n-jobs 16

A replicate whose blip_arms.json was produced with the same settings is
skipped, so a run can be resumed after an interruption.
"""

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import copy
import glob
import json
import platform
import re
import sys
import time
import traceback
from datetime import datetime
from multiprocessing import Pool

import numpy as np
import scipy
import cvxopt
import cvxpy
import pyblip
from pyblip import blip as pyblip_blip
from pyblip.create_groups import (MIN_PEP, CandidateGroup,
                                  _dist_matrices_to_groups, filter_by_purity,
                                  sequential_groups)

ARMS = {
    "A_capped":          dict(capped=True,  hier=False, min_purity=0.0),
    "B_uncapped":        dict(capped=False, hier=False, min_purity=0.0),
    "C_hier_capped":     dict(capped=True,  hier=True,  min_purity=0.0),
    "C_hier_uncapped":   dict(capped=False, hier=True,  min_purity=0.0),
    "D_purity_capped":   dict(capped=True,  hier=False, min_purity=0.5),
    "D_purity_uncapped": dict(capped=False, hier=False, min_purity=0.5),
    "D_purity08_uncapped": dict(capped=False, hier=False, min_purity=0.8),
}


def q_key(q):
    return f"q{q:g}"


def coverage_key(q):
    """Name of the exported SuSiE CS list paired with level q (0.05 -> '0.95')."""
    return f"{1 - q:g}"


def versions():
    return dict(pyblip=pyblip.__version__, cvxpy=cvxpy.__version__,
                cvxopt=cvxopt.__version__, numpy=np.__version__,
                scipy=scipy.__version__, python=platform.python_version(),
                installed_solvers=cvxpy.installed_solvers())


def settings_dict(opts):
    return dict(pip_threshold=opts.pip_threshold, max_pep=opts.max_pep,
                max_size=opts.max_size, uncapped_max_size=opts.uncapped_max_size,
                prenarrow=opts.prenarrow, solver=opts.solver,
                bin_solver=opts.bin_solver, weight_fn="inverse_size",
                error="fdr", deterministic=True, q_values=list(opts.q),
                arms=list(opts.arms))


def jsonable(obj):
    if isinstance(obj, dict):
        return {k: jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [jsonable(v) for v in obj]
    if isinstance(obj, (np.bool_, bool)):
        return bool(obj)
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    return obj


# ---- Candidate groups ----

def susie_own_sets(alpha, q, max_size):
    """SuSiE's per-component sets as formed in the authors' wrapper
    (susie_results): sort a component's alphas, cumulate the top max_size,
    keep the component when that sum exceeds 1 - q, report the smallest top-k
    set reaching it with PEP = 1 - cumsum. Indices are 0-based here."""
    sets = []
    for l in range(alpha.shape[0]):
        order = np.argsort(-alpha[l], kind="stable")[:max_size]
        cums = np.cumsum(alpha[l, order])
        hits = np.where(cums > 1 - q)[0]
        if hits.size == 0:
            continue
        k = int(hits.min())
        sets.append(dict(component=l + 1, group=order[:k + 1].tolist(),
                         pep=float(1 - cums[k])))
    return sets


def product_pep(alpha_rel, group):
    """P(no signal in group) under the SuSiE posterior."""
    iter_peps = 1 - alpha_rel[:, list(group)].sum(axis=1)
    iter_peps[iter_peps < MIN_PEP] = MIN_PEP
    return float(np.exp(np.log(iter_peps).sum()))


def is_contiguous(group):
    g = sorted(group)
    return g[-1] - g[0] == len(g) - 1


def sub_corr(corr_sub, pos, orig_idx):
    ii = np.array([pos[int(j)] for j in orig_idx], dtype=int)
    return corr_sub[np.ix_(ii, ii)]


def uncapped_sizes(opts, nrel, p, largest_cs):
    """(sequential max_size, SuSiE-own-set max_size) for the uncapped arms."""
    if opts.uncapped_max_size == "auto":
        m = max(100, largest_cs)
        return min(m, nrel), m
    m = int(opts.uncapped_max_size)
    if m <= 0:
        return nrel, p
    return min(m, nrel), m


def run_arm(alpha, retained, corr_sub, pos, q, cfg, opts, largest_cs, seed,
            set_purity=None):
    """set_purity maps a SuSiE component to the exported purity of its set at
    the paired coverage; purity arms use it to drop impure own sets, whose
    members the exported correlation block does not cover."""
    L, p = alpha.shape
    nrel = int(len(retained))
    t0 = time.time()
    set_purity = set_purity or {}
    if cfg["capped"]:
        seq_max, own_max = opts.max_size, opts.max_size
    else:
        seq_max, own_max = uncapped_sizes(opts, nrel, p, largest_cs)

    result = dict(groups=[], peps=[], n_retained=nrel, n_cand_groups=0,
                  n_own_sets=0, n_own_added=0, n_own_dropped_impure=0,
                  n_own_skipped_nocorr=0, n_hier_added=0,
                  seq_max_size=int(seq_max), own_max_size=int(own_max),
                  own_sets=[], problem_status={})
    if nrel == 0:
        result["elapsed_seconds"] = time.time() - t0
        return result

    alpha_rel = alpha[:, retained]
    cand = sequential_groups(susie_alphas=alpha_rel, q=q, max_pep=opts.max_pep,
                             max_size=max(int(seq_max), 1),
                             prenarrow=opts.prenarrow)

    # SuSiE's own sets on the full alphas, remapped into the retained index.
    # SNPs outside the retained set are appended after it, as in the wrapper.
    orig2new = {int(j): i for i, j in enumerate(retained)}
    appended = []
    own = susie_own_sets(alpha, q, own_max)
    n_own_added = n_own_dropped_impure = n_own_skipped_nocorr = 0
    for s in own:
        if cfg["min_purity"] > 0:
            # An impure set can never pass the filter; a set the correlation
            # block does not cover cannot be evaluated and is skipped.
            pur = set_purity.get(s["component"])
            if pur is not None and pur < cfg["min_purity"]:
                n_own_dropped_impure += 1
                continue
            if any(int(j) not in pos for j in s["group"]):
                n_own_skipped_nocorr += 1
                continue
        g = []
        for j in s["group"]:
            if j not in orig2new:
                orig2new[j] = nrel + len(appended)
                appended.append(j)
            g.append(orig2new[j])
        # A contiguous window of retained SNPs within the sequential size cap is
        # already a sequential group (with a PEP at most as large); anything
        # else, in particular a set reaching into the appended SNPs, is added.
        if all(j < nrel for j in g) and len(g) <= seq_max and is_contiguous(g):
            continue
        cand.append(CandidateGroup(group=set(g), pep=s["pep"]))
        n_own_added += 1
    new2orig = np.concatenate([np.asarray(retained, dtype=int),
                               np.asarray(appended, dtype=int)])

    n_hier_added = 0
    if cfg["hier"] and nrel >= 2:
        dist = np.abs(1 - sub_corr(corr_sub, pos, retained))
        for g in _dist_matrices_to_groups(dist.copy()):
            if len(g) > seq_max or is_contiguous(g):
                continue
            pep = product_pep(alpha_rel, g)
            if pep < opts.max_pep:
                cand.append(CandidateGroup(group=set(int(j) for j in g), pep=pep))
                n_hier_added += 1

    if cfg["min_purity"] > 0:
        corr_new = sub_corr(corr_sub, pos, new2orig)
        cand = filter_by_purity(cand, min_purity=cfg["min_purity"],
                                corr_matrix=corr_new)

    n_cand = len(cand)
    detections, status = [], {}
    if n_cand > 0:
        np.random.seed(seed)   # BLiP perturbs its weights with np.random
        detections, status = pyblip_blip.BLiP(
            cand_groups=copy.deepcopy(cand), weight_fn="inverse_size",
            error="fdr", q=q, max_pep=opts.max_pep, deterministic=True,
            verbose=False, solver=opts.solver, return_problem_status=True)

    result.update(
        groups=[sorted(int(new2orig[j]) + 1 for j in d.group) for d in detections],
        peps=[float(d.pep) for d in detections],
        n_cand_groups=int(n_cand), n_own_sets=len(own),
        n_own_added=n_own_added, n_own_dropped_impure=n_own_dropped_impure,
        n_own_skipped_nocorr=n_own_skipped_nocorr, n_hier_added=n_hier_added,
        own_sets=[dict(component=s["component"],
                       group=[int(j) + 1 for j in s["group"]], pep=s["pep"])
                  for s in own],
        problem_status=jsonable(status),
        elapsed_seconds=time.time() - t0,
    )
    return result


# ---- Per-replicate driver ----

def process_replicate(task):
    rep_dir, opts = task
    out_path = os.path.join(rep_dir, "blip_arms.json")
    settings = settings_dict(opts)
    if os.path.exists(out_path) and not opts.overwrite:
        try:
            with open(out_path) as fh:
                if json.load(fh).get("settings") == settings:
                    return rep_dir, "exists", 0.0, ""
        except (OSError, ValueError):
            pass

    t0 = time.time()
    try:
        with open(os.path.join(rep_dir, "meta.json")) as fh:
            meta = json.load(fh)
        if int(meta.get("index_base", 1)) != 1:
            raise ValueError("meta.json indices must be 1-based")
        alpha = np.load(os.path.join(rep_dir, "alpha.npy"))
        corr_sub = np.load(os.path.join(rep_dir, "corr_sub.npy")).astype(np.float64)
        pip = np.asarray(meta["pip"], dtype=float)
        corr_index = np.asarray(meta["corr_index"], dtype=int) - 1
        pos = {int(j): i for i, j in enumerate(corr_index)}

        retained = np.where(pip >= opts.pip_threshold)[0]
        missing = [int(j) for j in retained if int(j) not in pos]
        if missing:
            raise ValueError(
                f"{len(missing)} retained SNPs lie outside corr_sub.npy; re-run "
                f"export_blip_inputs.R with pip_floor <= {opts.pip_threshold}")
        largest = {cov: max([len(e["cs"]) for e in entries] or [0])
                   for cov, entries in meta.get("susie_cs", {}).items()}
        set_purity = {cov: {int(e["component"]): float(e["min_abs_corr"])
                            for e in entries if e.get("min_abs_corr") is not None}
                      for cov, entries in meta.get("susie_cs", {}).items()}

        arms_out = {}
        for arm in opts.arms:
            arms_out[arm] = {}
            for q in opts.q:
                arms_out[arm][q_key(q)] = run_arm(
                    alpha, retained, corr_sub, pos, q, ARMS[arm], opts,
                    largest.get(coverage_key(q), 0), seed=int(meta["replicate_id"]),
                    set_purity=set_purity.get(coverage_key(q), {}))

        out = dict(replicate_id=meta["replicate_id"], K=meta["K"], p=meta["p"],
                   L=meta["L"], index_base=1, n_retained=int(len(retained)),
                   largest_susie_cs=largest, settings=settings,
                   versions=versions(), arms=arms_out,
                   elapsed_seconds=time.time() - t0,
                   timestamp=datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        tmp = out_path + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(out, fh)
        os.replace(tmp, out_path)
        return rep_dir, "ok", time.time() - t0, ""
    except Exception:
        return rep_dir, "error", time.time() - t0, traceback.format_exc()


def find_replicates(input_dir, ks, start, end):
    reps = []
    for meta_path in sorted(glob.glob(os.path.join(input_dir, "K*", "rep_*", "meta.json"))):
        rep_dir = os.path.dirname(meta_path)
        mk = re.fullmatch(r"K(\d+)", os.path.basename(os.path.dirname(rep_dir)))
        mr = re.fullmatch(r"rep_(\d+)", os.path.basename(rep_dir))
        if not mk or not mr:
            continue
        k, rid = int(mk.group(1)), int(mr.group(1))
        if ks and k not in ks:
            continue
        if start is not None and rid < start:
            continue
        if end is not None and rid > end:
            continue
        reps.append(rep_dir)
    return reps


def _init_worker(bin_solver):
    pyblip_blip.BIN_SOLVER = bin_solver


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input-dir", required=True, help="root written by export_blip_inputs.R")
    ap.add_argument("--k", type=int, nargs="*", default=None, help="K values to run (default: all found)")
    ap.add_argument("--start-replicate", type=int, default=None)
    ap.add_argument("--end-replicate", type=int, default=None)
    ap.add_argument("--q", type=float, nargs="+", default=[0.05, 0.10], help="FDR levels (default 0.05 0.10)")
    ap.add_argument("--arms", nargs="+", default=list(ARMS), choices=list(ARMS))
    ap.add_argument("--pip-threshold", type=float, default=0.01,
                    help="retain SNPs with PIP >= this (authors: 1 - PIP <= 0.99)")
    ap.add_argument("--max-pep", type=float, default=0.25,
                    help="max PEP for candidate groups and inside BLiP")
    ap.add_argument("--max-size", type=int, default=25, help="group size cap for the capped arms")
    ap.add_argument("--uncapped-max-size", default="auto",
                    help="'auto' = max(100, largest SuSiE CS at the paired coverage); "
                         "0 = no cap at all; or an integer")
    ap.add_argument("--prenarrow", action="store_true", help="pyblip prenarrowing (authors: off)")
    ap.add_argument("--solver", default="GLPK", help="cvxpy LP solver for BLiP")
    ap.add_argument("--bin-solver", default="GLPK_MI", help="cvxpy MILP solver for binarization")
    ap.add_argument("--n-jobs", type=int, default=1, help="worker processes")
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--quiet", action="store_true")
    opts = ap.parse_args(argv)

    if opts.uncapped_max_size != "auto":
        try:
            int(opts.uncapped_max_size)
        except ValueError:
            ap.error("--uncapped-max-size must be 'auto' or an integer")
    installed = cvxpy.installed_solvers()
    for s in (opts.solver, opts.bin_solver):
        if s not in installed:
            ap.error(f"solver {s} is not installed in cvxpy (installed: {installed})")
    pyblip_blip.BIN_SOLVER = opts.bin_solver

    reps = find_replicates(opts.input_dir, opts.k, opts.start_replicate, opts.end_replicate)
    print("\n========================================")
    print("BLiP re-run: SuSiE+BLiP arms")
    print("========================================")
    print("Input dir:  ", opts.input_dir)
    print("Replicates: ", len(reps))
    print("Settings:   ", json.dumps(settings_dict(opts)))
    print("Versions:   ", json.dumps(versions()))
    sys.stdout.flush()
    if not reps:
        print("Nothing to do.")
        return 0

    job_start = time.time()
    counts = {"ok": 0, "exists": 0, "error": 0}
    tasks = [(r, opts) for r in reps]

    def report(res):
        rep_dir, status, elapsed, msg = res
        counts[status] += 1
        if not opts.quiet or status == "error":
            label = os.path.join(os.path.basename(os.path.dirname(rep_dir)), os.path.basename(rep_dir))
            print(f"  {label}: {status} ({elapsed:.1f}s)")
            if msg:
                print(msg)
            sys.stdout.flush()

    if opts.n_jobs > 1:
        with Pool(opts.n_jobs, initializer=_init_worker, initargs=(opts.bin_solver,)) as pool:
            for res in pool.imap_unordered(process_replicate, tasks):
                report(res)
    else:
        for task in tasks:
            report(process_replicate(task))

    print(f"\nDone in {(time.time() - job_start) / 60:.1f} min: "
          f"{counts['ok']} run, {counts['exists']} already done, {counts['error']} failed.")
    return 1 if counts["error"] else 0


if __name__ == "__main__":
    sys.exit(main())
