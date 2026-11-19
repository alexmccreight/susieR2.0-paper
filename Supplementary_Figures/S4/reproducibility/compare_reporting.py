import ast
from contextlib import redirect_stdout
import io
from pathlib import Path
import sys

import numpy as np
import pandas as pd


def sparsepro_class(path):
    # Load the original class without running the file's top-level CLI.
    tree = ast.parse(path.read_text())
    tree.body = [
        node for node in tree.body
        if isinstance(node, (ast.Import, ast.ImportFrom))
        or (isinstance(node, ast.ClassDef) and node.name == "SparsePro")
    ]
    module = {}
    exec(compile(tree, str(path), "exec"), module)
    return module["SparsePro"]


def pair_sets(work):
    keys = ["case", "K", "strength", "effect"]
    posteriors = pd.read_csv(work / "posteriors.csv")
    sparsepro = sparsepro_class(work / "sparsepro_zld.py")
    reported = []
    for case, fit in posteriors.groupby("case"):
        alpha = fit.filter(regex=r"^alpha_").to_numpy()
        effects, p = alpha.shape
        model = sparsepro(p, effects, np.ones(p), 0, 1)
        model.gamma = alpha.T
        with redirect_stdout(io.StringIO()):
            sets, *_ = model.get_effect()
        truth = np.array(str(fit["truth"].iloc[0]).split(), dtype=int) - 1
        for effect, members in sets.items():
            if len(members):
                reported.append((
                    case, fit["K"].iloc[0], fit["strength"].iloc[0], effect + 1,
                    alpha[effect, members].sum(), len(members), np.isin(members, truth).any(),
                ))
    columns = keys + ["mass", "size", "contains_truth"]
    paired = pd.DataFrame(reported, columns=columns).merge(
        pd.read_csv(work / "susier-sets.csv"), on=keys,
        suffixes=("_sparsepro", "_susier"), validate="one_to_one",
    )
    print(f"{posteriors['case'].nunique()} converged fits, {len(paired)} paired sets")
    return paired.sort_values(keys).reset_index(drop=True)


if __name__ == "__main__":
    work, reference = map(Path, sys.argv[1:])
    paired = pair_sets(work)
    pd.testing.assert_frame_equal(
        paired, pd.read_csv(reference), check_exact=False, rtol=1e-10, atol=1e-12
    )
    paired.to_csv(work / "paired-posterior-mass.csv", index=False)
    for column in ("mass_sparsepro", "mass_susier"):
        below = (paired[column] < .95 - 1e-12).sum()
        print(f"{column}: {below}/{len(paired)} below 95%")
