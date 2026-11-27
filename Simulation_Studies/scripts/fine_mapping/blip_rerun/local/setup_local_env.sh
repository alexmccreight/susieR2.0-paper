#!/usr/bin/env bash
# Supplementary Figure S7, SuSiE+BLiP re-run: private Python and R
# environments. The system R library and Python are left untouched.
#   $BLIP_ENV_DIR/venv   Python with pyblip 0.4.2, cvxopt (the GLPK and GLPK_MI
#                        solvers for cvxpy), cvxpy, numpy, scipy, networkx,
#                        tqdm. Created by uv with Python PY_VERSION (3.12) when
#                        uv is installed, otherwise by python3 -m venv.
#   $BLIP_ENV_DIR/Rlib   susieR built from the development checkout at
#                        SUSIER_ROOT (default: ../susieR next to this
#                        repository, as in R/paths.R), with its Imports and
#                        LinkingTo dependencies, jsonlite and dplyr installed
#                        from CRAN when missing
#
# Usage: bash setup_local_env.sh [env dir]      (default $HOME/.susieR2_blip_env)
# Safe to re-run: susieR is rebuilt from source each time. Also works on a
# cluster login node, given a susieR checkout there.
set -euo pipefail

BLIP_ENV_DIR="${1:-${BLIP_ENV_DIR:-$HOME/.susieR2_blip_env}}"
PY_VERSION="${PY_VERSION:-3.12}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../../../.." && pwd)"
SUSIER_ROOT="${SUSIER_ROOT:-$(dirname "$REPO")/susieR}"
[ -f "$SUSIER_ROOT/DESCRIPTION" ] || {
  echo "susieR source not found at $SUSIER_ROOT; set SUSIER_ROOT to the checkout" >&2; exit 1; }
PY_PKGS=("pyblip==0.4.2" cvxopt cvxpy numpy scipy networkx tqdm)
mkdir -p "$BLIP_ENV_DIR"

VENV="$BLIP_ENV_DIR/venv"
if [ ! -x "$VENV/bin/python" ]; then
  if command -v uv >/dev/null 2>&1; then
    uv venv --python "$PY_VERSION" "$VENV"
  else
    python3 -m venv "$VENV"
  fi
fi
if command -v uv >/dev/null 2>&1; then
  uv pip install --python "$VENV/bin/python" "${PY_PKGS[@]}"
else
  "$VENV/bin/python" -m pip install --upgrade pip
  "$VENV/bin/python" -m pip install "${PY_PKGS[@]}"
fi
"$VENV/bin/python" - <<'PY'
import cvxopt, cvxpy, numpy, pyblip, scipy
print(f"pyblip {pyblip.__version__} | cvxpy {cvxpy.__version__} | cvxopt {cvxopt.__version__} "
      f"| numpy {numpy.__version__} | scipy {scipy.__version__}")
missing = [s for s in ("GLPK", "GLPK_MI") if s not in cvxpy.installed_solvers()]
if missing:
    raise SystemExit(f"cvxpy does not see {missing}; is cvxopt installed?")
PY

RLIB="$BLIP_ENV_DIR/Rlib"
mkdir -p "$RLIB"
# Dependencies of the local susieR (Imports + LinkingTo) and of the pipeline.
SUSIER_ROOT="$SUSIER_ROOT" R_LIBS_USER="$RLIB" Rscript -e '
  lib <- Sys.getenv("R_LIBS_USER")
  .libPaths(c(lib, .libPaths()))
  d <- read.dcf(file.path(Sys.getenv("SUSIER_ROOT"), "DESCRIPTION"))
  deps <- unlist(lapply(c("Imports", "LinkingTo"), function(f)
    if (f %in% colnames(d) && !is.na(d[1, f])) strsplit(d[1, f], ",")[[1]] else character(0)))
  deps <- trimws(gsub("[(].*[)]", "", deps))
  need <- setdiff(unique(c(deps[nzchar(deps)], "jsonlite", "dplyr")), rownames(installed.packages(priority = "base")))
  have <- vapply(need, function(p) requireNamespace(p, quietly = TRUE), logical(1))
  if (any(!have)) install.packages(need[!have], lib = lib, repos = "https://cloud.r-project.org")'
R_LIBS_USER="$RLIB" R CMD INSTALL --no-docs --no-multiarch --library="$RLIB" "$SUSIER_ROOT"
R_LIBS_USER="$RLIB" Rscript -e '
  .libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
  for (p in c("susieR", "jsonlite", "dplyr"))
    cat(sprintf("%s %s (%s)\n", p, as.character(packageVersion(p)), dirname(find.package(p))))'

echo "Environment ready in $BLIP_ENV_DIR"
echo "  PYTHON=$VENV/bin/python"
echo "  R_LIBS_USER=$RLIB   (susieR built from $SUSIER_ROOT)"
