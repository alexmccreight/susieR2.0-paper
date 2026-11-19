#!/bin/sh
set -eu

if [ "$#" -ne 1 ]; then
  printf 'Usage: sh reproduce.sh OUTPUT_DIRECTORY\n' >&2
  exit 2
fi

src=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
mkdir -p -- "$1"
out=$(CDPATH= cd -- "$1" && pwd)
work=$(mktemp -d "$out/.s7e.XXXXXX")
trap 'rm -rf "$work"' EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

sparsepro=e13ac87d6e650d7682ff1f33f6f9898d66932d64
susier=ed6e491e3d0235ea5fa2e8545ea206b5d759eb6e
fetch() {
  curl -fsSL --retry 2 "$1" -o "$work/$2"
  printf '%s  %s\n' "$3" "$work/$2" | shasum -a 256 -c
}

fetch "https://raw.githubusercontent.com/zhwm/SparsePro/$sparsepro/sparsepro_zld.py" \
  sparsepro_zld.py b4db3c8c9e0dae199bfc8c1ef7fe3e2820c6f6f5311942c00aac2cb011850319
fetch "https://raw.githubusercontent.com/stephenslab/susieR/$susier/data/rss_mismatch_cases.RData" \
  rss_mismatch_cases.RData b3bb503f4b28feb578048422e3758ef8f9fda7a30a7c7da10e416830d8faa483
fetch "https://raw.githubusercontent.com/stephenslab/susieR/$susier/R/susie_get_functions.R" \
  susie_get_functions.R d1f7bcd0dc03d5e9b5cd0c0555dc9c9950d059b370870f6683c454518954ba7b

Rscript --vanilla "$src/simulate.R" "$work"
python3 -W ignore::SyntaxWarning "$src/compare_reporting.py" \
  "$work" "$src/../data/paired-posterior-mass.csv"
cp "$work/paired-posterior-mass.csv" "$out/"
cd "$src/../../.."
Rscript --vanilla Supplementary_Figures/S4/plot_S4.R \
  "$out/paired-posterior-mass.csv" "$out"
printf 'Wrote S4.pdf (manuscript S7), S4.png, and paired-posterior-mass.csv to %s\n' "$out"
