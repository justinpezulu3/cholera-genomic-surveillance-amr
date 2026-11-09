#!/usr/bin/env bash
# 07_sensitivity.sh
# Sensitivity analysis: repeat the core alignment, Gubbins masking and IQ-TREE steps
# WITHOUT the pre-assembled genomes (and the QC-excluded isolate), to check that the
# main conclusions do not depend on them.
#
# Usage (from the project root):   bash scripts/07_sensitivity.sh [threads]
# Optional: EXCLUDE_IDS="ID1 ID2 ..." overrides the default exclusion list.
# Output:   results/sensitivity_no_assemblies/   (same layout as results/phylogeny/)
#
# How it works: script 05 writes to results/phylogeny/, so the main results are
# backed up first, script 05 is run with EXCLUDE set (all Snippy results are reused),
# the new output is moved aside, and the main results are restored (also on failure).
set -euo pipefail
cd "$(dirname "$0")/.."

EXCL="${EXCLUDE_IDS:-TZA_2017_DarEsSalaam ZMB_2009_Lus38 ZMB_2016_Lus1033 ZMB_2018_Lus7}"
OUT=results/sensitivity_no_assemblies
THREADS="${1:-2}"

[ -d results/phylogeny ] || { echo "results/phylogeny not found: run script 05 first" >&2; exit 1; }
BK=$(mktemp -d)
cp -a results/phylogeny "$BK/phylogeny"
cp results/variants/variant_summary.tsv results/variants/watchlist_variants.tsv "$BK/"
[ -f logs/phylo_tool_versions.txt ] && cp logs/phylo_tool_versions.txt "$BK/"

restore() {
  rm -rf results/phylogeny
  cp -a "$BK/phylogeny" results/phylogeny
  cp "$BK/variant_summary.tsv" "$BK/watchlist_variants.tsv" results/variants/
  [ -f "$BK/phylo_tool_versions.txt" ] && cp "$BK/phylo_tool_versions.txt" logs/ || true
  rm -rf "$BK"
}
trap restore EXIT

PATH="$PATH:$HOME/miniconda3/envs/chol-phylo/bin" EXCLUDE="$EXCL" \
  conda run --no-capture-output -n chol-snippy bash scripts/05_variants_phylogeny.sh "$THREADS" \
  > logs/sens_console.log 2>&1

rm -rf "$OUT"
mv results/phylogeny "$OUT"
echo "Sensitivity results written to $OUT (excluded: $EXCL)"
