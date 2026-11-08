#!/usr/bin/env bash
# run_all.sh: run the full workflow in order.
#
# Usage:   bash run_all.sh [threads]
# Options (environment variables):
#   FROM=3                 start at step N (1-6) after fixing a problem
#   EXCLUDE="ID1 ID2"      isolates to omit from step 5
#   OUTGROUP="ID"          outgroup used to root the tree in step 5
#   CHECKM2DB, BAKTA_DB    database paths (see README)
#
# The script stops at the first error. Pause after steps 2 and 3 to review
# QC outputs before trusting later results (see docs/decisions.md).

set -euo pipefail

THREADS="${1:-4}"
FROM="${FROM:-1}"

run() {
  local n="$1"; shift
  if [[ "$n" -lt "$FROM" ]]; then echo "== skipping step $n =="; return; fi
  echo "== step $n: $* =="
  "$@"
}

mkdir -p logs

# each step runs inside its own conda environment (see envs/ and create_envs.sh)
E() { local env="$1"; shift; conda run --no-capture-output -n "$env" "$@"; }

run 1 E chol-qc    bash scripts/01_download_data.sh config/metadata.tsv "$THREADS"
run 2 E chol-qc    bash scripts/02_qc.sh "$THREADS"
run 3 E chol-asm   bash scripts/03_assembly.sh "$THREADS" "${RAM:-8}"
run 4 E chol-amr   bash scripts/04_amr_annotation.sh "$THREADS"
run 5 E chol-phylo bash scripts/05_variants_phylogeny.sh "$THREADS"
run 6 E chol-r     Rscript R/06_visualise.R

echo
echo "Finished. Review results/qc/multiqc, results/assembly_qc/assembly_qc_summary.tsv"
echo "and results/figures/, then record decisions in docs/decisions.md and render the report:"
echo "  quarto render report/report.qmd"
