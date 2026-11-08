#!/usr/bin/env bash
# setup_tools_and_databases.sh
# One-off setup: (1) verify that every tool is installed in its environment,
# (2) download the databases (CheckM2, Bakta light, AMRFinderPlus, ABRicate),
# (3) save the database paths in ~/.bashrc so the pipeline scripts can find them.
#
# Usage:   bash scripts/setup_tools_and_databases.sh [database_dir]      (default: ~/databases)
# Safe to rerun: finished downloads are skipped. Do this on a good connection (a few GB).
# Run from the project root, in a normal Ubuntu terminal (not inside another conda env).

set -uo pipefail

DB="${1:-$HOME/databases}"
mkdir -p "$DB"

echo "=== 1. Tool check ==="
MISSING=0
check() {
  local env="$1"; shift
  for t in "$@"; do
    if conda run -n "$env" which "$t" >/dev/null 2>&1; then
      printf "  %-14s %-18s OK\n" "$env" "$t"
    else
      printf "  %-14s %-18s MISSING\n" "$env" "$t"; MISSING=$((MISSING+1))
    fi
  done
}
check chol-qc     prefetch fasterq-dump datasets fastqc fastp multiqc esearch
check chol-asm    shovill spades.py quast.py fastANI mlst
check chol-checkm2 checkm2
check chol-amr    bakta bakta_db amrfinder abricate blastn
check chol-phylo  snippy snippy-core run_gubbins.py snp-sites snp-dists
check chol-r      R Rscript quarto
# iqtree may be installed as iqtree2 or iqtree
if conda run -n chol-phylo which iqtree2 >/dev/null 2>&1 || conda run -n chol-phylo which iqtree >/dev/null 2>&1; then
  echo "  chol-phylo     iqtree/iqtree2     OK"
else
  echo "  chol-phylo     iqtree/iqtree2     MISSING"; MISSING=$((MISSING+1))
fi
echo "Missing tools: $MISSING"
[[ "$MISSING" -gt 0 ]] && echo "Send me the MISSING lines before continuing."

echo
echo "=== 2. Databases ==="

# CheckM2 (about 3 GB)
CM2=$(find "$DB/checkm2" -name "*.dmnd" 2>/dev/null | head -n 1 || true)
if [[ -z "$CM2" ]]; then
  echo "Downloading CheckM2 database..."
  conda run --no-capture-output -n chol-checkm2 checkm2 database --download --path "$DB/checkm2"
  CM2=$(find "$DB/checkm2" -name "*.dmnd" 2>/dev/null | head -n 1 || true)
else
  echo "CheckM2 database already present: $CM2"
fi

# Bakta light database
BK="$DB/bakta/db-light"
if [[ ! -d "$BK" ]]; then
  echo "Downloading Bakta light database..."
  conda run --no-capture-output -n chol-amr bakta_db download --output "$DB/bakta" --type light
else
  echo "Bakta database already present: $BK"
fi

# AMRFinderPlus database (updates itself)
echo "Updating AMRFinderPlus database..."
conda run --no-capture-output -n chol-amr amrfinder -u

# ABRicate databases
for d in resfinder card; do
  if conda run -n chol-amr abricate --list 2>/dev/null | awk '{print $1}' | grep -qx "$d"; then
    echo "ABRicate database present: $d"
  else
    echo "Fetching ABRicate database: $d"
    conda run --no-capture-output -n chol-amr abricate-get_db --db "$d" --force
  fi
done

echo
echo "=== 3. Saving database paths ==="
add_var() {
  local name="$1" value="$2"
  [[ -z "$value" ]] && { echo "  WARNING: $name path not found; not saved"; return; }
  sed -i "/^export ${name}=/d" "$HOME/.bashrc"
  echo "export ${name}=${value}" >> "$HOME/.bashrc"
  echo "  $name=$value"
}
add_var CHECKM2DB "$CM2"
[[ -d "$BK" ]] && add_var BAKTA_DB "$BK"

echo
echo "Done. Run:  source ~/.bashrc    then check:  echo \$CHECKM2DB \$BAKTA_DB"
echo "Database versions are recorded by scripts 03 and 04; note the download date in docs/decisions.md."
