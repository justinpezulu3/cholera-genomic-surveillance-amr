#!/usr/bin/env bash
# 04_amr_annotation.sh
# Genome annotation (Bakta), AMR determinant detection (AMRFinderPlus + ABRicate),
# targeted virulence/lineage marker lookup, and merged determinant tables.
#
# Usage:   bash scripts/04_amr_annotation.sh [threads]
# Input:   results/assembly/fasta/<isolate>.fasta
#          results/assembly_qc/assembly_qc_summary.tsv  (optional; FLAG isolates are noted)
# Output:  results/annotation/<isolate>/         Bakta output
#          results/amr/amrfinder/<isolate>.tsv
#          results/amr/abricate/<isolate>.<db>.tsv
#          results/amr/targets/<isolate>.markers.tsv
#          results/amr/amr_long.tsv              every hit, every tool, with identity/coverage
#          results/amr/amr_matrix.tsv            isolate x determinant presence, with tool support
#          logs/amr.log, logs/amr_tool_versions.txt
#
# One-time database setup (needs internet):
#   bakta_db download --output ~/databases/bakta --type light
#   export BAKTA_DB=~/databases/bakta/db-light
#   amrfinder -u                      # updates AMRFinderPlus database
#   abricate-get_db --db resfinder --force ; abricate-get_db --db card --force
#
# IMPORTANT: a detected determinant is evidence of resistance-ASSOCIATED genotype,
# not proof of phenotypic resistance. Do not label isolates "resistant" from this output.

set -euo pipefail

THREADS="${1:-4}"
FASTA="results/assembly/fasta"
ANN="results/annotation"
AMR="results/amr"
LOG="logs/amr.log"

# ABRicate thresholds (report in Methods)
MINID=90
MINCOV=80

mkdir -p "$ANN" "$AMR"/{amrfinder,abricate,targets} logs
log() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }

# --- versions and database provenance -----------------------------------------
{
  echo "# AMR tool versions ($(date '+%F %T'))"
  bakta --version 2>&1
  amrfinder --version 2>&1
  amrfinder --database_version 2>&1 || true
  abricate --version 2>&1
  abricate --list 2>&1 || true
} > logs/amr_tool_versions.txt 2>&1 || true

shopt -s nullglob
FASTAS=("$FASTA"/*.fasta)
[[ ${#FASTAS[@]} -eq 0 ]] && { echo "ERROR: no assemblies in $FASTA. Run 03_assembly.sh first." >&2; exit 1; }

# Decide whether AMRFinderPlus supports the organism option (needed for point mutations)
AF_ORG=()
if amrfinder --list_organisms 2>/dev/null | grep -q "Vibrio_cholerae"; then
  AF_ORG=(--organism Vibrio_cholerae)
  log "AMRFinderPlus: using --organism Vibrio_cholerae (point mutations enabled)"
else
  log "WARN: Vibrio_cholerae not in 'amrfinder --list_organisms'; update with 'amrfinder -u'. Point mutations (e.g. gyrA/parC) will NOT be screened."
fi

for F in "${FASTAS[@]}"; do
  ID=$(basename "$F" .fasta)

  # ---- 1. Bakta annotation ----
  if [[ -s "${ANN}/${ID}/${ID}.tsv" ]]; then
    log "SKIP Bakta $ID"
  elif [[ -n "${BAKTA_DB:-}" ]]; then
    log "Bakta: $ID"
    bakta --db "$BAKTA_DB" --output "${ANN}/${ID}" --prefix "$ID" \
      --genus Vibrio --species cholerae --locus-tag "${ID//[^A-Za-z0-9]/}" \
      --threads "$THREADS" --force "$F" >> "$LOG" 2>&1
  else
    log "WARN: BAKTA_DB not set; skipping annotation for $ID (marker lookup will be limited)"
  fi

  # ---- 2. AMRFinderPlus (nucleotide mode) ----
  if [[ ! -s "${AMR}/amrfinder/${ID}.tsv" ]]; then
    log "AMRFinderPlus: $ID"
    amrfinder -n "$F" "${AF_ORG[@]}" --plus --name "$ID" \
      --threads "$THREADS" -o "${AMR}/amrfinder/${ID}.tsv" >> "$LOG" 2>&1
  fi

  # ---- 3. ABRicate (independent second opinion) ----
  for DB in resfinder card; do
    if abricate --list 2>/dev/null | awk '{print $1}' | grep -qx "$DB"; then
      [[ -s "${AMR}/abricate/${ID}.${DB}.tsv" ]] || \
        abricate --db "$DB" --minid "$MINID" --mincov "$MINCOV" --threads "$THREADS" \
          "$F" > "${AMR}/abricate/${ID}.${DB}.tsv" 2>> "$LOG"
    else
      log "WARN: ABRicate database '$DB' not installed; skipping"
    fi
  done

  # ---- 4. Marker lookup in Bakta annotation (supporting evidence only) ----
  # Gene names come from the annotation and can be missed or misnamed; confirm
  # important calls by BLAST/alignment before stating them in the report.
  OUT="${AMR}/targets/${ID}.markers.tsv"
  echo -e "marker\tmatches_in_annotation" > "$OUT"
  BTSV="${ANN}/${ID}/${ID}.tsv"
  for M in ctxA ctxB tcpA toxR toxT; do
    if [[ -s "$BTSV" ]]; then
      N=$(awk -F'\t' -v m="$M" 'BEGIN{IGNORECASE=1} !/^#/ && ($4==m){c++} END{print c+0}' "$BTSV")
    else N="NA"; fi
    echo -e "${M}\t${N}" >> "$OUT"
  done
  if [[ -s "$BTSV" ]]; then
    N=$(awk -F'\t' 'BEGIN{IGNORECASE=1} !/^#/ && ($7 ~ /SXT/ || $4 ~ /^int/){c++} END{print c+0}' "$BTSV")
  else N="NA"; fi
  echo -e "SXT_related_annotations\t${N}" >> "$OUT"

  # ---- 5. Optional: SXT/ICE coverage by BLAST if you supply a reference ----
  # Place an SXT/ICE reference (e.g. ICEVchInd5; choose and cite the accession) at
  # data/targets/sxt_ice.fasta to enable this.
  if [[ -s data/targets/sxt_ice.fasta ]] && command -v blastn >/dev/null 2>&1; then
    blastn -query data/targets/sxt_ice.fasta -subject "$F" \
      -outfmt "6 qseqid sseqid pident length qlen qstart qend" \
      > "${AMR}/targets/${ID}.sxt_blast.tsv" 2>> "$LOG" || true
  fi
done

# --- 6. Merge into long and wide tables ---------------------------------------
python3 - <<'PY'
import csv, glob, os, collections

def rows(path, delim="\t"):
    with open(path, newline="") as f:
        header = f.readline().lstrip("#").rstrip("\n").split(delim)
        for line in f:
            if line.startswith("#") or not line.strip(): continue
            yield dict(zip(header, line.rstrip("\n").split(delim)))

def pick(d, *names, default="NA"):
    for n in names:
        if n in d and d[n] != "": return d[n]
    return default

long_rows = []

# AMRFinderPlus (column names differ between versions, so try both)
for p in sorted(glob.glob("results/amr/amrfinder/*.tsv")):
    iso = os.path.basename(p)[:-4]
    for d in rows(p):
        long_rows.append([
            iso, "AMRFinderPlus",
            pick(d, "Element symbol", "Gene symbol"),
            pick(d, "Element name", "Sequence name"),
            pick(d, "Type", "Element type"),
            pick(d, "Subtype", "Element subtype"),
            pick(d, "Class"), pick(d, "Subclass"),
            pick(d, "Method"),
            pick(d, "% Identity to reference sequence", "% Identity to reference"),
            pick(d, "% Coverage of reference sequence", "% Coverage of reference"),
            pick(d, "Contig id"), pick(d, "Start"), pick(d, "Stop"),
        ])

# ABRicate
for p in sorted(glob.glob("results/amr/abricate/*.tsv")):
    base = os.path.basename(p)[:-4]
    iso, db = base.rsplit(".", 1)
    for d in rows(p):
        long_rows.append([
            iso, f"ABRicate:{db}",
            pick(d, "GENE"), pick(d, "PRODUCT"), "NA", "NA",
            pick(d, "RESISTANCE"), "NA", "NA",
            pick(d, "%IDENTITY"), pick(d, "%COVERAGE"),
            pick(d, "SEQUENCE"), pick(d, "START"), pick(d, "END"),
        ])

hdr = ["isolate", "tool", "symbol", "name", "type", "subtype", "class", "subclass",
       "method", "pct_identity", "pct_coverage", "contig", "start", "stop"]
with open("results/amr/amr_long.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t"); w.writerow(hdr); w.writerows(long_rows)

# Wide matrix: rows = isolate, columns = determinant symbol; value = tools supporting it.
# Matching is by case-insensitive symbol, so naming differences between tools can
# hide agreement (e.g. allele-level vs gene-level names): review manually.
isolates = sorted({r[0] for r in long_rows} | {os.path.basename(p)[:-4] for p in glob.glob("results/amr/amrfinder/*.tsv")})
support = collections.defaultdict(set)   # (iso, symbol) -> set(tools)
for r in long_rows:
    iso, tool, sym = r[0], r[1], r[2]
    if sym in ("NA", ""): continue
    support[(iso, sym.lower())].add("AF" if tool == "AMRFinderPlus" else "AB")
syms = sorted({s for (_, s) in support})
with open("results/amr/amr_matrix.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t"); w.writerow(["isolate"] + syms)
    for iso in isolates:
        w.writerow([iso] + ["+".join(sorted(support.get((iso, s), []))) or "-" for s in syms])

print(f"Wrote results/amr/amr_long.tsv ({len(long_rows)} hits) and results/amr/amr_matrix.tsv "
      f"({len(isolates)} isolates x {len(syms)} determinants)")
print("Matrix legend: AF = AMRFinderPlus, AB = ABRicate, AF+AB = both, - = not detected")
PY

log "Done. Next: inspect amr_long.tsv for identity/coverage and partial hits before plotting."
log "Remember: determinants are resistance-associated genotype, not confirmed phenotype."
