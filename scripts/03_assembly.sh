#!/usr/bin/env bash
# 03_assembly.sh
# De novo assembly (SPAdes --isolate on reads subsampled to ~100x), assembly QC (QUAST, CheckM2), species/identity
# checks (fastANI vs reference, MLST) and a pass/flag table.
#
# Usage:   bash scripts/03_assembly.sh [threads] [ram_gb]
# Input:   data/trimmed/<isolate>_{1,2}.fastq.gz  (from 02_qc.sh)
#          data/reference/*/ncbi_dataset/data/*/*.fna  (from 01_download_data.sh)
# Output:  results/assembly/<isolate>/           SPAdes output (spades/contigs.fasta etc.)
#          results/assembly/fasta/<isolate>.fasta  renamed final contigs
#          results/assembly_qc/quast/            QUAST report
#          results/assembly_qc/checkm2/          CheckM2 report
#          results/assembly_qc/fastani.tsv       ANI vs reference
#          results/assembly_qc/mlst.tsv          MLST calls
#          results/assembly_qc/assembly_qc_summary.tsv   metrics + PASS/FLAG
#          logs/assembly.log, logs/assembly_tool_versions.txt
#
# One-time setup for CheckM2 (needs internet, ~3 GB):
#   checkm2 database --download --path ~/databases/checkm2
#   export CHECKM2DB=~/databases/checkm2/CheckM2_database/uniref100.KO.1.dmnd

set -euo pipefail

THREADS="${1:-4}"
# SPAdes memory cap in GB (-m). It is a safety limit, not a requirement: a 4 Mb genome at 100x
# normally needs 2-3 GB. 8 is fine on a 5 GB WSL machine because swap backs it.
# (Shovill was dropped: it refuses --ram < 8 AND --ram > physical RAM, so it cannot run on this laptop.)
RAM="${2:-8}"
TRIM="data/trimmed"
ASM="results/assembly"
FASTA="${ASM}/fasta"
QCDIR="results/assembly_qc"
LOG="logs/assembly.log"

# Assembly parameters (record any change in docs/decisions.md)
MINLEN=200        # drop contigs shorter than this
MINCOV=5          # drop contigs with lower coverage than this
DEPTH=100         # subsample reads to this depth: lowers memory and time, 100x is plenty for a 4 Mb genome
GENOME_SIZE=4030000   # approx. V. cholerae genome size (N16961), used only to compute subsampling depth
SEED=100          # fixed seed so subsampling is reproducible

mkdir -p "$ASM" "$FASTA" "$QCDIR"/{quast,checkm2} logs
log() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }

# Show WHY a step failed (tool output goes to the log file, so it would otherwise be hidden)
trap 'rc=$?; echo "ERROR: 03_assembly.sh failed (exit $rc) near line $LINENO. Last lines of $LOG:" >&2; tail -n 30 "$LOG" >&2' ERR

{
  echo "# Assembly tool versions ($(date '+%F %T'))"
  spades.py --version 2>&1 | head -n 2
  seqtk 2>&1 | sed -n 2p
  quast.py --version 2>&1 | head -n 1
  checkm2 --version 2>&1 | head -n 1
  fastANI --version 2>&1 | head -n 1
  mlst --version 2>&1
} > logs/assembly_tool_versions.txt 2>&1 || true

REF=$(ls data/reference/*/ncbi_dataset/data/*/*.fna 2>/dev/null | head -n 1 || true)
[[ -z "$REF" ]] && { echo "ERROR: reference .fna not found under data/reference/. Run 01_download_data.sh." >&2; exit 1; }
log "Reference: $REF"

shopt -s nullglob
R1S=("$TRIM"/*_1.fastq.gz)
[[ ${#R1S[@]} -eq 0 ]] && { echo "ERROR: no trimmed reads in $TRIM. Run 02_qc.sh first." >&2; exit 1; }

# --- 1. Assembly ---------------------------------------------------------------
for R1 in "${R1S[@]}"; do
  ID=$(basename "$R1" _1.fastq.gz)
  R2="${TRIM}/${ID}_2.fastq.gz"

  if [[ -s "${FASTA}/${ID}.fasta" ]]; then
    log "SKIP $ID (assembly exists)"; continue
  fi

  OUT="${ASM}/${ID}"
  rm -rf "$OUT"; mkdir -p "$OUT"

  # (a) Subsample to ~DEPTH x using the base count from fastp (reads were already trimmed in 02_qc.sh)
  BASES=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["summary"]["after_filtering"]["total_bases"])' "results/qc/fastp/${ID}.json")
  FRAC=$(awk -v b="$BASES" -v d="$DEPTH" -v g="$GENOME_SIZE" 'BEGIN{f=d*g/b; if(f>1)f=1; printf "%.4f", f}')
  log "Assembling $ID: ${BASES} bases after trimming, subsampling fraction ${FRAC} (target ${DEPTH}x)"
  if [[ "$FRAC" == "1.0000" ]]; then
    S1="$R1"; S2="$R2"
  else
    S1="${OUT}/sub_1.fastq.gz"; S2="${OUT}/sub_2.fastq.gz"
    seqtk sample -s "$SEED" "$R1" "$FRAC" | gzip > "$S1"   # same seed + fraction on both files keeps pairs in sync
    seqtk sample -s "$SEED" "$R2" "$FRAC" | gzip > "$S2"
  fi

  # (b) SPAdes in isolate mode (high-coverage single-colony Illumina data)
  spades.py --isolate -1 "$S1" -2 "$S2" -o "${OUT}/spades" -t "$THREADS" -m "$RAM" >> "$LOG" 2>&1

  # (c) Keep contigs >= MINLEN bp and >= MINCOV x (SPAdes headers: NODE_n_length_L_cov_C); rename tidily
  python3 - "$ID" "${OUT}/spades/contigs.fasta" "${FASTA}/${ID}.fasta" "$MINLEN" "$MINCOV" <<'PY' | tee -a "$LOG"
import sys, re
iso, src, dst, minlen, mincov = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), float(sys.argv[5])
def records(path):
    h, s = None, []
    for line in open(path):
        if line.startswith(">"):
            if h is not None: yield h, "".join(s)
            h, s = line[1:].strip(), []
        else:
            s.append(line.strip())
    if h is not None: yield h, "".join(s)
kept = dropped = 0
with open(dst, "w") as out:
    for h, seq in records(src):
        m = re.search(r"length_(\d+)_cov_([\d.]+)", h)
        length, cov = (int(m.group(1)), float(m.group(2))) if m else (len(seq), 0.0)
        if length < minlen or cov < mincov:
            dropped += 1; continue
        kept += 1
        out.write(f">{iso}_c{kept:03d} len={length} cov={cov:.1f}\n")
        for i in range(0, len(seq), 80): out.write(seq[i:i+80] + "\n")
print(f"{iso}: kept {kept} contigs, dropped {dropped} (<{minlen} bp or <{mincov}x)")
PY

  # (d) Free disk: remove subsampled reads and bulky SPAdes intermediates (keeps contigs, scaffolds, graph, log)
  rm -f "${OUT}/sub_1.fastq.gz" "${OUT}/sub_2.fastq.gz"
  rm -rf "${OUT}/spades/corrected" "${OUT}/spades/misc" "${OUT}/spades/tmp" "${OUT}"/spades/K*
done

# --- 1b. Pre-assembled genomes (data_type=assembly) join the same QC set -----------
for A in data/assemblies/*.fasta; do
  [[ -e "$A" ]] || continue
  AID=$(basename "$A" .fasta)
  if [[ ! -s "${FASTA}/${AID}.fasta" ]]; then
    cp "$A" "${FASTA}/${AID}.fasta"; log "Added pre-assembled genome $AID (not assembled here)"
  fi
done

FASTAS=("$FASTA"/*.fasta)
IDS=(); for f in "${FASTAS[@]}"; do IDS+=("$(basename "$f" .fasta)"); done

# --- 2. QUAST ------------------------------------------------------------------
log "QUAST"
quast.py "${FASTAS[@]}" -r "$REF" -l "$(IFS=,; echo "${IDS[*]}")" \
  -t "$THREADS" --min-contig "$MINLEN" -o "$QCDIR/quast" --silent >> "$LOG" 2>&1

# --- 3. CheckM2 ----------------------------------------------------------------
log "CheckM2"
# CheckM2 normally lives in its own environment (chol-checkm2); use it if not on PATH
CHECKM2=()
if command -v checkm2 >/dev/null 2>&1; then CHECKM2=(checkm2)
elif command -v conda >/dev/null 2>&1 && conda env list | awk '{print $1}' | grep -qx chol-checkm2; then
  CHECKM2=(conda run --no-capture-output -n chol-checkm2 checkm2)
fi
if [[ ${#CHECKM2[@]} -gt 0 && -n "${CHECKM2DB:-}" ]]; then
  # --lowmem reduces RAM use (slower); appropriate for a ~5 GB machine
  "${CHECKM2[@]}" predict --threads "$THREADS" --input "$FASTA" --extension fasta \
    --output-directory "$QCDIR/checkm2" --database_path "$CHECKM2DB" --lowmem --force >> "$LOG" 2>&1
else
  log "WARN: CheckM2 or CHECKM2DB not available; completeness/contamination will be NA"
fi

# --- 4. Species/identity check: ANI vs reference -------------------------------
log "fastANI vs reference"
printf '%s\n' "${FASTAS[@]}" > "$QCDIR/query_list.txt"
fastANI --ql "$QCDIR/query_list.txt" -r "$REF" -t "$THREADS" -o "$QCDIR/fastani.tsv" >> "$LOG" 2>&1 || true

# --- 5. MLST (V. cholerae scheme; limited resolution for 7th-pandemic isolates) -
log "MLST"
mlst --scheme vcholerae "${FASTAS[@]}" > "$QCDIR/mlst.tsv" 2>> "$LOG" || \
  log "WARN: MLST failed; check 'mlst --list' for the V. cholerae scheme name"

# --- 6. Summary table with transparent thresholds ------------------------------
python3 - <<'PY'
import csv, os

# Thresholds: adjust only with justification, and report them in the Methods.
SIZE_MIN, SIZE_MAX = 3.8e6, 4.4e6   # expected V. cholerae genome ~4.0 Mb
MAX_CONTIGS        = 300
MIN_N50            = 50_000
MIN_COMPLETENESS   = 95.0
MAX_CONTAMINATION  = 5.0
MIN_ANI            = 95.0

def read_tsv(path):
    if not os.path.exists(path): return []
    with open(path) as f: return list(csv.DictReader(f, delimiter="\t"))

quast = {r["Assembly"]: r for r in read_tsv("results/assembly_qc/quast/transposed_report.tsv")}
checkm = {r["Name"]: r for r in read_tsv("results/assembly_qc/checkm2/quality_report.tsv")}

ani = {}
if os.path.exists("results/assembly_qc/fastani.tsv"):
    for line in open("results/assembly_qc/fastani.tsv"):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 3:
            ani[os.path.basename(p[0]).replace(".fasta", "")] = float(p[2])

mlst = {}
if os.path.exists("results/assembly_qc/mlst.tsv"):
    for line in open("results/assembly_qc/mlst.tsv"):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 3:
            mlst[os.path.basename(p[0]).replace(".fasta", "")] = f"{p[1]}:ST{p[2]}"

def num(x):
    try: return float(x)
    except (TypeError, ValueError): return None

hdr = ["isolate", "total_length", "contigs", "N50", "GC_pct", "completeness",
       "contamination", "ANI_vs_ref", "MLST", "status", "flags"]
rows = []
for iso, q in quast.items():
    length = num(q.get("Total length"))
    contigs = num(q.get("# contigs"))
    n50 = num(q.get("N50"))
    gc = q.get("GC (%)", "NA")
    c = checkm.get(iso, {})
    comp, cont = num(c.get("Completeness")), num(c.get("Contamination"))
    a = ani.get(iso)

    flags = []
    if length is None or not (SIZE_MIN <= length <= SIZE_MAX): flags.append("genome_size")
    if contigs is None or contigs > MAX_CONTIGS: flags.append("fragmented")
    if n50 is None or n50 < MIN_N50: flags.append("low_N50")
    if comp is None: flags.append("completeness_NA")
    elif comp < MIN_COMPLETENESS: flags.append("low_completeness")
    if cont is None: flags.append("contamination_NA")
    elif cont > MAX_CONTAMINATION: flags.append("contamination")
    if a is None or a < MIN_ANI: flags.append("ANI_below_species_threshold")

    rows.append([iso, length, contigs, n50, gc, comp, cont, a, mlst.get(iso, "NA"),
                 "PASS" if not flags else "FLAG", ";".join(flags) or "-"])

with open("results/assembly_qc/assembly_qc_summary.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t"); w.writerow(hdr); w.writerows(rows)

print("Wrote results/assembly_qc/assembly_qc_summary.tsv")
for r in rows: print("\t".join(map(str, r)))
PY

log "Assembly complete. Review assembly_qc_summary.tsv: do not carry FLAG isolates forward without a documented decision."
