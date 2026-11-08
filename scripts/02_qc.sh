#!/usr/bin/env bash
# 02_qc.sh
# Read quality control: FastQC (raw) -> fastp trimming -> FastQC (trimmed) -> MultiQC.
#
# Usage:   bash scripts/02_qc.sh [threads]
# Input:   data/raw/<isolate>_{1,2}.fastq.gz
# Output:  results/qc/raw_fastqc/         FastQC on raw reads
#          data/trimmed/<isolate>_{1,2}.fastq.gz
#          results/qc/fastp/<isolate>.{json,html}
#          results/qc/trimmed_fastqc/     FastQC on trimmed reads
#          results/qc/multiqc/            combined report
#          results/qc/read_summary.tsv    pairs before/after, Q30, GC
#          logs/qc.log, logs/qc_tool_versions.txt
#
# Trimming parameters are deliberately modest and are recorded in the log.
# Adjust them only after inspecting the raw FastQC/MultiQC report.

set -euo pipefail

THREADS="${1:-4}"
RAW="data/raw"
TRIM="data/trimmed"
QC="results/qc"
LOG="logs/qc.log"

# fastp parameters (document any change in docs/decisions.md)
MIN_Q=20          # sliding-window / per-base quality cutoff
MIN_LEN=50        # discard reads shorter than this after trimming
WINDOW=4          # sliding window size

mkdir -p "$TRIM" "$QC"/{raw_fastqc,trimmed_fastqc,fastp,multiqc} logs
log() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }

# Show WHY a step failed instead of exiting silently (tool output goes to the log file)
trap 'rc=$?; echo "ERROR: 02_qc.sh failed (exit $rc) near line $LINENO. Last lines of $LOG:" >&2; tail -n 25 "$LOG" >&2' ERR

{
  echo "# QC tool versions ($(date '+%F %T'))"
  fastqc --version
  fastp --version 2>&1
  multiqc --version
} > logs/qc_tool_versions.txt 2>&1

shopt -s nullglob
R1S=("$RAW"/*_1.fastq.gz)
[[ ${#R1S[@]} -eq 0 ]] && { echo "ERROR: no FASTQ files in $RAW. Run 01_download_data.sh first." >&2; exit 1; }

# --- 0. Integrity check: a download interrupted mid-compression leaves a truncated
#        .fastq.gz that looks present but crashes FastQC. Find these before QC. -------
BAD=()
if [[ "${SKIP_GZ_CHECK:-0}" == "1" ]]; then
  log "Skipping gzip integrity check (SKIP_GZ_CHECK=1)"
else
  log "Checking gzip integrity of raw FASTQ files (a few minutes; skip on reruns with SKIP_GZ_CHECK=1)"
  for f in "$RAW"/*.fastq.gz; do
    if ! gzip -t "$f" 2>>"$LOG"; then
      log "CORRUPT or TRUNCATED: $f"; BAD+=("$f")
    fi
  done
fi
for f in "$RAW"/*.fastq; do
  [[ -e "$f" ]] && log "WARNING: leftover uncompressed file $f (interrupted download?)"
done
if [[ ${#BAD[@]} -gt 0 ]]; then
  echo "ERROR: ${#BAD[@]} damaged file(s). Delete them and rerun 01_download_data.sh:" >&2
  printf '  rm %s\n' "${BAD[@]}" >&2
  exit 1
fi

# --- 1. FastQC on raw reads ----------------------------------------------------
log "FastQC on raw reads"
fastqc -t "$THREADS" --outdir "$QC/raw_fastqc" "$RAW"/*.fastq.gz >> "$LOG" 2>&1

# --- 2. fastp trimming per isolate --------------------------------------------
for R1 in "${R1S[@]}"; do
  ID=$(basename "$R1" _1.fastq.gz)
  R2="${RAW}/${ID}_2.fastq.gz"
  [[ ! -s "$R2" ]] && { log "SKIP $ID: missing R2"; continue; }

  if [[ -s "${TRIM}/${ID}_1.fastq.gz" && -s "${QC}/fastp/${ID}.json" ]]; then
    log "SKIP $ID (already trimmed)"; continue
  fi

  log "fastp: $ID"
  fastp \
    --in1 "$R1" --in2 "$R2" \
    --out1 "${TRIM}/${ID}_1.fastq.gz" --out2 "${TRIM}/${ID}_2.fastq.gz" \
    --detect_adapter_for_pe \
    --cut_right --cut_right_window_size "$WINDOW" --cut_right_mean_quality "$MIN_Q" \
    --length_required "$MIN_LEN" \
    --thread "$THREADS" \
    --json "${QC}/fastp/${ID}.json" --html "${QC}/fastp/${ID}.html" \
    --report_title "$ID" >> "$LOG" 2>&1
done

# --- 3. FastQC on trimmed reads ------------------------------------------------
log "FastQC on trimmed reads"
fastqc -t "$THREADS" --outdir "$QC/trimmed_fastqc" "$TRIM"/*.fastq.gz >> "$LOG" 2>&1

# --- 4. MultiQC ----------------------------------------------------------------
log "MultiQC"
multiqc "$QC" -o "$QC/multiqc" --force >> "$LOG" 2>&1

# --- 5. Per-isolate summary table from fastp JSON -----------------------------
python3 - <<'PY'
import json, glob, os
rows = []
for j in sorted(glob.glob("results/qc/fastp/*.json")):
    d = json.load(open(j))
    iso = os.path.basename(j).replace(".json", "")
    b, a = d["summary"]["before_filtering"], d["summary"]["after_filtering"]
    rows.append([
        iso,
        b["total_reads"] // 2, a["total_reads"] // 2,
        round(b["q30_rate"] * 100, 2), round(a["q30_rate"] * 100, 2),
        round(a["gc_content"] * 100, 2),
        round(d["duplication"]["rate"] * 100, 2),
    ])
hdr = ["isolate", "pairs_before", "pairs_after", "q30_before_pct",
       "q30_after_pct", "gc_after_pct", "duplication_pct"]
with open("results/qc/read_summary.tsv", "w") as f:
    f.write("\t".join(hdr) + "\n")
    for r in rows:
        f.write("\t".join(map(str, r)) + "\n")
print(f"Wrote results/qc/read_summary.tsv ({len(rows)} isolates)")
PY

log "QC complete. Open results/qc/multiqc/multiqc_report.html and check read_summary.tsv."
log "Review before assembly: low Q30, very low read counts, GC far from ~47-48%, or adapter content."
