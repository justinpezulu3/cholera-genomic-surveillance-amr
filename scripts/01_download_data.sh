#!/usr/bin/env bash
# 01_download_data.sh
# Download reads (SRA) and the reference genome listed in config/metadata.tsv.
#
# Usage:   bash scripts/01_download_data.sh [config/metadata.tsv] [threads]
# Needs:   sra-tools (prefetch, fasterq-dump, vdb-validate), pigz or gzip, and
#          NCBI "datasets" CLI (for the reference assembly).
# Output:  data/raw/<isolate_id>_{1,2}.fastq.gz
#          data/reference/<assembly_accession>/
#          logs/download.log, logs/tool_versions.txt, data/raw/checksums.md5

set -euo pipefail

META="${1:-config/metadata.tsv}"
THREADS="${2:-4}"
RAW_DIR="data/raw"
REF_DIR="data/reference"
LOG_DIR="logs"
LOG="${LOG_DIR}/download.log"

mkdir -p "$RAW_DIR" "$REF_DIR" "$LOG_DIR" data/assemblies tmp

log() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }

# --- record tool versions for reproducibility ---------------------------------
{
  echo "# Tool versions ($(date '+%F %T'))"
  prefetch --version 2>&1 | head -n 2 || true
  fasterq-dump --version 2>&1 | head -n 2 || true
  datasets --version 2>&1 | head -n 1 || true
} > "${LOG_DIR}/tool_versions.txt"

# --- locate columns by header name so column order can change -----------------
col() { head -n 1 "$META" | tr '\t' '\n' | grep -nx "$1" | cut -d: -f1; }
C_ID=$(col isolate_id); C_RUN=$(col sra_run); C_ASM=$(col assembly_accession); C_ROLE=$(col role)
C_TYPE=$(col data_type || true)   # 'reads' (default) or 'assembly'

[[ -z "$C_ID" || -z "$C_RUN" || -z "$C_ASM" || -z "$C_ROLE" ]] && {
  echo "ERROR: metadata is missing one of isolate_id, sra_run, assembly_accession, role" >&2; exit 1; }

COMPRESS="gzip -f"
command -v pigz >/dev/null 2>&1 && COMPRESS="pigz -f -p ${THREADS}"

tail -n +2 "$META" | while IFS=$'\t' read -r -a F; do
  ID="${F[$((C_ID-1))]}"; RUN="${F[$((C_RUN-1))]}"
  ASM="${F[$((C_ASM-1))]}"; ROLE="${F[$((C_ROLE-1))]}"

  [[ -z "$ID" || "$ID" == \#* ]] && continue
  TYPE="reads"; [[ -n "$C_TYPE" ]] && TYPE="${F[$((C_TYPE-1))]}"

  # ---- reference genome: download assembly, not reads ----
  if [[ "$ROLE" == "reference" ]]; then
    if [[ -d "${REF_DIR}/${ASM}" ]]; then
      log "SKIP reference ${ASM} (already downloaded)"; continue
    fi
    log "Downloading reference ${ASM}"
    datasets download genome accession "$ASM" \
      --include genome,gff3,gbff --filename "${REF_DIR}/${ASM}.zip" >> "$LOG" 2>&1
    unzip -oq "${REF_DIR}/${ASM}.zip" -d "${REF_DIR}/${ASM}"
    rm -f "${REF_DIR}/${ASM}.zip"
    continue
  fi

  # ---- pre-assembled genomes (data_type = assembly): no reads, enter at QC/AMR/phylogeny ----
  if [[ "$TYPE" == "assembly" ]]; then
    if [[ -s "data/assemblies/${ID}.fasta" ]]; then
      log "SKIP assembly ${ID} (already downloaded)"; continue
    fi
    if [[ -z "$ASM" || "$ASM" == "TODO" || "$ASM" == "NA" ]]; then
      log "SKIP ${ID}: no assembly accession filled in yet"; continue
    fi
    log "Downloading assembly ${ID} (${ASM})"
    datasets download genome accession "$ASM" --include genome \
      --filename "tmp/${ID}.zip" >> "$LOG" 2>&1
    unzip -oq "tmp/${ID}.zip" -d "tmp/${ID}_asm"
    cat tmp/"${ID}"_asm/ncbi_dataset/data/*/*.fna > "data/assemblies/${ID}.fasta"
    rm -rf "tmp/${ID}.zip" "tmp/${ID}_asm"
    continue
  fi

  # ---- study/comparator/outgroup reads ----
  if [[ -z "$RUN" || "$RUN" == "TODO" || "$RUN" == "NA" ]]; then
    log "SKIP ${ID}: no SRA run accession filled in yet"; continue
  fi
  if [[ -s "${RAW_DIR}/${ID}_1.fastq.gz" && -s "${RAW_DIR}/${ID}_2.fastq.gz" ]]; then
    log "SKIP ${ID} (${RUN}): FASTQ already present"; continue
  fi

  log "Downloading ${ID} (${RUN})"
  prefetch "$RUN" --output-directory tmp --max-size 20G >> "$LOG" 2>&1
  vdb-validate "tmp/${RUN}" >> "$LOG" 2>&1 || { log "WARN: validation issue for ${RUN}"; }

  fasterq-dump "tmp/${RUN}" --split-files --threads "$THREADS" \
    --outdir tmp --temp tmp >> "$LOG" 2>&1

  # Paired-end check: fail loudly rather than silently continue with singles
  if [[ ! -s "tmp/${RUN}_1.fastq" || ! -s "tmp/${RUN}_2.fastq" ]]; then
    log "ERROR: ${RUN} did not produce paired FASTQ files; check layout in metadata"
    exit 1
  fi

  mv "tmp/${RUN}_1.fastq" "${RAW_DIR}/${ID}_1.fastq"
  mv "tmp/${RUN}_2.fastq" "${RAW_DIR}/${ID}_2.fastq"
  $COMPRESS "${RAW_DIR}/${ID}_1.fastq" "${RAW_DIR}/${ID}_2.fastq"
  rm -rf "tmp/${RUN}" "tmp/${RUN}"*.fastq 2>/dev/null || true

  # Read-count sanity check (reads in R1 must equal reads in R2)
  R1=$(( $(zcat "${RAW_DIR}/${ID}_1.fastq.gz" | wc -l) / 4 ))
  R2=$(( $(zcat "${RAW_DIR}/${ID}_2.fastq.gz" | wc -l) / 4 ))
  [[ "$R1" -ne "$R2" ]] && { log "ERROR: ${ID} R1=${R1} R2=${R2} mismatch"; exit 1; }
  log "OK ${ID}: ${R1} read pairs"
done

# --- checksums so downloads can be verified later -----------------------------
( cd "$RAW_DIR" && ls *.fastq.gz >/dev/null 2>&1 && md5sum *.fastq.gz > checksums.md5 ) || true
log "Done. Checksums in ${RAW_DIR}/checksums.md5"
