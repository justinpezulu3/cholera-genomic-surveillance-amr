#!/usr/bin/env bash
# 05_variants_phylogeny.sh
# Variant calling against the reference (Snippy), core SNP alignment (snippy-core),
# recombination masking (Gubbins), maximum-likelihood tree (IQ-TREE), SNP distances.
#
# Usage:   bash scripts/05_variants_phylogeny.sh [threads]
# Optional environment variables:
#   EXCLUDE="ID1 ID2"   isolates to leave out (e.g. FLAGged assemblies, or a distant outgroup)
#   OUTGROUP="OUT_01"   isolate ID used to root the tree in IQ-TREE (must be included)
# Input:   data/trimmed/<isolate>_{1,2}.fastq.gz, data/reference/*/ncbi_dataset/data/*/
# Output:  results/variants/snippy/<isolate>/            per-isolate variants
#          results/variants/variant_summary.tsv           SNP/indel counts per isolate
#          results/variants/watchlist_variants.tsv        variants in selected genes
#          results/phylogeny/core/                        snippy-core alignment + stats
#          results/phylogeny/gubbins/                     recombination-masked alignment
#          results/phylogeny/iqtree/cholera.treefile      final ML tree (+ bootstrap)
#          results/phylogeny/snp_dists_{raw,masked}.tsv   pairwise SNP distances
#          logs/phylo.log, logs/phylo_tool_versions.txt
#
# Notes on interpretation (state these in the report):
#  - Variants are called against ONE reference, so a distant isolate (e.g. a
#    non-O1/O139 outgroup) shrinks the shared core and can distort distances.
#    Consider running with and without it (use EXCLUDE).
#  - Trees from a handful of public isolates show relatedness among THOSE isolates.
#    They do not describe the cholera epidemic as a whole.

set -euo pipefail

THREADS="${1:-4}"
EXCLUDE="${EXCLUDE:-}"
OUTGROUP="${OUTGROUP:-}"

TRIM="data/trimmed"
VAR="results/variants"
PHY="results/phylogeny"
LOG="logs/phylo.log"

# Snippy parameters (report in Methods)
MINCOV=10
MINFRAC=0.9
MINQUAL=100

mkdir -p "$VAR/snippy" "$PHY"/{core,gubbins,iqtree} logs
log() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }

IQ=$(command -v iqtree2 || command -v iqtree || true)
[[ -z "$IQ" ]] && { echo "ERROR: iqtree2/iqtree not found in PATH." >&2; exit 1; }

{
  echo "# Variant/phylogeny tool versions ($(date '+%F %T'))"
  snippy --version 2>&1
  run_gubbins.py --version 2>&1 | head -n 1
  "$IQ" --version 2>&1 | head -n 1
  snp-sites -V 2>&1 | head -n 1 || true
  snp-dists -v 2>&1 | head -n 1 || true
} > logs/phylo_tool_versions.txt 2>&1 || true

# --- reference (prefer GenBank so variants are annotated) ----------------------
REF_GBK=$(ls data/reference/*/ncbi_dataset/data/*/genomic.gbff 2>/dev/null | head -n 1 || true)
REF_FNA=$(ls data/reference/*/ncbi_dataset/data/*/*.fna 2>/dev/null | head -n 1 || true)
REF="${REF_GBK:-$REF_FNA}"
[[ -z "$REF" ]] && { echo "ERROR: reference not found under data/reference/." >&2; exit 1; }
log "Reference for Snippy: $REF"

shopt -s nullglob
R1S=("$TRIM"/*_1.fastq.gz)
[[ ${#R1S[@]} -eq 0 ]] && { echo "ERROR: no trimmed reads in $TRIM." >&2; exit 1; }

IDS=()
for R1 in "${R1S[@]}"; do
  ID=$(basename "$R1" _1.fastq.gz)
  if [[ " $EXCLUDE " == *" $ID "* ]]; then log "EXCLUDED $ID"; continue; fi
  IDS+=("$ID")
done
# Pre-assembled genomes (no reads) join through snippy --ctgs
for A in data/assemblies/*.fasta; do
  [[ -e "$A" ]] || continue
  AID=$(basename "$A" .fasta)
  [[ -s "$TRIM/${AID}_1.fastq.gz" ]] && continue          # reads take priority
  if [[ " $EXCLUDE " == *" $AID "* ]]; then log "EXCLUDED $AID"; continue; fi
  # the reference genome itself is not a study isolate
  [[ "$AID" == "N16961_ref" ]] && continue
  IDS+=("$AID")
done
[[ ${#IDS[@]} -lt 3 ]] && { echo "ERROR: need at least 3 isolates plus the reference for a tree." >&2; exit 1; }

# --- 1. Snippy per isolate ------------------------------------------------------
for ID in "${IDS[@]}"; do
  if [[ -s "$VAR/snippy/$ID/snps.tab" ]]; then log "SKIP Snippy $ID"; continue; fi
  if [[ -s "$TRIM/${ID}_1.fastq.gz" ]]; then
    log "Snippy (reads): $ID"
    snippy --outdir "$VAR/snippy/$ID" --ref "$REF" \
      --R1 "$TRIM/${ID}_1.fastq.gz" --R2 "$TRIM/${ID}_2.fastq.gz" \
      --cpus "$THREADS" --mincov "$MINCOV" --minfrac "$MINFRAC" --minqual "$MINQUAL" \
      --force >> "$LOG" 2>&1
  else
    # assembled contigs (e.g. long-read assemblies): variants are called from contigs,
    # so depth/frequency filters do not apply. Record the platform in the report.
    log "Snippy (assembled contigs): $ID"
    snippy --outdir "$VAR/snippy/$ID" --ref "$REF" --ctgs "data/assemblies/${ID}.fasta" \
      --cpus "$THREADS" --force >> "$LOG" 2>&1
  fi
done

# --- 2. Core alignment ----------------------------------------------------------
log "snippy-core"
DIRS=(); for ID in "${IDS[@]}"; do DIRS+=("$VAR/snippy/$ID"); done
( cd "$PHY/core" && snippy-core --ref "$OLDPWD/results/variants/snippy/${IDS[0]}/ref.fa" --prefix core \
    $(for d in "${DIRS[@]}"; do echo "$OLDPWD/$d"; done) >> "$OLDPWD/$LOG" 2>&1 )

CORE_FULL="$PHY/core/core.full.aln"
snippy-clean_full_aln "$CORE_FULL" > "$PHY/core/clean.full.aln" 2>> "$LOG"

NSEQ=$(grep -c '>' "$PHY/core/clean.full.aln")
log "Alignment has $NSEQ sequences (isolates + reference)"

# --- 3. Recombination masking ---------------------------------------------------
if [[ "$NSEQ" -ge 4 ]]; then
  log "Gubbins"
  ( cd "$PHY/gubbins" && run_gubbins.py --prefix gubbins --threads "$THREADS" \
      "$OLDPWD/$PHY/core/clean.full.aln" >> "$OLDPWD/$LOG" 2>&1 )
  MASKED="$PHY/gubbins/gubbins.filtered_polymorphic_sites.fasta"
else
  log "WARN: fewer than 4 sequences; Gubbins skipped (recombination NOT masked)"
  MASKED="$PHY/core/clean.full.aln"
fi

snp-sites -c -o "$PHY/core/clean.core.aln" "$MASKED" >> "$LOG" 2>&1

# --- 4. SNP distances (before and after recombination masking) ------------------
snp-dists "$PHY/core/core.aln" > "$PHY/snp_dists_raw.tsv" 2>> "$LOG"
snp-dists "$PHY/core/clean.core.aln" > "$PHY/snp_dists_masked.tsv" 2>> "$LOG"

# --- 5. IQ-TREE (ascertainment-bias correction because only variable sites kept) -
log "IQ-TREE"
OG_ARGS=(); [[ -n "$OUTGROUP" ]] && OG_ARGS=(-o "$OUTGROUP")
"$IQ" -s "$PHY/core/clean.core.aln" -m MFP+ASC -B 1000 --alrt 1000 \
  -T "$THREADS" --seed 12345 --prefix "$PHY/iqtree/cholera" "${OG_ARGS[@]}" -redo >> "$LOG" 2>&1

# --- 6. Variant summaries -------------------------------------------------------
python3 - <<'PY'
import csv, glob, os, collections

WATCH = {"gyrA", "gyrB", "parC", "parE", "ctxB", "toxR"}   # edit as needed
summ, watch = [], []

for p in sorted(glob.glob("results/variants/snippy/*/snps.tab")):
    iso = os.path.basename(os.path.dirname(p))
    types = collections.Counter()
    with open(p, newline="") as f:
        for d in csv.DictReader(f, delimiter="\t"):
            types[d.get("TYPE", "NA")] += 1
            gene = d.get("GENE", "")
            if gene in WATCH:
                watch.append([iso, gene, d.get("CHROM"), d.get("POS"), d.get("TYPE"),
                              d.get("REF"), d.get("ALT"), d.get("EFFECT"), d.get("AA_POS", "")])
    summ.append([iso, types.get("snp", 0), types.get("ins", 0), types.get("del", 0),
                 types.get("complex", 0), types.get("mnp", 0), sum(types.values())])

with open("results/variants/variant_summary.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["isolate", "snp", "ins", "del", "complex", "mnp", "total"])
    w.writerows(summ)

with open("results/variants/watchlist_variants.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["isolate", "gene", "chrom", "pos", "type", "ref", "alt", "effect", "aa_pos"])
    w.writerows(watch)

print("Wrote variant_summary.tsv and watchlist_variants.tsv")
print("Core genome stats per isolate: results/phylogeny/core/core.txt")
PY

log "Done. Check core.txt (aligned fraction), Gubbins recombination GFF, and bootstrap support before interpreting the tree."
log "Root the tree with the outgroup in R (ggtree) and keep support values visible on the figure."
