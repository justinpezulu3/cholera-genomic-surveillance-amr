#!/usr/bin/env bash
# 00_find_candidates.sh
# List publicly available V. cholerae genome ASSEMBLIES from a given country (or countries),
# with collection date, sequencing technology and assembly quality, to help choose isolates.
#
# Usage:   bash scripts/00_find_candidates.sh "Zambia" [tsv_name]
#          bash scripts/00_find_candidates.sh "Zambia,Kenya,Tanzania,Malawi" comparators
# Needs:   NCBI datasets CLI (environment chol-qc), python3, internet
# Output:  results/candidates/vibrio_cholerae_assemblies.jsonl   (downloaded once, tens of MB)
#          results/candidates/<tsv_name>.tsv                      (filtered candidates)
#
# Notes:
#  - Metadata (country, date) is whatever the submitter wrote in the BioSample; check it.
#  - 'sequencing_tech' tells you the platform (e.g. Illumina, Oxford Nanopore). Record it
#    in config/metadata.tsv: platform differences are a limitation to report.
#  - Do not pick isolates by looking only at this table. Read the source publication.

set -euo pipefail

COUNTRIES="${1:-Zambia}"
NAME="${2:-candidates}"
OUT="results/candidates"
JSONL="${OUT}/vibrio_cholerae_assemblies.jsonl"

mkdir -p "$OUT"

if [[ ! -s "$JSONL" ]]; then
  echo "Downloading V. cholerae (taxid 666) GenBank assembly summaries (one-off)..."
  datasets summary genome taxon 666 --assembly-source GenBank --as-json-lines > "$JSONL"
fi

python3 - "$COUNTRIES" "$JSONL" "${OUT}/${NAME}.tsv" <<'PY'
import csv, json, sys, collections

countries = [c.strip().lower() for c in sys.argv[1].split(",") if c.strip()]
jsonl, out_tsv = sys.argv[2], sys.argv[3]

def pick(d, *path, default="NA"):
    for p in path:
        if isinstance(d, dict) and p in d: d = d[p]
        else: return default
    return default if d in (None, "") else d

rows, scanned = [], 0
with open(jsonl) as f:
    for line in f:
        line = line.strip()
        if not line: continue
        r = json.loads(line); scanned += 1
        ai = r.get("assembly_info", {}) or {}
        bs = ai.get("biosample", {}) or {}
        attrs = {a.get("name"): a.get("value") for a in (bs.get("attributes") or []) if isinstance(a, dict)}
        geo  = bs.get("geo_loc_name")  or attrs.get("geo_loc_name")  or "NA"
        date = bs.get("collection_date") or attrs.get("collection_date") or "NA"
        src  = bs.get("isolation_source") or attrs.get("isolation_source") or "NA"
        if countries and not any(c in str(geo).lower() for c in countries):
            continue
        inf = pick(r, "organism", "infraspecific_names", default={})
        strain = (inf.get("strain") or inf.get("isolate") or "NA") if isinstance(inf, dict) else "NA"
        st = r.get("assembly_stats", {}) or {}
        rows.append({
            "accession": r.get("accession", "NA"),
            "strain_or_isolate": strain,
            "country": geo,
            "collection_date": date,
            "isolation_source": src,
            "sequencing_tech": ai.get("sequencing_tech", "NA") or "NA",
            "assembly_level": ai.get("assembly_level", "NA"),
            "contigs": st.get("number_of_contigs", "NA"),
            "contig_N50": st.get("contig_n50", "NA"),
            "total_length": st.get("total_sequence_length", "NA"),
            "biosample": bs.get("accession", "NA"),
            "bioproject": ai.get("bioproject_accession", "NA"),
            "release_date": ai.get("release_date", "NA"),
        })

rows.sort(key=lambda x: (str(x["country"]), str(x["collection_date"])))
cols = list(rows[0].keys()) if rows else ["accession"]
with open(out_tsv, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
    w.writeheader(); w.writerows(rows)

print(f"Scanned {scanned} assemblies; {len(rows)} match {countries}. Wrote {out_tsv}")
if rows:
    print("\nBy sequencing technology:")
    for k, v in collections.Counter(r["sequencing_tech"] for r in rows).most_common():
        print(f"  {v:4d}  {k}")
    print("\nBy collection year (first 4 characters of the date):")
    for k, v in sorted(collections.Counter(str(r["collection_date"])[:4] for r in rows).items()):
        print(f"  {k}: {v}")
PY

echo
echo "Next: open the TSV (e.g. 'column -t -s\$'\\t' ${OUT}/${NAME}.tsv | less -S'), then choose a small set."
echo "Check each candidate's BioProject and publication before using it."
