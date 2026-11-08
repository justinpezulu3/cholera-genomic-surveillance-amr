#!/usr/bin/env bash
# 00b_illumina_runs.sh
# For a BioProject, list Illumina paired-end WGS runs with country, collection date and
# estimated depth, so you can pick comparator isolates that suit the read-based pipeline.
#
# Usage:   bash scripts/00b_illumina_runs.sh PRJNA616030
# Needs:   entrez-direct (esearch, efetch) in environment chol-qc, python3, internet
# Output:  results/candidates/<BioProject>_illumina_runs.tsv
#
# Estimated depth = bases / 4.0 Mb. Aim for roughly >=50x for assembly. Very low depth
# runs will assemble badly and should be avoided.
# Optional: export NCBI_API_KEY=... for faster queries (free from your NCBI account).

set -euo pipefail

BP="${1:?Usage: bash scripts/00b_illumina_runs.sh PRJNAxxxxxx}"
OUT="results/candidates"
RUNINFO="${OUT}/${BP}_runinfo.csv"
mkdir -p "$OUT"

if [[ ! -s "$RUNINFO" ]]; then
  echo "Fetching run table for ${BP}..."
  esearch -db sra -query "${BP}[BioProject]" | efetch -format runinfo > "$RUNINFO"
fi

python3 - "$BP" "$RUNINFO" "${OUT}/${BP}_illumina_runs.tsv" <<'PY'
import csv, subprocess, sys, time
import xml.etree.ElementTree as ET

bp, runinfo, out_tsv = sys.argv[1:4]

with open(runinfo, newline="") as f:
    runs = [r for r in csv.DictReader(f) if r.get("Run")]
print(f"{len(runs)} runs in {bp}")

keep = [r for r in runs
        if r.get("Platform", "").upper() == "ILLUMINA"
        and r.get("LibraryLayout", "").upper() == "PAIRED"
        and r.get("LibraryStrategy", "").upper() in ("WGS", "OTHER", "")]
print(f"{len(keep)} are Illumina paired-end WGS")

cache = {}
def biosample_meta(bs):
    if bs in cache: return cache[bs]
    meta = {"geo_loc_name": "NA", "collection_date": "NA", "isolation_source": "NA"}
    try:
        res = subprocess.run(["efetch", "-db", "biosample", "-id", bs, "-format", "xml"],
                             capture_output=True, text=True, timeout=60)
        root = ET.fromstring(res.stdout)
        for a in root.iter("Attribute"):
            for key in (a.get("harmonized_name"), a.get("attribute_name")):
                if key in meta and a.text and meta[key] == "NA":
                    meta[key] = a.text.strip()
    except Exception as e:
        print(f"  warning: could not read {bs}: {e}", file=sys.stderr)
    cache[bs] = meta
    time.sleep(0.4)   # be polite to NCBI
    return meta

rows = []
for r in keep:
    bs = r.get("BioSample", "NA")
    m = biosample_meta(bs) if bs != "NA" else {"geo_loc_name": "NA", "collection_date": "NA", "isolation_source": "NA"}
    try: depth = round(float(r["bases"]) / 4.0e6, 1)
    except Exception: depth = "NA"
    rows.append([r["Run"], bs, r.get("Model", "NA"), r.get("spots", "NA"), r.get("bases", "NA"),
                 depth, m["geo_loc_name"], m["collection_date"], m["isolation_source"]])

rows.sort(key=lambda x: (str(x[6]), str(x[7])))
with open(out_tsv, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["run", "biosample", "instrument", "spots", "bases", "est_depth_x",
                "country", "collection_date", "isolation_source"])
    w.writerows(rows)
print(f"Wrote {out_tsv}")
PY
