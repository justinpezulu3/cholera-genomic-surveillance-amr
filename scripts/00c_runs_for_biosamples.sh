#!/usr/bin/env bash
# 00c_runs_for_biosamples.sh
# For genome assemblies already listed in a candidates table (from 00_find_candidates.sh),
# check whether public raw SRA reads exist for each assembly's BioSample, and what they are.
# This answers: "can I run this isolate through the full read-based pipeline?"
#
# Usage:   bash scripts/00c_runs_for_biosamples.sh TABLE.tsv BIOPROJECT [REGEX] [PER_GROUP]
#   TABLE.tsv   output of 00_find_candidates.sh, e.g. results/candidates/Zambia.tsv
#   BIOPROJECT  assembly BioProject to examine, e.g. PRJNA701309
#   REGEX       optional filter on the country/province column, e.g. "Lusaka|Eastern"
#   PER_GROUP   keep at most N assemblies per (province, year), best (fewest contigs) first.
#               Default 2. Use 0 to query all.
# Needs:   entrez-direct (esearch, efetch, xtract) in chol-qc, python3, internet
# Output:  results/candidates/<BIOPROJECT>_biosample_runs.tsv
# Tip:     export NCBI_API_KEY=... (free from your NCBI account) for fewer connection errors.
#
# Queries are retried when the connection drops. A sample reported "FAILED" means the
# connection failed, NOT that reads do not exist. Rerun the script to retry those.

set -euo pipefail

TABLE="${1:?Usage: bash scripts/00c_runs_for_biosamples.sh TABLE.tsv BIOPROJECT [REGEX] [PER_GROUP]}"
BP="${2:?BioProject accession required}"
REGEX="${3:-}"
PER_GROUP="${4:-2}"
OUT="results/candidates/${BP}_biosample_runs.tsv"
mkdir -p results/candidates

python3 - "$TABLE" "$BP" "$REGEX" "$PER_GROUP" "$OUT" <<'PY'
import collections, csv, re, subprocess, sys, time

table, bp, regex, per_group, out = sys.argv[1:6]
per_group = int(per_group)
rx = re.compile(regex) if regex else None

def to_int(x, default=10**9):
    try: return int(x)
    except Exception: return default

rows = []
with open(table, newline="") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if r.get("bioproject") != bp: continue
        if r.get("strain_or_isolate", "").startswith("ena-yuan"): continue   # duplicate re-submissions
        if rx and not rx.search(r.get("country", "")): continue
        rows.append(r)

rows.sort(key=lambda r: (r["country"], str(r["collection_date"])[:4], to_int(r["contigs"])))
if per_group > 0:
    kept, seen = [], collections.Counter()
    for r in rows:
        key = (r["country"], str(r["collection_date"])[:4])
        if seen[key] < per_group:
            kept.append(r); seen[key] += 1
    rows = kept
print(f"{len(rows)} BioSamples to query in {bp}", flush=True)

def sh(cmd, timeout=120):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=timeout)

def count_runs(bs):
    """Return number of SRA records for the BioSample, or None if the connection failed."""
    for attempt in range(4):
        try:
            r = sh(f'esearch -db sra -query "{bs}[BioSample]" | xtract -pattern ENTREZ_DIRECT -element Count')
            txt = r.stdout.strip().split("\n")[-1].strip() if r.stdout.strip() else ""
            if txt.isdigit(): return int(txt)
        except Exception:
            pass
        time.sleep(2 * (attempt + 1))
    return None

def get_runinfo(bs):
    for attempt in range(4):
        try:
            r = sh(f'esearch -db sra -query "{bs}[BioSample]" | efetch -format runinfo')
            lines = [l for l in r.stdout.splitlines() if l.strip()]
            if lines and lines[0].startswith("Run,"):
                return list(csv.DictReader(lines))
        except Exception:
            pass
        time.sleep(2 * (attempt + 1))
    return None

result = []
for i, r in enumerate(rows, 1):
    bs = r["biosample"]
    n = count_runs(bs)
    if n is None:
        status, runs = "FAILED (connection)", []
    elif n == 0:
        status, runs = "no_sra_reads", []
    else:
        info = get_runinfo(bs)
        if info is None:
            status, runs = "FAILED (connection)", []
        else:
            status, runs = "reads_found", info
    plat = ",".join(sorted({x.get("Platform", "") for x in runs})) or "NA"
    lay = ",".join(sorted({x.get("LibraryLayout", "") for x in runs})) or "NA"
    try: depth = round(sum(float(x["bases"]) for x in runs) / 4.0e6, 1) if runs else "NA"
    except Exception: depth = "NA"
    result.append([r["accession"], r["strain_or_isolate"], r["country"], r["collection_date"],
                   bs, status, ",".join(x["Run"] for x in runs) or "NA", plat, lay, depth,
                   r["contigs"], r["contig_N50"], r["total_length"]])
    print(f"[{i}/{len(rows)}] {bs} {r['strain_or_isolate']}: {status} {','.join(x['Run'] for x in runs)}", flush=True)
    time.sleep(0.4)

with open(out, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["assembly", "strain", "country", "collection_date", "biosample", "status",
                "sra_runs", "platform", "layout", "est_depth_x", "asm_contigs", "asm_N50", "asm_length"])
    w.writerows(result)

c = collections.Counter(x[5] for x in result)
print("\nSummary:", dict(c))
print(f"Wrote {out}")
if c.get("FAILED (connection)"):
    print("Some queries failed on the connection. Rerun the same command to retry.")
PY
