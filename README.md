# Cholera Genomic Surveillance and Antimicrobial Resistance Analysis

A reproducible bacterial genomics workflow applied to a small, deliberately selected set of publicly available *Vibrio cholerae* genomes.

> **Status:** template. Fill the bracketed fields after you have chosen isolates and run the pipeline.

## Research question

What genomic differences distinguish selected *Vibrio cholerae* isolates, and what resistance-associated genetic determinants can be identified in their genomes?

## Design in one paragraph

[N] public isolates were selected from NCBI SRA/ENA: [n] from Zambia, [n] from other African countries, [n] from other regions, plus the reference genome N16961 and [an outgroup / no outgroup]. Selection criteria and per-isolate rationale are in `config/metadata.tsv` and `docs/decisions.md`. The sample is a convenience sample of public data and **does not represent all cholera cases in Zambia or elsewhere.**

## Headline result

`results/figures/fig_tree_amr.png`: [one-sentence description of what the figure shows, written after you have looked at it]

## Workflow

| Step | Script | Main tools | Key outputs |
|---|---|---|---|
| 1. Retrieve data | `scripts/01_download_data.sh` | SRA Toolkit, NCBI datasets | `data/raw/`, `data/reference/`, checksums |
| 2. Read QC | `scripts/02_qc.sh` | FastQC, fastp, MultiQC | `results/qc/` |
| 3. Assembly + QC | `scripts/03_assembly.sh` | Shovill/SPAdes, QUAST, CheckM2, fastANI, MLST | `results/assembly/`, `results/assembly_qc/` |
| 4. Annotation + AMR | `scripts/04_amr_annotation.sh` | Bakta, AMRFinderPlus, ABRicate | `results/amr/` |
| 5. Variants + phylogeny | `scripts/05_variants_phylogeny.sh` | Snippy, Gubbins, IQ-TREE 2 | `results/variants/`, `results/phylogeny/` |
| 6. Figures | `R/06_visualise.R` | ggtree, pheatmap | `results/figures/`, `results/tables/` |
| 7. Report | `report/report.qmd` | Quarto | `report/report.html` |

The workflow was first explored with Galaxy's bacterial genomics and AMR tutorials (Galaxy Training Network) and then reproduced with the command-line tools above.

## Quick start

```bash
# 1. environments (one small environment per pipeline stage)
bash create_envs.sh

# 2. one-time databases (see headers of scripts 03 and 04)
conda run --no-capture-output -n chol-checkm2 checkm2 database --download --path ~/databases/checkm2
conda run --no-capture-output -n chol-amr bakta_db download --output ~/databases/bakta --type light
conda run --no-capture-output -n chol-amr amrfinder -u
export CHECKM2DB=~/databases/checkm2/CheckM2_database/uniref100.KO.1.dmnd
export BAKTA_DB=~/databases/bakta/db-light

# 3. fill in config/metadata.tsv (accessions, country, year, source)

# 4. run
bash run_all.sh 4          # 4 threads
quarto render report/report.qmd
```

Useful options: `FROM=3 bash run_all.sh 4` resumes at step 3; `EXCLUDE="ID1 ID2"` and `OUTGROUP="ID"` are passed to step 5.

## Repository layout

```
config/metadata.tsv     isolate table (single source of truth)
scripts/                shell steps 01-05
R/                      figures and tables
report/report.qmd       scientific report
docs/decisions.md       isolate selection, thresholds, deviations
environment.yml         software environment (freeze to environment.lock.yml)
logs/                   run logs and tool/database versions
results/                derived outputs (figures and small tables are tracked)
data/                   raw reads and references (NOT tracked; regenerate from accessions)
```

## How to interpret the outputs

- A detected gene or mutation is **evidence of a resistance-associated genotype**, not proof of phenotypic resistance. Phenotype requires susceptibility testing; where public phenotype data exist it is compared explicitly, otherwise this is stated.
- Absence of a determinant in a short-read assembly does not prove absence from the isolate; mobile elements such as the SXT/ICE often fragment.
- The tree shows relatedness among the selected isolates only.
- All thresholds are listed in script headers and `docs/decisions.md`.

## Limitations

[Write after analysis. At minimum: small convenience sample; public-data sampling bias; short-read assembly limits; single-reference variant calling; genotype versus phenotype; database versions.]

## Reproducibility

- Software: `environment.lock.yml`, `logs/*_tool_versions.txt`, `logs/R_sessionInfo.txt`
- Databases: versions recorded by scripts 03 and 04 (record the download date)
- Data: accession list in `config/metadata.tsv`; checksums in `data/raw/checksums.md5`
- Random seeds: IQ-TREE seed fixed in script 05

## Data sources

All sequence data are public (NCBI SRA/ENA, RefSeq). Cite each BioProject and the original publication where one exists. See `config/metadata.tsv`.

## Citation and licence

Code: MIT (add a `LICENSE` file). Cite the tools you used: [fill from the report's Methods].

## Author

Justin, [institution], [year]. [Contact / portfolio link]
