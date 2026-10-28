# Decisions log

Record every choice that could change a result. Date each entry. Reviewers value an honest log more than a perfect pipeline.

## 1. Isolate selection

- **Search performed:** [date, database (NCBI SRA / ENA / Pathogen Detection), query terms, filters]
- **Candidates found:** [number]
- **Inclusion criteria:** Illumina paired-end; collection year and country available; [others]
- **Exclusion criteria and what was excluded:** [list accession and reason]
- **Final set:** see `config/metadata.tsv` (one-line rationale per isolate)
- **Known sampling biases:** [e.g. outbreak isolates over-represented; single laboratory; years missing]

## 2. Reference genome

- N16961, RefSeq GCF_000006745.1. Reason: [7th-pandemic El Tor reference, widely used]. Implication: variants are relative to this genome; divergent isolates align poorly.

## 3. Quality thresholds

| Item | Value used | Source / justification |
|---|---|---|
| fastp quality / min length | Q20 window, 50 bp | [after inspecting raw MultiQC] |
| Genome size window | 3.8-4.4 Mb | expected *V. cholerae* size |
| Max contigs | 300 | |
| Min N50 | 50 kb | |
| CheckM2 completeness / contamination | >=95% / <=5% | |
| fastANI vs reference | >=95% | species threshold |
| Snippy mincov / minfrac / minqual | 10 / 0.9 / 100 | |
| ABRicate min identity / coverage | 90% / 80% | |

## 4. Flagged isolates

| Isolate | Flag(s) | Decision (keep / exclude / keep with caveat) | Reason |
|---|---|---|---|
| | | | |

## 5. Phylogeny choices

- Outgroup used: [yes/no, which]. Run with and without outgroup? [result of comparison]
- Recombination masked with Gubbins: [yes/no; number of sequences]
- Model selected by ModelFinder: [from `results/phylogeny/iqtree/cholera.iqtree`]
- Bootstrap: ultrafast, 1000 replicates; SH-aLRT 1000

## 6. AMR interpretation

- Databases and versions: AMRFinderPlus [version, DB version]; ABRicate resfinder [date], card [date]
- Disagreements between tools and how they were resolved: [list]
- Phenotype data available for any isolate? [yes/no; source]. If yes, genotype-phenotype comparison: [summary]

## 7. Deviations from the plan

- [date] [what changed and why]
