# Decisions log

Record every choice that could change a result. Date each entry. Reviewers value an honest log more than a perfect pipeline.

## 1. Isolate selection

- **Search performed:** [date, database (NCBI SRA / ENA / Pathogen Detection), query terms, filters]
- **Candidates found:** [number]
- **Inclusion criteria:** Illumina paired-end; collection year and country available; [others]
- **Exclusion criteria and what was excluded:** [list accession and reason]
- **Final set:** see `config/metadata.tsv` (one-line rationale per isolate)
- **Known sampling biases:** [e.g. outbreak isolates over-represented; single laboratory; years missing]

## 1b. Hybrid data design (decided)

Most public Zambian genomes are long-read (Nanopore) data, while the pipeline's read-level steps assume Illumina paired-end reads. Design chosen:

- **Illumina comparators (data_type = reads):** full pipeline from raw reads (QC, assembly, QC gates, AMR, variants).
- **Zambian genomes (data_type = assembly):** pre-assembled GenBank genomes that enter at assembly QC, AMR detection and phylogeny (variants called from contigs).
- Consequences to state in the report: read QC and assembly were not performed by this project for the assembled genomes; platform differences (Illumina vs Nanopore) can bias SNP calls and AMR hit fragmentation; assemblies were produced by the original authors with their own methods.
- Platform for every isolate is recorded in `config/metadata.tsv`.
- Sensitivity check planned: compare tree topology and SNP distances with and without the assembled genomes.

## 1c. Isolate shortlist (decided 29 Oct 2026)

Candidate search: 295 Zambian GenBank assemblies found with `00_find_candidates.sh` (taxid 666); read availability checked per BioSample with `00c_runs_for_biosamples.sh`; East African Illumina runs listed with `00b_illumina_runs.sh` (PRJNA616030).

Findings that drove the design:
- PRJNA701309 (Lusaka/Chiengi 2009-2018, 18 assemblies checked): **no public SRA reads linked to any BioSample**, so these enter as pre-assembled genomes only.
- PRJNA514245 (2020-2024, 29 isolates checked): Illumina paired-end reads available for all, 71-614x depth, so these run through the full read-based pipeline.
- PRJNA1228025 (2023, 32 complete Nanopore genomes): not used in the main analysis (platform mismatch, near-identical isolates).
- PRJEB71012: duplicate re-submissions of PRJNA701309 assemblies (`ena-yuan-...`); excluded to avoid double counting.

Selection rules: spread across years and provinces; fewest contigs among assemblies; depth about 70-250x for reads (very deep runs avoided to limit download size); one isolate per province-time point.

Excluded older assemblies: 44_Lusaka_2009 (2,844 contigs, 4.58 Mb), 58 and 59_Lusaka_2010 (about 2,100-2,200 contigs, about 4.54-4.59 Mb), Chiengi_2017 (1,514 contigs): sizes and fragmentation suggest contamination or poor data.

Tree rooting: no outgroup is used; the tree is rooted on N16961 (1975), an older 7th-pandemic reference. Limitation: rooting on an in-lineage strain gives relative, not absolute, ordering.

Sampling caveat to report: isolates were chosen for spread and data availability, not at random, from public submissions.

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

| TZA_2017_DarEsSalaam | genome_size; fragmented; low_N50; contamination (14.1%) | exclude | Likely mixed culture: 4.47 Mb, ST unresolved; AMR/SNP calls unreliable |
| UGA_2018_Kampala | low_N50; contamination (7.8%) | keep with caveat | ANI 99.98%; sensitivity check with/without in tree |
| ZMB_2018_Lus7 | fragmented; contamination (5.2%) | keep with caveat | Pre-assembled by authors; marginal |
| KEN_2007_HomaBay | low_N50 (31 kb) | keep with caveat | ~56x depth; contamination 2.4%, ANI 99.97% |
| ZMB_2009_Lus38, ZMB_2016_Lus1033 | fragmented | keep with caveat | Pre-assembled by authors; N50 about 199 kb |

## 5. Phylogeny choices

- Outgroup used: [yes/no, which]. Run with and without outgroup? [result of comparison]
- Recombination masked with Gubbins: [yes/no; number of sequences]
- Model selected by ModelFinder: [from `results/phylogeny/iqtree/cholera.iqtree`]
- Bootstrap: ultrafast, 1000 replicates; SH-aLRT 1000

## 6. AMR interpretation

- Databases and versions: AMRFinderPlus [version, DB version]; ABRicate resfinder [date], card [date]
- Disagreements between tools and how they were resolved: [list]
- Phenotype data available for any isolate? [yes/no; source]. If yes, genotype-phenotype comparison: AMRFinderPlus 3.12.8 (DB 2024-07-22.1, now outdated) and ABRicate (ResFinder, CARD; 2026-Apr-3). Partial hits (<80% coverage) are not counted. Calls are resistance-associated genotype, not phenotype. No phenotype data available.

## 7. Deviations from the plan

- [04/11/2026]
- Bakta 1.9.4 crashed with pyhmmer 0.10.15 (hit.name is text, not bytes); patched bakta/features/orf.py to accept both.

## Decisions added during assembly, AMR and phylogeny (Nov 2026)

### Assembly
- Shovill was dropped. It requires `--ram` >= 8 and <= physical RAM, which cannot be satisfied on a laptop where WSL2 has about 4.8 GB.
- Replaced by: seqtk subsampling to about 100x (seed 100, same fraction and seed for R1 and R2), then SPAdes `--isolate` with `-m 8`. Contigs < 200 bp or < 5x coverage are removed. See `scripts/03_assembly.sh`.
- `setuptools` was added to `chol-asm` (SPAdes and QUAST import `distutils`, removed in Python 3.12+).

### Isolate selection and exclusions
- 16 isolates entered the analysis: 13 read-based and 3 pre-assembled Lusaka genomes (PRJNA701309), plus the N16961 reference (GCF_000006745.1).
- TZA_2017_DarEsSalaam was excluded after assembly QC flagged it. [Fill in the exact reason and metrics from results/assembly_qc/assembly_qc_summary.tsv.]
- KEN_2007_HomaBay was kept but flagged: low depth after trimming (about 56x) and many heterozygous calls (8,808).
- Selection is a convenience set of public isolates. It does not represent all Zambian cholera.

### Bakta
- Bakta 1.9.4 on Python 3.10 crashed in the CDS step: `hit.name.decode()` fails when pyhmmer returns text instead of bytes.
- Local patch to `bakta/features/orf.py` (4 replacements) and `cds.py` (3 replacements): each `.decode()` call is wrapped in an `isinstance(..., bytes)` check. The original is kept as `orf.py.bak`.

### AMR interpretation
- Hits are called resistance-associated determinants. Genotype is not phenotype, and no phenotype data were available.
- Partial hits are not counted. AMRFinderPlus `PARTIAL_CONTIG_END` hits are excluded from figures and conclusions. Example: blaZ in ZMB_2016_Lus1033 (70.11% coverage, 100% identity) and blaI in ZMB_2009_Lus38 (78.57% coverage, 98.99% identity), both on small contigs of public draft assemblies.
- ctxA and ctxB are virulence markers and are not shown in the AMR panels.
- ABRicate reports one gene under several database entries (for example `sul2` and `sul2_2`). Columns differing only by a trailing `_<number>` are merged, and CARD's `vibrio_cholerae_varG` is merged with `varG`.
- gyrA S83I and parC S85L are point mutations: ABRicate does not call them, so they appear as "AMRFinderPlus only". That is not a disagreement between tools.
- Database versions: AMRFinderPlus 3.12.8 with database 2024-07-22.1 (outdated, a stated limitation); ABRicate resfinder and card, 2026-Apr-3.

### Variants and phylogeny
- Snippy 4.6.0 runs from its own environment (`chol-snippy`); script 05 appends `chol-phylo/bin` to PATH.
- Reference: N16961 (GCF_000006745.1) GenBank file. `snippy-core` must be given Snippy's own `ref.fa` (sequence names `NC_002505`, `NC_002506`). The NCBI `.fna` uses versioned names (`NC_002505.1`) and caused "has entry ... but reference doesn't".
- Gubbins 3.4.3 masks recombination; IQ-TREE 2 uses `-m MFP+ASC` (best fit K3P+ASC) with 1000 ultrafast bootstraps and 1000 SH-aLRT replicates.
- The tree is rooted on N16961 for display. It is a pragmatic root (the reference is itself a 7th-pandemic El Tor isolate), not a true outgroup.
- Ultrafast bootstrap labels are shown only when >= 70.

### Storage
- Raw FASTQ files were deleted to recover disk space on D: (the WSL virtual disk). They can be regenerated with `scripts/01_download_data.sh`.
