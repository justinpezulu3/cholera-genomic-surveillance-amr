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
