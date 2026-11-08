# Results (draft, check every number against results/ before using)

## Dataset
Sixteen public *Vibrio cholerae* isolates were analysed against the N16961 reference: 13 from sequencing reads and 3 pre-assembled Lusaka genomes. One further isolate (TZA_2017_DarEsSalaam) was excluded after assembly QC. The set is a convenience sample and does not represent Zambian cholera as a whole.

## Core-genome alignment
Each isolate aligned to about 96-97% of the 4,033,464 bp reference (3.86-3.92 Mb). After recombination masking with Gubbins, the alignment contained 247 variable sites, of which 63 were parsimony-informative. IQ-TREE selected K3P+ASC.

## Relationships (SNP distances, recombination masked)
- Isolates were 100-133 SNPs from N16961.
- A tight older group (KEN_2007, KEN_2009, KEN_2015, GLB_HTI_2010, TZA_2016, UGA_2018 and the Zambian isolates from 2016, 2018 and 2020) differed by mostly 0-9 SNPs (up to 18 including KEN_2009).
- The 2023-24 Zambian isolates fell into two separate groups: EAS_3712, NOR_3779 and CEN_3757 (1-4 SNPs apart) and LUS_3598 and COP_3752 (2 SNPs apart). The two groups were about 59-63 SNPs apart, and each was about 28-43 SNPs from the older group.
- ZMB_2009_Lus38 lay between them (about 21-33 SNPs from the older group). GLB_BGD_2009 was the most distant isolate (37-70 SNPs from the others).
- Interpretation: the 2023-24 Zambian isolates come from at least two distinct lineages. With five isolates chosen by us, this says nothing about how common either lineage is.

## Tree support
Only a few clades were well supported (ultrafast bootstrap 100 for each 2023-24 group, 98 for ZMB_2018_Lus7 and ZMB_2020_ZTox). Most other internal nodes had low support, which is expected with 63 informative sites, so the branching order within the older group should not be interpreted.

## Resistance-associated determinants (genotype only)
- Present in all isolates: almEFG, catB9, dfrA1, varG and the gyrA S83I mutation.
- parC S85L was present in every isolate except ZMB_2009_Lus38, which carries the wild-type base at that position (confirmed as a real base call in the alignment).
- aph(3'')-Ib, aph(6)-Id, sul2 and floR were found in KEN_2007, KEN_2009, GLB_HTI_2010, ZMB_2023_LUS_3598, ZMB_2024_COP_3752 and ZMB_2009_Lus38. GLB_BGD_2009 carried a different set (aph(3'')-Ib, aph(6)-Id, sul2, dfrA31, qnrVC1, tet(59)).
- Absence of a determinant in an assembled genome may reflect a fragmented mobile element.
- Two partial hits (blaZ, blaI) were not counted.
- No phenotypic data were available, so none of these findings establishes resistance.

## Limitations
- Small, non-random sample; five Zambian isolates from 2023-24.
- Low resolution: 247 variable sites.
- N16961 is a pragmatic root, not a true outgroup.
- Three genomes entered as assemblies, not reads.
- KEN_2007 had low depth and many heterozygous calls.
- AMRFinderPlus database from July 2024.
