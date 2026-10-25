# OncodriveFML test fixtures

These files are used only by `modules/local/oncodrivefml/tests/main.nf.test`.

| File | Content | Source | Terms |
| --- | --- | --- | --- |
| `cadd_subset.tsv.gz`, `cadd_subset.tsv.gz.tbi` | CADD v1.7 GRCh38 whole-genome SNV scores over the six synthetic fixture elements on chr5 (10,836 rows) | Extracted unchanged from `whole_genome_SNVs.tsv.gz`, CADD v1.7, <https://cadd.gs.washington.edu> | **Not covered by this repository's MIT licence.** CADD scores are freely available for non-commercial applications; commercial use requires a licence from the University of Washington. Cite Schubach et al. 2024, _Nucleic Acids Res_ 52:D1143–D1154, doi:10.1093/nar/gkad989. |
| `genome_hg38_chr5_head/chr5.txt` | First 1,400,000 bases of the GRCh38 (hg38) chromosome 5 sequence, in the per-chromosome text format that bgdata uses | bgdata genome reference hg38-20161209 | Public reference genome sequence |
| `oncodrivefml_test.conf` | Test copy of `assets/oncodrivefml.CADD_GRCh38.conf` (score file name and one core changed) | This repository | MIT |

The mutation table and element catalogue the test uses (`tests/data/activedriverwgs/`) are synthetic and contain no patient data.
