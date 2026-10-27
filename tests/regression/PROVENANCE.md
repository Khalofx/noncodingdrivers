# Regression baseline: provenance

The regression gates compare v1 runs with the frozen tables described here. Everything is
identified by content (checksums), not by repository or tag. The frozen tables are kept on the
project cluster under `noncodingdrivers/data/baseline/` (read-only); there is no public deposit.

## Inputs

| Input                                                                                                               | Checksum                                                                                                                             |
| ------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------ |
| PCAWG public somatic MAF (hg19)                                                                                     | MD5 `009b9f88ec10d07cf826174c3dea5b78`                                                                                               |
| Liver-HCC-public canonical mutation table (GRCh38; 3,379,254 rows)                                                  | SHA-256 `39bf4a0a3fe322a78abfb1e6b0a5b8d45800ecea7ae0caf1223c14111798f0cf`                                                           |
| LIRI-JP canonical mutation table (GRCh38)                                                                           | SHA-256 `4000e4fe22e2c9642bc158d2d0ead0e554927097ba785f76591538111557742d`                                                           |
| PACA-CA canonical mutation table (GRCh38)                                                                           | SHA-256 `5da1f75690a79b04193d9958f9f3100d406eb46555aecc63fdac24c4ff93612d`                                                           |
| Element catalogue: Ensembl Regulatory Build 110, CDS removed against GENCODE v50 basic (595,725 rows, 544,737 ids)  | MD5 `39f838e6c81fdb0b152504f2d51939b5`                                                                                               |
| Regulatory Build 110 GFF (20221007)                                                                                 | MD5 `ee3b05c9020f7da4244d563995a579bb` (re-checked in Step 5; the recipe in docs/references.md rebuilds the catalogue byte for byte) |
| GRCh38 FASTA, chr-named, 25 contigs (`.fai` MD5 `963c862ff26b0d914956faaf7192e54b`)                                 | MD5 `838cd9f56ca7d7dd345cb025d188c991`                                                                                               |
| GRCh38 sequence dictionary `derived/hg38.chr.dict` (samtools 1.23.1 `dict`; Picard LiftoverVcf)                     | MD5 `eee72059c0f8b6407202baa38ba3be8a`                                                                                               |
| GENCODE v50 basic GTF                                                                                               | MD5 `11e77cf1d4b78b43f57d608e3f7f2202`                                                                                               |
| JASPAR 2026 CORE vertebrates (1,019 PFMs)                                                                           | MD5 `c63218acf88eb5fe10255b7c8ef23d48`                                                                                               |
| hg19ToHg38 liftover chain                                                                                           | MD5 `35887f73fe5e2231656504d1f6430900`                                                                                               |
| UCSC hg19 FASTA `hg19/hg19.fa`, 93 contigs (`hg19.fa.gz` MD5 `806c02398f5ac5da8ffd6da2d1d5d1a9`, as published)      | MD5 `530d89d3ef07fdb2a9b3c701fb4ca486`; `.fai` MD5 `e700bec606168273af0920bd1dd49eb9`                                                |
| Contig rename map `assets/references/GRCh38.rename_chrs.map` (479 rows: Ensembl names + RefSeq map of the old repo) | MD5 `91ac4e7b7e10ae19a8ffbe5270642711`                                                                                               |
| CADD GRCh38 v1.7 whole-genome SNVs, `whole_genome_SNVs.tsv.gz` (user-supplied, not redistributed)                   | MD5 `88577a55f1cd519d44e0f415ba248eb9`; `.tbi` MD5 `347df8fac17ea374c4598f4f44c7ce8b`                                                |
| Ensembl VEP 116 cache, `homo_sapiens_vep_116_GRCh38.tar.gz`                                                         | MD5 `6e10959391069a9eebb8d8d165fd29be`                                                                                               |
| OncodriveFML genome: official bgdata `genomereference/hg38-20161209`                                                | copied unchanged into `data/refs/oncodrivefml_genome/`                                                                               |
| Variant-flag mask: ENCODE Blacklist v2 hg38 `masks/hg38-blacklist.v2.bed.gz` (Boyle-Lab/Blacklist commit 61a04d2c)  | MD5 `83fe6bf8187a64dee8079b80f75ba289`                                                                                               |
| Variant-flag mask: GIAB v3.6 GRCh38 segmental duplications `masks/GRCh38_segdups.bed.gz`                            | MD5 `72ec3f7db18ef9318c64b1df271703fe`                                                                                               |
| Variant-flag mask: GIAB v3.6 GRCh38 low mappability (all) `masks/GRCh38_lowmappabilityall.bed.gz`                   | MD5 `8a78766bb99abf1d25a5b307d49a0324`                                                                                               |
| Variant-flag mask: Umap k100 single-read mappability hg38 `masks/k100.umap.bed.gz`                                  | MD5 `1ef87107f449ccaf7e537de9740692fb`                                                                                               |
| gnomAD v4.1 genomes echtvar archive `gnomad/gnomad.v4.1.genomes.echtvar.zip` (echtvar 0.2.4, 759,336,320 variants; field configuration MD5 `06f2ebddcf38deebfbd8873ed5e40986`) | MD5 `7bc89a78568a424d7694e2d351617bcd` |

hg19 MAF to GRCh38 liftover of Liver-HCC-public (the baseline canonical table): 3,380,218 records in,
3,379,254 out; 875 failed (176 unmapped, 2 split mappings, 697 REF mismatches); 65,116 insertions
carried without REF validation. Variant-type labels in the table: SNP, DEL, INS, DNP, TNP, ONP.

## Tools

ActiveDriverWGS 1.2.1 (R 4.5.3 site library, no container; BSgenome.Hsapiens.UCSC.hg38 1.4.5);
OncodriveFML 2.5.0 (seed 123); Nextflow 25.10.4; python 3.11, pandas 2.3.3, numpy 1.26.4, scipy 1.17.1;
Ensembl VEP 116.2; MultiQC 1.21. MutEnricher and motifbreakR did not contribute to the frozen tables.

Images used by the baseline runs (from the task scripts of the 2026-10-06 runs), by SIF SHA-256:

| Role                 | SHA-256                                                            |
| -------------------- | ------------------------------------------------------------------ |
| utils (python tools) | `ff62ef3b67216a5a1fc4acc7705409a5225cf215f24e02b38981c12632eb677b` |
| OncodriveFML 2.5.0   | `fd60e15863c41fc75c614f5616ca3a5711174719fe7d5c43026f003512661a5d` |
| Ensembl VEP 116.2    | `e0db557d89fef548d9c12df9e48637c5becd78d8de7c858c9ece8b5e31c578c5` |

Recorded for the Step 2 equivalence checks but not used by the frozen tables: utils (earlier build)
`64d63a83a17c2c8b1b335e7bd978d5b3a311ca2022260ae87ff3683d7a0c6f54`, MutEnricher 1.3.3 `75ff50a0a1eb2011496ae7495ad4f719977170b6fde85820deba32a8bebf7b4a`, motifbreakR 2.24.0 `104a28d3eb3bffb83c9224c500c0396e1cf76d8bde02c90efad646bac9abef4f`.
Package lists are in `images/<role>.packages.txt`.

## Runs

| Cohort                                                     | Donors | Run date                                                | Frozen consensus table SHA-256                                                                                                                     |
| ---------------------------------------------------------- | ------ | ------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------- |
| Liver-HCC-public                                           | 261    | 2026-10-06 (values identical to the 2026-09-01 results) | `f624a5ba1016f430eba7b0a6257b9534a9fb8c958e4d120e0fb05831dac2f36d`                                                                                 |
| LIRI-JP                                                    | 261    | 2026-10-06                                              | `a44a2ac207db7728323985c3258e088fe5da67b1f922b3cfcd5e8b81d6bb8f3e`                                                                                 |
| PACA-CA                                                    | 147    | 2026-10-06                                              | `e04820539b801ef8920fc751845a953223a800f0b9a8695c2430eafad3b6bb53`                                                                                 |
| SYNTH-NC (60 donors, 5 planted drivers + 1 spectrum decoy) | 60     | 2026-09-23                                              | `fa379a4fe79d963b1ed188e98dd9660949f4ae8d7d5e8cd31b2b7d46848d5802` (truth JSON `4c91d76f597697e35ac31903c3208a8e696ab071c92f3d25178d02ef22b26a16`) |

`baseline.sha256` holds the same checksums in `sha256sum -c` format, with paths relative to
`data/baseline/` in the project folder (`<cohort>/<cohort>.consensus.tsv`, `SYNTH-NC/SYNTH-NC.truth.json`).

SYNTH-NC at q < 0.1: ACAT recovers 5 of 5 planted drivers, EBM 4 of 5 (misses E5, 3 donors, designed
at the detection limit; ACAT q 0.073, EBM q 0.31). The spectrum decoy E6 is significant under both and,
with the A/T spectrum control removed, is no longer flagged.

## Combination check (2026-10-06 tables)

`combination_check.R` (ACAT 0.91 @ e737d7d, EmpiricalBrownsMethod 1.34.0, RobustRankAggreg 1.2.1)
on the three frozen tables; copies are kept with the baseline in `data/baseline/combination_check/`.

| File                                  | SHA-256                                                            |
| ------------------------------------- | ------------------------------------------------------------------ |
| `combination_check.R`                 | `fdbb21732b7d105e75c0c91a7084cc89d76fd4c1135aa5f47d54e955df905d9e` |
| `combination_method_comparison.csv`   | `50dd89859dc278a1dc8d12ce0cbca98818b4f76c5f6fef71d18307c5a0872530` |
| `combination_discordant_elements.csv` | `65590b2652c21aa32eeb71262e76895468d68313d58badc38d416af0dfda64e5` |

## Rebuilt images compared with the baseline images (Step 2.4, 2026-10-22)

Published images (ghcr.io/khalofx/noncodingdrivers-*) compared with `images/*.packages.txt` and on the committed fixtures:

| Image                 | Package differences                                                                                                                                     | Fixture output                                                                                                                                                                                                                                                                                    |
| --------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| oncodrivefml 2.5.0    | none (pip, setuptools, wheel are not listed by the newer `pip freeze`)                                                                                  | identical to the baseline image, md5 e766fa4e1a7a3b115b5ecd6a3de57bfd                                                                                                                                                                                                                             |
| mutenricher 1.3.3     | none (pip, setuptools, wheel as above)                                                                                                                  | not run yet (Step 3.5)                                                                                                                                                                                                                                                                            |
| activedriverwgs 1.2.1 | none: R 4.5.3, ActiveDriverWGS 1.2.1, BSgenome.Hsapiens.UCSC.hg38 1.4.5                                                                                 | counts identical; pp_element and fdr_element differ by at most 1.1e-13 (relative) from the site-R baseline runtime, because the image uses OpenBLAS 0.3.26 and site R OpenBLAS 0.3.32. Regression comparisons of ActiveDriverWGS p-values therefore use a relative tolerance, not md5             |
| utils 1.0.0           | pysam and pyliftover removed on purpose (no v1.0 user); python-dateutil, pytz and tzdata had floated, so utils 1.0.1 pins them to the baseline versions | no fixture (no module uses utils yet)                                                                                                                                                                                                                                                             |
| combinepvalues 1.0.0  | new image (the baseline combination ran in a local R environment)                                                                                       | combination_check.R in this image reproduces combination_method_comparison.csv and combination_discordant_elements.csv byte for byte; bin/combine_pvalues.R gives the same significant sets (EBM 31 / 36 / 1, ACAT 28 / 33 / 0), discordant sets and TERT q (EBM 3.07e-05, rank 5; ACAT 9.54e-06) |

## Annotation-layer baseline (Step 4, frozen 2026-10-26)

Cluster folder `data/baseline/step4_annotation/` (read-only; lines also in `data/baseline/baseline.sha256`):

| File                                             | Origin                                                                                                                                                                                                                                                    | SHA-256                                                            |
| ------------------------------------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------ |
| `Liver-HCC-public.candidate_mutations.tsv`       | release run, `variants/` (266 rows, 28 ACAT elements)                                                                                                                                                                                                     | `f7a2ea634dc733a89d9688849a30a3214998fd2f67faf07cb052abb821bc2311` |
| `Liver-HCC-public.candidate_mutations.vep.tsv`   | release run, VEP summary columns                                                                                                                                                                                                                          | `e0d645c29b9a8a233a63bc0313f437d3ecb602661c0fcce3807a2a95423cce55` |
| `Liver-HCC-public.vep.raw.tsv`                   | release run, VEP tab output: **ensembl-vep 116.2**, cache 116                                                                                                                                                                                             | `b51b7afd2ea5241d40cfb5106ab527404972109f902a16e5606e37c2f3315646` |
| `Liver-HCC-public.variants.versions.yml`         | release run                                                                                                                                                                                                                                               | `6efbc07dc67ae76bba52c107573a57b03de6ff0d7a301f122f3ffa87ddc4c4c6` |
| `Liver-HCC-public.motifbreakr.pilot_effects.tsv` | motifbreakR pilot, 10 Oct 2026 (motifbreakR 2.24.0, JASPAR 2026 CORE vertebrates, p < 1e-4, uniform background)                                                                                                                                           | `452534bc257ce400617ff732cf4cd5dcee611eace4694d1f8f74ee1bcfeb3936` |
| `SYNTH-NC.regulatory_context.v0.1.tsv`           | v0.1 `regulatory_context.py` on every catalogue element x 12 tissues (SYNTH-NC run; the same 49 ENCODE files, `refs/encode/PROVENANCE.tsv`). The release ran no regulatory-context layer, so this is the reference for Step 4 gate 1. It has no chrY rows | `e8a0c7fc19d256bd07ddf2f5913dc0de6c9f18bbd9816e5f8968f58033984fba` |

The pipeline runs ensembl-vep **116.3** (nf-core `ensemblvep/vep`); differences from 116.2 are listed in `docs/validation.md`.

motifbreakR image: `containers/motifbreakr/` (explicit conda lock, 285 packages). SIF built on Karakoram from `motifbreakr.def`
with `apptainer build --fakeroot` (Apptainer 1.5.3): `noncodingdrivers-motifbreakr-2.24.0.sif`, SHA-256
`73f3101aea2ce5820a9ee136344822a96975a937996fa55bd2b7da32ad85b455`; installed packages equal the lock. In it, the
universalmotif route reproduces the pilot matrices exactly with `convert_type(pseudocount = 3.2)` (0.8 per base; 0.8 differs by up to 0.126), and TERT C228T x ELK4 MA0076.3 gives ref p 9.994507e-04, alt p 2.288818e-05 with both parsers.
