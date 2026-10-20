# Regression baseline: provenance

The regression gates compare v1 runs with the frozen tables described here. Everything is
identified by content (checksums), not by repository or tag. The frozen tables are kept on the
project cluster under `noncodingdrivers/data/baseline/` (read-only); there is no public deposit.

## Inputs

| Input | Checksum |
|---|---|
| PCAWG public somatic MAF (hg19) | MD5 `009b9f88ec10d07cf826174c3dea5b78` |
| Liver-HCC-public canonical mutation table (GRCh38; 3,379,254 rows) | SHA-256 `39bf4a0a3fe322a78abfb1e6b0a5b8d45800ecea7ae0caf1223c14111798f0cf` |
| LIRI-JP canonical mutation table (GRCh38) | SHA-256 `4000e4fe22e2c9642bc158d2d0ead0e554927097ba785f76591538111557742d` |
| PACA-CA canonical mutation table (GRCh38) | SHA-256 `5da1f75690a79b04193d9958f9f3100d406eb46555aecc63fdac24c4ff93612d` |
| Element catalogue: Ensembl Regulatory Build 110, CDS removed against GENCODE v50 basic (595,725 rows, 544,737 ids) | MD5 `39f838e6c81fdb0b152504f2d51939b5` |
| Regulatory Build 110 GFF (20221007) | MD5 `ee3b05c9020f7da4244d563995a579bb` (recorded earlier; not re-checked in Step 0) |
| GRCh38 FASTA, chr-named, 25 contigs (`.fai` MD5 `963c862ff26b0d914956faaf7192e54b`) | MD5 `838cd9f56ca7d7dd345cb025d188c991` |
| GENCODE v50 basic GTF | MD5 `11e77cf1d4b78b43f57d608e3f7f2202` |
| JASPAR 2026 CORE vertebrates (1,019 PFMs) | MD5 `c63218acf88eb5fe10255b7c8ef23d48` |
| hg19ToHg38 liftover chain | MD5 `35887f73fe5e2231656504d1f6430900` |
| CADD GRCh38 v1.7 whole-genome SNVs, `whole_genome_SNVs.tsv.gz` (user-supplied, not redistributed) | MD5 `88577a55f1cd519d44e0f415ba248eb9`; `.tbi` MD5 `347df8fac17ea374c4598f4f44c7ce8b` |
| Ensembl VEP 116 cache, `homo_sapiens_vep_116_GRCh38.tar.gz` | MD5 `6e10959391069a9eebb8d8d165fd29be` |
| OncodriveFML genome: official bgdata `genomereference/hg38-20161209` | copied unchanged into `data/refs/oncodrivefml_genome/` |

hg19 MAF to GRCh38 liftover of Liver-HCC-public (the baseline canonical table): 3,380,218 records in,
3,379,254 out; 875 failed (176 unmapped, 2 split mappings, 697 REF mismatches); 65,116 insertions
carried without REF validation. Variant-type labels in the table: SNP, DEL, INS, DNP, TNP, ONP.

## Tools

ActiveDriverWGS 1.2.1 (R 4.5.3 site library, no container; BSgenome.Hsapiens.UCSC.hg38 1.4.5);
OncodriveFML 2.5.0 (seed 123); Nextflow 25.10.4; python 3.11, pandas 2.3.3, numpy 1.26.4, scipy 1.17.1;
Ensembl VEP 116.2; MultiQC 1.21. MutEnricher and motifbreakR did not contribute to the frozen tables.

Images used by the baseline runs (from the task scripts of the 2026-10-06 runs), by SIF SHA-256:

| Role | SHA-256 |
|---|---|
| utils (python tools) | `ff62ef3b67216a5a1fc4acc7705409a5225cf215f24e02b38981c12632eb677b` |
| OncodriveFML 2.5.0 | `fd60e15863c41fc75c614f5616ca3a5711174719fe7d5c43026f003512661a5d` |
| Ensembl VEP 116.2 | `e0db557d89fef548d9c12df9e48637c5becd78d8de7c858c9ece8b5e31c578c5` |

Recorded for the Step 2 equivalence checks but not used by the frozen tables: utils (earlier build)
`64d63a83a17c2c8b1b335e7bd978d5b3a311ca2022260ae87ff3683d7a0c6f54`, MutEnricher 1.3.3 `75ff50a0a1eb2011496ae7495ad4f719977170b6fde85820deba32a8bebf7b4a`, motifbreakR 2.24.0 `104a28d3eb3bffb83c9224c500c0396e1cf76d8bde02c90efad646bac9abef4f`.
Package lists are in `images/<role>.packages.txt`.

## Runs

| Cohort | Donors | Run date | Frozen consensus table SHA-256 |
|---|---|---|---|
| Liver-HCC-public | 261 | 2026-10-06 (values identical to the 2026-09-01 results) | `f624a5ba1016f430eba7b0a6257b9534a9fb8c958e4d120e0fb05831dac2f36d` |
| LIRI-JP | 261 | 2026-10-06 | `a44a2ac207db7728323985c3258e088fe5da67b1f922b3cfcd5e8b81d6bb8f3e` |
| PACA-CA | 147 | 2026-10-06 | `e04820539b801ef8920fc751845a953223a800f0b9a8695c2430eafad3b6bb53` |
| SYNTH-NC (60 donors, 5 planted drivers + 1 spectrum decoy) | 60 | 2026-09-23 | `fa379a4fe79d963b1ed188e98dd9660949f4ae8d7d5e8cd31b2b7d46848d5802` (truth JSON `4c91d76f597697e35ac31903c3208a8e696ab071c92f3d25178d02ef22b26a16`) |

`baseline.sha256` holds the same checksums in `sha256sum -c` format, with paths relative to
`data/baseline/` in the project folder (`<cohort>/<cohort>.consensus.tsv`, `SYNTH-NC/SYNTH-NC.truth.json`).

SYNTH-NC at q < 0.1: ACAT recovers 5 of 5 planted drivers, EBM 4 of 5 (misses E5, 3 donors, designed
at the detection limit; ACAT q 0.073, EBM q 0.31). The spectrum decoy E6 is significant under both and,
with the A/T spectrum control removed, is no longer flagged.

## Combination check (2026-10-06 tables)

`combination_check.R` (ACAT 0.91 @ e737d7d, EmpiricalBrownsMethod 1.34.0, RobustRankAggreg 1.2.1)
on the three frozen tables; copies are kept with the baseline in `data/baseline/combination_check/`.

| File | SHA-256 |
|---|---|
| `combination_check.R` | `fdbb21732b7d105e75c0c91a7084cc89d76fd4c1135aa5f47d54e955df905d9e` |
| `combination_method_comparison.csv` | `50dd89859dc278a1dc8d12ce0cbca98818b4f76c5f6fef71d18307c5a0872530` |
| `combination_discordant_elements.csv` | `65590b2652c21aa32eeb71262e76895468d68313d58badc38d416af0dfda64e5` |
