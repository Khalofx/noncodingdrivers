# khalofx/noncodingdrivers: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [unreleased<!-- TODO nf-core: replace with date on release -->]

Initial release of khalofx/noncodingdrivers, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Driver core: `mutations`-tier samplesheet (`assets/schema_input.json`), per-cohort table merge
  (nf-core `csvtk/concat` + `gawk`), ActiveDriverWGS, OncodriveFML (licence-gated:
  `--accept_oncodrivefml_licence`), harmonisation (`HARMONIZECALLS`), anchoring on `element_id`
  (nf-core `csvtk/join`) and p-value combination (`COMBINEPVALUES`: empirical Brown's method by
  default, ACAT, robust rank aggregation, votes, Benjamini–Hochberg), in the `DRIVER_DISCOVERY`
  subworkflow.
- MutEnricher as a QC-only caller (`--run_mutenricher`): its native region table is published
  unmodified and never enters the combination; it reads per-sample VCFs.
- Run-summary MultiQC sections (`RUNSUMMARYMQC`, `bin/run_summary_mqc.py`): elements tested,
  per-caller coverage and caller output rows, without candidate ranking.
- Element coordinates as fragment spans with `n_fragments` for catalogue ids split by the
  coding-sequence exclusion.
- `benchmarks/`: regression gate tools (`compare_to_baseline.py`, `validate_liver_hcc.py`,
  `table_to_vcfs.py`), Rheinbay et al. 2020 Table 5 concordance and COSMIC audits; pytest and ruff
  CI (`.github/workflows/pytest.yml`).
- Test profile on the companion repository `Khalofx/noncodingdrivers-test-datasets` (synthetic
  TERT-window cohort), Karakoram site profile.
- Annotation layers in the `EXPLAIN` subworkflow; none of them filters or re-ranks a candidate:
  patient-resolved candidate mutations on the fragments of the significant elements (nf-core
  `bedtools/intersect`), nearest gene from the element midpoint to GENCODE gene TSSs with a
  strand-aware signed distance (nf-core `bedtools/closest`; `--gencode_gtf` or `--gene_tss_bed`),
  ENCODE regulatory context per tissue and layer (`--run_regulatory_context`; nf-core
  `bedtools/coverage` and `bedtools/map`), Ensembl VEP on candidate mutations (`--run_vep`; nf-core
  `ensemblvep/vep`, offline cache, never `--regulatory`) and motifbreakR on candidate SNVs
  (`--run_motifbreakr`; local `MOTIFBREAKR`, JASPAR 2026 CORE vertebrates bundled).
- Pipeline-owned motifbreakR image `ghcr.io/khalofx/noncodingdrivers-motifbreakr:2.24.0`, built
  from an explicit conda lock (`containers/motifbreakr/`), with an Apptainer recipe for HPC builds.
- Optional cohort-level signature QC (`--run_sigprofiler`; nf-core `sigprofiler`).
- Reference datasheet `assets/references/GRCh38.yml` (schema `assets/schema_references.json`):
  the validated GRCh38 reference set with checksums, paths relative to `--reference_base`;
  explicit parameters take precedence. Asset recipes and checksums in `docs/references.md`.
- `vcf`, `sarek` and `maf` entry tiers through the `INPUT_HARMONISATION` subworkflow: contig
  renaming (nf-core `bcftools/annotate`), MAF split per tumour (gawk) and conversion with
  maf2vcf.pl (local `MAF2VCF`, vcf2maf 1.6.22), sorting, normalisation, PASS and primary-chromosome
  filtering and the canonical table (nf-core `bcftools/sort`, `bcftools/norm`, `bcftools/view`,
  `bcftools/query`, gawk). The per-sample PASS VCFs feed MutEnricher.
- Assembly check per samplesheet row (local `ASSEMBLYCHECK`, EBI vcf_assembly_checker 0.10.2):
  REF alleles against the FASTA of the declared build; below `--assembly_min_match` (0.99) the run
  stops with `ASSEMBLY_MISMATCH` before any caller.
- Opt-in GRCh37/hg19 liftover (`--allow_hg19_liftover`; nf-core `picard/liftovervcf`), rejected
  records published; GRCh37/hg19 rows without the opt-in stop with `ASSEMBLY_NOT_ENABLED`.

### `Fixed`

### `Dependencies`

| Tool                  | Version |
| --------------------- | ------- |
| ActiveDriverWGS       | 1.2.1   |
| OncodriveFML          | 2.5.0   |
| MutEnricher           | 1.3.3   |
| EmpiricalBrownsMethod | 1.34.0  |
| ACAT                  | 0.91    |
| RobustRankAggreg      | 1.2.1   |
| csvtk                 | 0.37.0  |
| gawk                  | 5.3.1   |
| MultiQC               | 1.35    |
| bedtools              | 2.31.1  |
| Ensembl VEP           | 116.3   |
| motifbreakR           | 2.24.0  |
| universalmotif        | 1.28.0  |
| SigProfilerExtractor  | 1.2.2   |

### `Deprecated`
