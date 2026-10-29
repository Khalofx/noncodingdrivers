# khalofx/noncodingdrivers: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0 - [unreleased]

Initial release of khalofx/noncodingdrivers, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Entry tiers `vcf`, `sarek` and `maf` (`assets/schema_input.json`) through the
  `INPUT_HARMONISATION` subworkflow: contig renaming (nf-core `bcftools/annotate`), MAF split per
  tumour and conversion with maf2vcf.pl (local `MAF2VCF`, vcf2maf 1.6.22), sorting, normalisation,
  PASS and primary-chromosome filtering and the canonical cohort table (nf-core `bcftools/sort`,
  `bcftools/norm`, `bcftools/view`, `bcftools/query`, gawk); per-sample PASS VCFs published. The
  canonical `mutations` tier is kept for the test profile and the regression gates only.
- Assembly check of every row against the FASTA of its declared build (local `ASSEMBLYCHECK`, EBI
  vcf_assembly_checker 0.10.2): below `--assembly_min_match` (0.99) the run stops with
  `ASSEMBLY_MISMATCH` before any caller.
- Opt-in GRCh37/hg19 liftover (`--allow_hg19_liftover`; nf-core `picard/liftovervcf`) with the
  rejected records published; hg19 rows without the opt-in stop with `ASSEMBLY_NOT_ENABLED`.
- Reference datasheet `assets/references/GRCh38.yml` (schema `assets/schema_references.json`): the
  validated GRCh38 reference set with checksums, paths relative to `--reference_base`; explicit
  parameters take precedence. Asset recipes in `docs/references.md`.
- Driver discovery (`DRIVER_DISCOVERY`) on a CDS-excluded Ensembl Regulatory Build catalogue:
  ActiveDriverWGS 1.2.1 (sharded, FDR over all shards) and OncodriveFML 2.5.0 with CADD v1.7
  (licence-gated: `--accept_oncodrivefml_licence`); MutEnricher 1.3.3 as a QC-only caller
  (`--run_mutenricher`). Element coordinates as fragment spans with `n_fragments`.
- P-value combination (`COMBINEPVALUES`): the empirical Brown's method (default) or ACAT sets
  `consensus_q` (Benjamini–Hochberg, `--fdr_threshold`); ACAT, robust rank aggregation and caller
  votes are published as sensitivity columns; `--skip_acat` where the ACAT package is unavailable.
- Variant flags (`VARIANT_FLAGS`, on by default; flag, never filter): ENCODE Blacklist v2, GIAB v3.6
  segmental duplications and low mappability, Umap k100 uniqueness (nf-core `bedtools/intersect`,
  `bedtools/coverage`, `bedtools/subtract`) and gnomAD v4.1 genomes germline flags (nf-core
  `echtvar/anno`; archive built with `echtvar/encode`); element-level fractions and callable
  fraction.
- PCAWG-style post-filter flags (`POSTFILTERFLAGS`: minimum support, mappability, maximum sample
  fraction; `--postfilter_mode pcawg` adds `q_postfilter`) and the long-form evidence table
  (`WRITEDRIVEREVIDENCE`, `evidence/<cohort>.driver_evidence.tsv`, layout
  `assets/schema_driver_evidence.json`). Evidence tiers are assigned downstream from this table.
- Null calibration (`NULL_CALIBRATION`, `--run_null_calibration`, off by default): neutral copies
  of each cohort (each SNV moved within +/- 50 kb to the same trinucleotide, same patient) re-run
  through the callers; inflation factor, null discoveries and QQ points per caller and
  combination, with a one-sided `--null_lambda_band` gate. Reported only.
- Annotation layers of the significant elements (`EXPLAIN`; none filters or re-ranks): candidate
  mutations per patient on the element fragments (nf-core `bedtools/intersect`), nearest GENCODE
  gene TSS with a strand-aware signed distance (nf-core `bedtools/closest`), ENCODE regulatory
  context per tissue and layer (`--run_regulatory_context`; nf-core `bedtools/coverage`,
  `bedtools/map`), Ensembl VEP 116.3 on candidate mutations (`--run_vep`; nf-core
  `ensemblvep/vep`, offline cache) and motifbreakR 2.24.0 on candidate SNVs (`--run_motifbreakr`;
  JASPAR 2026 CORE vertebrates bundled).
- Optional cohort-level signature QC (`--run_sigprofiler`; nf-core `sigprofiler`).
- MultiQC 1.35 report with a per-cohort run summary (`RUNSUMMARYMQC`).
- Pipeline-owned images `ghcr.io/khalofx/noncodingdrivers-{utils,activedriverwgs,mutenricher,combinepvalues,motifbreakr}`,
  built from `containers/` by `.github/workflows/build-containers.yml` and pinned by digest; the
  OncodriveFML image is built locally from `containers/oncodrivefml/`.
- Profiles `test` (synthetic TERT-window cohort), `test_full` (SYNTH-NC, synthetic genome-wide
  cohort with planted drivers) on `Khalofx/noncodingdrivers-test-datasets`, and the `karakoram`
  site profile.
- Regression baseline records (`tests/regression/`) and gate tools (`benchmarks/`); pytest and ruff
  CI.

### `Dependencies`

| Tool                  | Version |
| --------------------- | ------- |
| ActiveDriverWGS       | 1.2.1   |
| OncodriveFML          | 2.5.0   |
| MutEnricher           | 1.3.3   |
| EmpiricalBrownsMethod | 1.34.0  |
| ACAT                  | 0.91    |
| RobustRankAggreg      | 1.2.1   |
| bcftools              | 1.25    |
| bedtools              | 2.31.1  |
| csvtk                 | 0.37.0  |
| echtvar               | 0.2.4   |
| gawk                  | 5.3.1   |
| Picard                | 3.5.0   |
| vcf2maf               | 1.6.22  |
| vcf-validator         | 0.10.2  |
| Ensembl VEP           | 116.3   |
| motifbreakR           | 2.24.0  |
| universalmotif        | 1.28.0  |
| SigProfilerExtractor  | 1.2.2   |
| MultiQC               | 1.35    |
