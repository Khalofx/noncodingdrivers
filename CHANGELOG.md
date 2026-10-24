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

### `Deprecated`
