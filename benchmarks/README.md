# benchmarks/

Gate tooling and validation scripts. Nextflow never runs anything in this folder; the
regression gates of the build guide run these by hand or from `local/` on the cluster.

| Script                                        | Purpose                                                                                                                                                   |
| --------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `compare_to_baseline.py`                      | Compare a v1 consensus table with `tests/regression/` and the frozen combination-check CSVs (significant sets, ranks, per-caller p, EBM/ACAT discordance) |
| `validate_liver_hcc.py`                       | Liver-HCC-public acceptance checks (`--layer core`: significant set, loci, TERT; `--layer full`: explanation-layer outputs)                               |
| `table_to_vcfs.py`                            | Project a canonical mutation table into per-sample PASS VCFs (fixtures and the VCF-parity gate)                                                           |
| `rheinbay_concordance.py`                     | Rheinbay et al. 2020 (doi:10.1038/s41586-020-1965-x) Table 5 recovery                                                                                     |
| `cosmic_compare.py`, `cosmic_sample_audit.py` | COSMIC overlap audits (user-supplied COSMIC files)                                                                                                        |
| `make_cohort_fixture.sh`                      | Provenance of the `tests/data/cohort/` fixtures                                                                                                           |

Tests for these scripts live in `benchmarks/tests/` and run in `.github/workflows/pytest.yml`.
