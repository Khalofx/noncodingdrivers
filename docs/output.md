# khalofx/noncodingdrivers: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Input harmonisation](#input-harmonisation) - Assembly checks, liftover rejects and per-sample PASS VCFs of vcf, sarek and maf rows
- [Variant flags](#variant-flags) - Mask and germline flags of every mutation row, element-level fractions and callable fraction
- [Post-filter flags and evidence table](#post-filter-flags-and-evidence-table) - PCAWG-style post-filter flags and the long-form `driver_evidence.tsv`
- [Null calibration](#null-calibration) - Per-caller inflation factor on neutral copies of each cohort (opt-in)
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### Input harmonisation

<details markdown="1">
<summary>Output files</summary>

- `harmonisation/assembly_check/`
  - `<row>.assembly_check.tsv`: one row per samplesheet row: records checked against the FASTA of the declared assembly, records matching, match fraction, records on contigs outside that FASTA, threshold and status (`pass`, `low_evidence`; a row with status `ASSEMBLY_MISMATCH` stops the run).
  - `<row>.assembly_mismatches.txt`: the first 1,000 mismatching records reported by vcf_assembly_checker.
- `harmonisation/liftover_rejects/`
  - `<sample>.unlifted.vcf.gz`: records of GRCh37/hg19 rows that Picard LiftoverVcf could not lift (FILTER names the reason, e.g. `NoTarget`, `MismatchedRefAllele`).
- `harmonisation/maf2vcf_skipped/`
  - `<sample>.skipped.tsv`: MAF rows that maf2vcf.pl dropped because their reference allele disagrees with the FASTA (only written when there are any).
- `harmonisation/vcf/`
  - `<sample>.pass.vcf.gz` (+ `.tbi`): per-sample VCF on GRCh38 after normalisation, PASS (or missing FILTER) records on chr1-22, X and Y; the input of MutEnricher.
- `mutations/`
  - `<cohort>.mutations.tsv`: the cohort mutation table (`sample chrom pos ref alt vtype`) that every caller reads: the mutations rows as given, plus the harmonised vcf, sarek and maf rows.

</details>

### Variant flags

Flag, never filter: the mutation table, the callers and the consensus are unchanged
(`--run_variant_flags`, on by default).

<details markdown="1">
<summary>Output files</summary>

- `variant_flags/`
  - `<cohort>.variant_flags.tsv`: one row per row of the cohort mutation table, in the same order. `in_blacklist_v2`, `in_segdup`, `in_low_mappability` (`true`/`false`): the variant's span overlaps the ENCODE Blacklist v2, the GIAB v3.6 segmental duplications or the GIAB v3.6 low-mappability regions (an insertion is tested on the base before it). `umap_unique_fraction`: covered fraction of the span by Umap k100 uniquely mappable sequence. `gnomad_filter`, `gnomad_af`, `gnomad_fafmax_faf95_max`: the gnomAD v4.1 genomes record of the same allele (indels anchored on the GRCh38 FASTA; several FILTER values comma-joined; NA when gnomAD does not carry the variant). `germline_pop_common`: `true` when the gnomAD FILTER is PASS and `fafmax_faf95_max` (grpmax FAF95) is above `--germline_faf_threshold`; NA when the germline flag is skipped.
  - `<cohort>.element_flags.tsv`: one row per element id of the catalogue. `n_mutations`, `n_patients`, `max_sample_fraction` (rows assigned by the point rule of the candidate table); `frac_in_blacklist_v2`, `frac_in_segdup`, `frac_in_low_mappability`, `mean_umap_unique_fraction`, `frac_germline_pop_common` over the element's mutations; `element_bp`, `callable_bp` and `callable_fraction`: the element's territory in Umap k100 unique sequence minus the ENCODE Blacklist v2.
- `references/`
  - `callable_territory.bed`: Umap k100 unique sequence minus the ENCODE Blacklist v2.

</details>

### Post-filter flags and evidence table

`--run_driver_evidence` (on by default). Evidence tiers are not assigned by the pipeline: they are
derived downstream from `driver_evidence.tsv`.

<details markdown="1">
<summary>Output files</summary>

- `postfilter/`
  - `<cohort>.postfilter.tsv`: one row per element of the consensus table with a value, a threshold and `pass`, `fail` or `not_applied` for minimum support (`--postfilter_min_mutations`, `--postfilter_min_patients`), mappability (share of the element's mutations outside GIAB low mappability above `--postfilter_min_mappable_fraction`) and maximum sample fraction (`--postfilter_max_sample_fraction`), and `postfilter_status`. With `--postfilter_mode pcawg`, `q_postfilter`: Benjamini-Hochberg after failing elements are set to p = 1; `consensus_q` is never changed.
- `evidence/`
  - `<cohort>.driver_evidence.tsv`: long form, one row per cohort x element x evidence item, with `value`, `unit`, `status`, `threshold`, `source`, `source_version` and `native_or_lifted` (layout: `assets/schema_driver_evidence.json`; sources and versions: `assets/evidence_sources.tsv`). Items: the consensus and per-caller p- and q-values, the element flags (mask fractions are `flagged` above 0.5), the post-filter flags and, with `--run_null_calibration`, one cohort-level row per method (`element_id` `.`, `null_lambda_<method>`). `--evidence_scope nominal` (default) writes the elements with any p < 0.01 or `consensus_q` < 0.25; `all` writes every element.

</details>

### Null calibration

`--run_null_calibration` (off by default; in-house, reported only). Each replicate re-runs the
callers on a neutral copy of the cohort and costs about as much as the observed run. Read the
number of elements at q < `--fdr_threshold` in the null runs first: any element there is a false
positive under the null. A method is `excluded_for_inflation` when a null replicate's `lambda` is
above the upper bound of `--null_lambda_band`; a `lambda` below the lower bound is reported as
`conservative` (less power, no false positives). ActiveDriverWGS, and the EBM and ACAT combinations
that include it, have a null `lambda` near 0.5 by construction (one-sided count test; see
`docs/validation.md`).

<details markdown="1">
<summary>Output files</summary>

- `null_calibration/`
  - `<cohort>.null_calibration.tsv`: per run (`observed`, `null1` ... `null<n>`) and method (each combined caller, `ebm`, `acat`): the element set (ActiveDriverWGS expected count >= 1), median p, inflation factor `lambda` (median chi-squared quantile / 0.4549) and the number of elements with q < `--fdr_threshold`; the `null_summary` rows give the median, minimum and maximum null `lambda` and the band status (`within_band`, `excluded_for_inflation` above `--null_lambda_band`, `conservative` below it).
  - `<cohort>.null_qq.tsv`: expected and observed -log10 p per run and method (all points with observed -log10 p >= 2 and up to 500 evenly spaced ranks).
  - `simulation/<cohort>.null<r>.null_summary.tsv`: rows carried, non-SNV rows not carried, SNVs dropped (REF mismatch, undefined context) and duplicates collapsed per replicate.

</details>

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://docs.seqera.io/platform-cloud/reports/overview) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
