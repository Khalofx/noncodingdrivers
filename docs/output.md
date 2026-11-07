# khalofx/noncodingdrivers: Output

## Introduction

This document describes the output produced by the pipeline. `<cohort>` is the samplesheet `cohort` of a row; every table is tab-separated with a header row.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Input harmonisation](#input-harmonisation) - Assembly checks, liftover rejects, per-sample PASS VCFs and the cohort mutation table
- [Driver discovery](#driver-discovery) - Native result tables of ActiveDriverWGS, OncodriveFML and MutEnricher
- [Consensus](#consensus) - Per-element combined p- and q-values and the consensus call
- [Variant flags](#variant-flags) - Mask and germline flags of every mutation row, element-level fractions and callable fraction
- [Post-filter flags and evidence table](#post-filter-flags-and-evidence-table) - PCAWG-style post-filter flags and the long-form `driver_evidence.tsv`
- [Null calibration](#null-calibration) - Per-caller inflation factor on neutral copies of each cohort (opt-in)
- [Candidates and annotation layers](#candidates-and-annotation-layers) - Mutations, nearest gene, regulatory context, VEP (with AlphaGenome AVI scores), motifbreakR and OpenCRAVAT for the significant elements
- [Context layers](#context-layers) - COSMIC matches, cancer gene-list labels, ENCODE-rE2G enhancer-gene links and SCREEN V4 cCREs (each opt-in)
- [Signature QC](#signature-qc) - SigProfiler on the cohort (opt-in)
- [MultiQC](#multiqc) - Run summary report
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

### Driver discovery

<details markdown="1">
<summary>Output files</summary>

- `callers/activedriverwgs/`
  - `<cohort>.activedriverwgs.tsv`: the ActiveDriverWGS result table of all shards (`id`, `pp_element`, `element_muts_obs`, `element_muts_exp`, `element_enriched`, the site-level columns, `fdr_element`, `fdr_site`, `has_site_mutations`); `fdr_element` is recomputed over all shards.
  - `<cohort>.activedriverwgs.untested.txt`: element ids ActiveDriverWGS could not test (for example an empty background window). They are recorded, not dropped, and cannot be significant.
- `callers/oncodrivefml/`
  - `<cohort>.oncodrivefml.tsv`: the native OncodriveFML table (`GENE_ID` is the element id; `P_VALUE`, `Q_VALUE` and their `_NEG` counterparts, mutation and sample counts).
  - `<cohort>.oncodrivefml.status.txt`: the exit code of the OncodriveFML run.
- `callers/oncodriveclustl/` (`--run_oncodriveclustl`)
  - `<cohort>.oncodriveclustl.tsv`: the native OncodriveCLUSTL element table (`ENSID` is the element id; `TOTAL_MUT`, `CLUSTERED_MUT`, `CLUSTERS`, `SIM_CLUSTERS`, `SCORE`, `P_`/`Q_EMPIRICAL`, `P_`/`Q_ANALYTICAL`, `P_`/`Q_TOPCLUSTER`; `CGC` is `Non Available` on hg38). One row per element with at least `--oncodriveclustl_element_mutations` SNVs; p-values are empty for elements without an observed cluster (not tested). q-values are Benjamini–Hochberg over the tested elements of all shards. `P_ANALYTICAL`/`Q_ANALYTICAL` are harmonised as `p_oncodriveclustl`/`q_oncodriveclustl`.
  - `<cohort>.oncodriveclustl.clusters.tsv`: one row per cluster (`RANK` of its element, coordinates, `MAX_COORD`, `WIDTH`, `N_MUT`, `N_SAMPLES`, `FRA_UNIQ_SAMPLES`, `SCORE`, `P`).
  - `<cohort>.oncodriveclustl.log`: the OncodriveCLUSTL log of every shard.
- `callers/mutenricher/` (`--run_mutenricher`): the native MutEnricher region table (SNVs only, sample-count statistic, local background). QC only: it never enters the combination.

</details>

### Consensus

<details markdown="1">
<summary>Output files</summary>

- `consensus/`
  - `<cohort>.consensus.tsv`: one row per element tested by at least one combined caller (`--consensus_callers`). `p_<caller>`, `q_<caller>`: per-caller p- and q-values, of the combined callers and of every caller reported alongside (carried, not combined); `chrom`, `start`, `end`, `n_fragments`: the element span (a catalogue id split by the coding-sequence exclusion spans several fragments); `n_tested_callers` (combined callers that tested the element); `ebm_p`, `ebm_q` (empirical Brown's method) and `acat_p`, `acat_q` (ACAT), each with Benjamini–Hochberg q-values; `rra_score`, `rra_q` (robust rank aggregation) and `vote_count` (combined callers with their own q <= `--fdr_threshold`) as sensitivity columns; `combination_method`, `consensus_q` (the q-value of the chosen combination), `consensus_rank` and `is_consensus` (`True` when `consensus_q` < `--fdr_threshold`).
  - `<cohort>.consensus.ebm_parameters.tsv`: the callers combined, the number of elements and of complete elements, the sum of the off-diagonal covariance, and Brown's scale `brown_c` and degrees of freedom `brown_df`.
- `consensus/intermediate/`
  - `elements.coords.tsv`: span and fragment count of every catalogue id.
  - `<cohort>.<caller>.harmonized.tsv`: each caller's result as `element_id`, `p_<caller>`, `q_<caller>`.
  - `<cohort>.matched.tsv`: the per-caller columns joined on `element_id`, the input of the combination.

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
  - `<cohort>.driver_evidence.tsv`: long form, one row per cohort x element x evidence item, with `value`, `unit`, `status`, `threshold`, `source`, `source_version` and `native_or_lifted` (layout: `assets/schema_driver_evidence.json`; sources and versions: `assets/evidence_sources.tsv`). Items: the consensus and per-caller p- and q-values (combined callers and callers reported alongside; with `--run_oncodriveclustl` also `oncodriveclustl_n_clusters` and `oncodriveclustl_clustered_mutations` from its native table, `not_tested` for elements below its mutation cutoff), the element flags (mask fractions are `flagged` above 0.5), the post-filter flags and, with `--run_null_calibration`, one cohort-level row per method (`element_id` `.`, `null_lambda_<method>`; for OncodriveCLUSTL `null_fpr_oncodriveclustl` instead, the share of its tested null elements at p <= 0.05 over the replicates, threshold `0.05`). When a caller is reported alongside, one cohort-level row per caller records its role (`caller_role_<caller>`: value `1`, status `combined`, or value `0`, status `reported_alongside`); with every caller combined (the default) these rows are not written. `--evidence_scope nominal` (default) writes the elements with any p < 0.01 or `consensus_q` < 0.25; `all` writes every element.
  - Layer summaries: for each written element with `is_consensus` `True` and each annotation or context layer that is on, a fixed set of element-level items read from that layer's own table (no per-variant or per-link rows; written after the element's other items). A significant element missing from a layer's table gets that layer's items with value `NA`, status `not_applied`. With every layer off (the default) no such row is written and the table is unchanged. None of these items filters or re-ranks an element.

    | Layer (option)                            | Items                                                                                                                                                                                                                                                                                          | Read from                                          |
    | ----------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------- |
    | AlphaGenome AVI (`--run_alphagenome_avi`) | `avi_max_phred` (unit `score`; `NA`, `not_applied` when no candidate variant is scored), `avi_n_scored`, `avi_n_unscored` (candidate variants with and without an AVI score)                                                                                                                   | `candidates/avi/<cohort>.avi.element_summary.tsv`  |
    | OpenCRAVAT (`--run_opencravat`)           | `opencravat_n_ccre_screen`, `opencravat_n_repeat`, `opencravat_n_clinvar`: candidate variants of the element with at least one non-NA column of that annotator                                                                                                                                 | `opencravat/<cohort>.opencravat.tsv`               |
    | COSMIC (`--run_cosmic_context`)           | counts only: `cosmic_n_matched_mutations` (distinct candidate mutations with an exact COSMIC match), `cosmic_n_samples` (COSMIC samples carrying one), `cosmic_n_sample_name_overlaps`; `source_version` `COSMIC <release>`; no COSMIC identifier or record                                    | `cosmic/<cohort>.cosmic.element_summary.tsv`       |
    | Cancer gene list (`--cancer_gene_list`)   | `cancer_gene_nearest_in_list`, `cancer_gene_nearest_pc_in_list` (value `1`, status `in_list`, or `0`, `not_in_list`; `NA`, `not_applied` without a gene), `cancer_gene_re2g_targets_in_list` (count; `NA` without ENCODE-rE2G); `source_version` names the list and its MD5                    | `cancer_genes/<cohort>.cancer_gene_context.tsv`    |
    | ENCODE-rE2G (`--run_enhancer_gene_links`) | `re2g_n_distal_targets` (distinct distal target genes over the requested tissues), `re2g_top_distal_target` (unit `label`: the gene with the highest distal score over those tissues; status `top_is_nearest`, `nearest_in_targets`, `nearest_not_in_targets`, `no_link` or `no_nearest_gene`) | `enhancer_gene/<cohort>.enhancer_gene_summary.tsv` |
    | SCREEN V4 (`--run_ccre_overlay`)          | `ccre_v4_dominant_class` (unit `label`; `none` without a cCRE), `ccre_v4_covered_fraction`                                                                                                                                                                                                     | `ccre/<cohort>.ccre_v4_summary.tsv`                |

    Items computed from the candidate mutations (AVI, OpenCRAVAT, COSMIC) carry the cohort's coordinate provenance; the element-level overlaps (cancer gene list, ENCODE-rE2G, SCREEN V4) are `native`. Only the two label items carry a name in `value`; the schema keeps every other value numeric or `NA`.

</details>

### Null calibration

`--run_null_calibration` (off by default; in-house, reported only). Each replicate re-runs the
callers on a neutral copy of the cohort and costs about as much as the observed run. Read the
number of elements at q < `--fdr_threshold` in the null runs first: any element there is a false
positive under the null. A method is `excluded_for_inflation` when a null replicate's `lambda` is
above the upper bound of `--null_lambda_band`; a `lambda` below the lower bound is reported as
`conservative` (less power, no false positives). ActiveDriverWGS, and the EBM and ACAT combinations
that include it, have a null `lambda` near 0.5 by construction (one-sided count test; see
`docs/validation.md`). OncodriveCLUSTL (`--run_oncodriveclustl`) tests only elements with an
observed cluster, a data-dependent subset, so the median-p `lambda` is not defined for it: it has no
row in `null_calibration.tsv` and no QQ points. It gets a false-positive-rate check per null
replicate instead (`null_fpr.tsv`, below).

<details markdown="1">
<summary>Output files</summary>

- `null_calibration/`
  - `<cohort>.null_calibration.tsv`: per run (`observed`, `null1` ... `null<n>`) and method (each combined caller, `ebm`, `acat`): the element set (ActiveDriverWGS expected count >= 1), median p, inflation factor `lambda` (median chi-squared quantile / 0.4549) and the number of elements with q < `--fdr_threshold`; the `null_summary` rows give the median, minimum and maximum null `lambda` and the band status (`within_band`, `excluded_for_inflation` above `--null_lambda_band`, `conservative` below it).
  - `<cohort>.null_qq.tsv`: expected and observed -log10 p per run and method (all points with observed -log10 p >= 2 and up to 500 evenly spaced ranks).
  - `<cohort>.null_fpr.tsv` (with `--run_oncodriveclustl`): OncodriveCLUSTL's false-positive-rate check per run (`observed`, `null1` ... `null<n>`): `n_analysed` (elements in its native table, at least `--oncodriveclustl_element_mutations` SNVs), `n_tested` (elements with a p-value), `min_p`, `n_p_le_0.05`, `n_p_le_0.01`, `n_q_lt_fdr` and `fpr_p_le_0.05` (= `n_p_le_0.05` / `n_tested`). The `null_summary` row sums the counts over the replicates (pooled `fpr_p_le_0.05`), keeps the smallest p and the largest per-replicate `n_q_lt_fdr`, and is `pass` when no null replicate has more than one element at q < `--fdr_threshold`, otherwise `fail`.
  - `simulation/<cohort>.null<r>.null_summary.tsv`: rows carried, non-SNV rows not carried, SNVs dropped (REF mismatch, undefined context) and duplicates collapsed per replicate.

</details>

### Candidates and annotation layers

Written only for cohorts with at least one element with `is_consensus` true. None of these files
filters or re-ranks an element.

<details markdown="1">
<summary>Output files</summary>

- `candidates/`
  - `<cohort>.candidates.tsv`: one row per patient x mutation x significant element: `sample chrom pos ref alt vtype`, `element_id` and the fragment the mutation lies in (`fragment_start`, `fragment_end`). Mutations between the fragments of an element (excluded coding bases) are not candidates.
  - `<cohort>.nearest_gene.tsv`: one row per significant element: the nearest GENCODE gene TSS to the element midpoint (`nearest_gene`, `nearest_gene_ensembl_id`, `gene_type`, `nearest_tss`, `gene_strand`), the strand-aware signed `distance_to_tss` and its absolute value, `relative_position` (`overlapping_tss`, `upstream` or `downstream`) and `nearest_gene_status`.
  - `intermediate/<cohort>.significant.sorted.bed`: the fragments of the significant elements.
  - `vep/<cohort>.candidates.vep.tab.gz` and `_summary.html` (`--run_vep`): Ensembl VEP 116.3 consequences of the candidate mutations (tab format, offline cache; `--regulatory` is not used). With `--run_alphagenome_avi` the table has two more columns, `AVI_PHRED` and `AVI_RAW` (VEP AVI plugin; `-` when not scored).
  - `avi/` (`--run_alphagenome_avi`): AlphaGenome Variant Impact scores, annotation and hypothesis only.
    - `<cohort>.avi.tsv`: one row per distinct candidate variant: `variant_id` (`chrom:pos:ref:alt`), `chrom`, `pos`, `ref`, `alt`, `vtype`, `element_ids` (comma-separated), `n_samples`, `avi_phred`, `avi_raw` (NA when not scored) and `avi_status`: `scored`, `not_scored_non_snv` (insertion, deletion or multi-base substitution: the scores cover SNVs only), `not_in_avi_file` (an SNV the score file does not hold) or `not_in_vep_output`.
    - `<cohort>.avi.element_summary.tsv`: one row per element with candidate variants: `n_candidate_variants`, `n_snv`, `n_scored`, `max_avi_phred` with its `top_variant` (NA when nothing is scored) and `n_phred_ge_20` (`AVI_PHRED` 20 = top 1 % of all possible substitutions; a description, not a threshold).
- `regulatory_context/` (`--run_regulatory_context`)
  - `<cohort>.regulatory_context.tsv`: one row per significant element x staged ENCODE file: `tissue`, `layer`, `accession`, `n_fragments`, `n_peaks`, `covered_bp`, `element_bp`, `covered_frac`, `max_signal` with its `signal_source` (narrowPeak `signalValue` or DNase hotspot `max_density`; comparable within a layer only) and `status` (`overlap` or `no_overlap`). A tissue x layer that was not staged has no row.
  - `per_file/`: the bedtools coverage and signal BEDs behind each row.
- `motifbreakr/` (`--run_motifbreakr`)
  - `<cohort>.motifbreakr.tsv`: the native motifbreakR result for the candidate SNVs.
  - `<cohort>.motifbreakr.per_tf.tsv`: one row per SNV x transcription factor: the best motif, `ref_pvalue`, `alt_pvalue`, `best_pvalue`, `direction` (`gain`, `loss` or `none`), `fold_change`, and `p_below_resolution` for p-values below the resolution of `calculatePvalue()` (flagged, never floored; their `fold_change` is NA).
- `opencravat/` (`--run_opencravat`): OpenCRAVAT 3.2.4 on the candidate mutations, annotation only.
  - `<cohort>.opencravat.tsv`: one row per distinct candidate variant: `variant_id`, `chrom`, `pos`, `ref`, `alt`, `vtype`, `element_ids`, `n_samples`, `oc_status` (`annotated` or `not_in_report`), then the variant-level OpenCRAVAT columns named `<module>__<column>`: `base__*` (OpenCRAVAT's normalised hg38 position and alleles, gene, transcript, sequence ontology, all mappings, sample count), `ccre_screen__*` (SCREEN cCRE accession and `classification`, e.g. PLS, pELS, dELS, CTCF-only), `clinvar__*` (ClinVar significance, review status, disease names, oncogenicity and somatic classifications) and `repeat__*` (RepeatMasker class, family and name). NA = no annotation. No hg19 coordinate is produced.
  - `<cohort>.opencravat.network.tsv`: every connection attempt and host lookup of the run (`call`, `host`, `port`, `n_calls`, `is_local`); a non-local target fails the run.
  - `<cohort>.opencravat.telemetry.tsv`: the run's telemetry record: the configured metrics, store and publish addresses (`http://127.0.0.1:9`), the proxy, the number of non-local targets (0) and the number of metrics payloads OpenCRAVAT saved locally, which the task deleted.
- `references/`
  - `gencode.v50.basic.annotation.tss.bed`: the gene TSSs used for the nearest-gene table.

</details>

### Context layers

Each layer is off by default and writes only its own folder, for cohorts with at least one
significant element. None of them filters, re-ranks or weights an element; the tables have stable
column names, one row per element (or per element and tissue) for the summaries, and each layer's
summary is also written to `driver_evidence.tsv` (see [the evidence table](#post-filter-flags-and-evidence-table)).

<details markdown="1">
<summary>Output files</summary>

- `cosmic/` (`--run_cosmic_context`): counts only, no COSMIC identifiers; cite COSMIC and the release when reporting them.
  - `<cohort>.cosmic.mutation_summary.tsv`: one row per distinct candidate variant: `cohort cosmic_release chrom pos ref alt vtypes element_ids n_elements n_cohort_samples cohort_samples`, `cosmic_match` (`yes`/`no`), `n_cosmic_mutation_ids`, `n_cosmic_samples`, `n_confirmed_somatic_samples` (status "Confirmed somatic variant"), `n_cosmic_studies`, `n_pubmed_pmids`, `n_exact_sample_name_overlaps` and `independence_status` (`no_match`, `possible_self_overlap` when a COSMIC sample name or id equals a cohort sample id exactly, otherwise `not_assessed`).
  - `<cohort>.cosmic.element_summary.tsv`: one row per significant element: `cohort cosmic_release element_id consensus_rank chrom start end consensus_q`, `n_candidate_sample_events`, `n_distinct_candidate_mutations`, `n_cosmic_matched_mutations`, `n_confirmed_somatic_mutations`, `n_cosmic_samples`, `n_exact_sample_name_overlaps`, `cosmic_match`.
  - `<cohort>.cosmic.provenance.json`: release, assembly and name of the COSMIC file, its computed and expected MD5 (`md5_check`: `match` or `not_requested`), the id column and line endings found, rows scanned, exact rows, collapsed observations, the candidate and match totals and the matching rules.
  - `restricted/<cohort>.cosmic.exact_matches.tsv`: **COSMIC records** - one row per collapsed COSMIC observation of a candidate variant (transcript rows collapsed by variant, COSV, COSS or sample name, phenotype, study, PMID and somatic status) with its COSMIC identifiers, sample name, study, PMID, status, HGVSG, genes and transcripts, `source_transcript_rows`, `exact_sample_name_overlap` and `independence_status`.
  - `restricted/NOTICE.txt`: the COSMIC licence notice. Do not redistribute `restricted/` (see [usage](usage.md#cosmic)).
- `cancer_genes/` (`--cancer_gene_list`)
  - `<cohort>.cancer_gene_context.tsv`: one row per significant element: `nearest_gene nearest_gene_ensembl_id nearest_gene_type` and its `nearest_in_list` (`yes`, `no`, `NA` without a nearest gene), `nearest_match_by` (`ensembl_id` or `symbol`), `nearest_list_tier`, `nearest_list_role`; the same for the nearest protein-coding gene (`nearest_protein_coding_gene`, `_ensembl_id`, `nearest_protein_coding_distance` in bp, `nearest_pc_in_list`, `nearest_pc_match_by`, `nearest_pc_list_tier`, `nearest_pc_list_role`); `re2g_targets_in_list` and `re2g_distal_targets_in_list` (ENCODE-rE2G targets in the list; `NA` when that layer is off); `list_name`, `list_md5`, `list_genes`.
- `enhancer_gene/` (`--run_enhancer_gene_links`)
  - `<cohort>.enhancer_gene_links.tsv`: one row per element x rE2G file x rE2G element x target gene, sorted by element, tissue and score: `element_id tissue biosample accession re2g_element re2g_class` (`promoter`, `genic`, `intergenic`) `target_gene target_gene_ensembl_id target_gene_tss is_self_promoter re2g_score overlap_bp` (bases of the element's fragments covered by the rE2G element).
  - `<cohort>.enhancer_gene_summary.tsv`: one row per significant element x requested tissue: `nearest_gene`, `nearest_gene_ensembl_id`, `n_files` staged and `n_files_with_link`, `n_link_rows`, `self_promoter_genes`, `n_distal_targets` and `distal_targets_by_max_score` (`GENE:score`, maximum over the tissue's files), `top_distal_target`/`top_distal_score`, `n_targets`, `top_target_any`/`top_any_score` (self-promoter links included), `nearest_gene_status` on the distal targets and `nearest_gene_status_any` on all targets (`top_is_nearest`, `nearest_in_targets`, `nearest_not_in_targets`, `no_link`, or `no_nearest_gene`), and `threshold` (`NA`: the files are pre-thresholded) with `threshold_source`.
- `ccre/` (`--run_ccre_overlay`)
  - `<cohort>.ccre_v4.tsv`: one row per element x SCREEN V4 cCRE: `element_id ccre_accession ccre_class overlap_bp element_bp`.
  - `<cohort>.ccre_v4_summary.tsv`: one row per significant element: `element_bp n_ccre covered_bp covered_frac dominant_class` (class with most bases; `none` without a cCRE) and `ccre_v4_classes_bp` (`class:bases`, most bases first).
  - `<cohort>.ccre_v4_celltype.tsv` (`--ccre_celltype_tissues`): one row per element x biosample file x cCRE: `element_id tissue biosample accession ccre_accession celltype_class overlap_bp`.

</details>

### Signature QC

<details markdown="1">
<summary>Output files</summary>

- `qc/sigprofiler/` (`--run_sigprofiler`): the SigProfiler output (mutation matrices and extracted signatures) of the cohort's SNVs; QC only, nothing else reads it.

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

The report holds the run summary of each cohort (elements tested, per-caller coverage and output rows, without candidate ranking; when a caller is reported alongside, a caller-roles table labels each caller combined or reported alongside), the VEP summary when `--run_vep` is set, and the software versions of the run. For more information about how to use MultiQC reports, see <http://multiqc.info>.

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
