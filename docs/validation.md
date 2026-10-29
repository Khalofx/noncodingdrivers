# khalofx/noncodingdrivers: Validation

Every result below compares a run of this pipeline with the frozen regression baseline: the
per-element result tables of three PCAWG cohorts (Liver-HCC-public, LIRI-JP and PACA-CA, frozen on
2026-10-06) and of the synthetic cohort SYNTH-NC (2026-09-23), with the inputs, images, tool
versions and checksums recorded in [`tests/regression/PROVENANCE.md`](../tests/regression/PROVENANCE.md),
the checksums in `sha256sum -c` form in `tests/regression/baseline.sha256`, and the expected values
in `tests/regression/baseline.json`. The cohort data are not redistributed (the reuse terms of the
PCAWG tables are unverified), so the regression runs are made on the Karakoram cluster, where the
baseline folders are kept read-only. All runs are GRCh38.

## Release runs (v1.0.0)

Pipeline commit `d171e8f` on Karakoram, with every container pinned by digest, ActiveDriverWGS in
its image instead of a site R installation, and the updated nf-core modules (bcftools 1.25).

| Run                                                   | Check                                                           | Result                                                                                                                                                                                                                                                                                                                                         |
| ----------------------------------------------------- | --------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `test_full` (SYNTH-NC, 60 donors, from test-datasets) | Planted drivers at q < 0.1 (`baseline.json`: EBM 4/5, ACAT 5/5) | EBM 4/5 (E1-E4), ACAT 5/5; ranks E1 1, E2 2, E3 4, E4 5, E5 8 (EBM q 0.31); spectrum decoy E6 rank 3                                                                                                                                                                                                                                           |
|                                                       | Against the SYNTH-NC run of the flags release                   | mutation table, OncodriveFML table, variant and element flags, callable territory, post-filters, candidates and nearest gene identical; evidence table identical apart from the 4 null-calibration rows of that run; ActiveDriverWGS p within 1.44e-12 and every consensus value within 3.15e-12 (relative)                                    |
|                                                       | Evidence table                                                  | 61,950 rows, 0 schema errors                                                                                                                                                                                                                                                                                                                   |
|                                                       | Cost                                                            | 44 min, 28.2 CPU-h, 403 tasks                                                                                                                                                                                                                                                                                                                  |
| Liver-HCC-public, 261 donors, annotation layers on    | Consensus against the annotation-layer run                      | 541,936 elements; EBM 31 and ACAT 28 significant, the same sets; TERT promoter rank 5 under both; `is_consensus` and `consensus_rank` identical; every q-value within 1.4e-14; ActiveDriverWGS p within 1.06e-10, EBM p within 7.0e-11, ACAT p within 1.01e-8 (below)                                                                          |
|                                                       | Annotation layers                                               | 315 candidate rows from 129 donors and 31 nearest-gene rows identical; regulatory context 1,519 rows identical; Ensembl VEP for 213/213 variants, consequence table identical (the TP53-adjacent CTCF site keeps its 3 splice variants); TERT C228T ELK4 (MA0076.3) gain; motifbreakR 699 SNV x TF rows at 175 SNVs identical; MultiQC renders |
|                                                       | Flags and evidence against the flags run                        | variant flags (3,379,254 rows), element flags, callable territory and post-filters identical; evidence table identical apart from the 4 null-calibration rows of that run                                                                                                                                                                      |
|                                                       | Cost                                                            | 4 h 20 min, 180.9 CPU-h, 191 tasks                                                                                                                                                                                                                                                                                                             |

**ACAT near p = 1.** The largest ACAT difference is at `enhancer_ENSR00001170012` (ACAT p
0.0571933325 against 0.0571933319, q 1 in both). Its two caller p-values are 0.9983
(ActiveDriverWGS) and 0.00159 (OncodriveFML), so their Cauchy terms nearly cancel and a change of
1e-12 in the ActiveDriverWGS p moves the combined p in the eighth digit. No significant set, rank
or q-value changes.

## Driver core: callers and combination

The cohorts were run from their canonical mutation tables with ActiveDriverWGS and OncodriveFML
(CADD v1.7) and both combinations. The checks of `benchmarks/compare_to_baseline.py` and the
regression gate script passed on all four runs (85 table checks).

| Check                                     | Liver-HCC-public (261 donors)                                                                                                                                                             | LIRI-JP (261 donors)                                                                                   | PACA-CA (147 donors)                                                                                                       |
| ----------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------ | -------------------------------------------------------------------------------------------------------------------------- |
| Cohort mutation table                     | SHA-256 equal to the baseline (3,379,254 rows)                                                                                                                                            | SHA-256 equal to the baseline                                                                          | SHA-256 equal to the baseline                                                                                              |
| Per-caller p-values                       | identical for all 541,936 tested elements (ActiveDriverWGS and OncodriveFML; OncodriveFML seed 123 reproduces exactly)                                                                    | identical for all 542,138 elements                                                                     | identical for all 541,422 elements                                                                                         |
| ACAT (`--combination_method acat`)        | 28 significant, the baseline set; TERT promoter rank 5, q 9.54468254e-06 (baseline 9.54468554e-06, relative 3.1e-7)                                                                       | 33 significant, the baseline set; minimum q 6.928e-10 (baseline 7.223e-10)                             | 0 significant (minimum q 0.361), as in the baseline                                                                        |
| EBM (`--combination_method ebm`, default) | 31 significant: 25 shared with ACAT, 6 EBM-only, 3 ACAT-only (both discordant sets equal to the baseline combination check); TERT promoter rank 5, q 3.07e-05; Brown c 1.01380, df 3.9455 | 36 significant (28 shared, 8 EBM-only, 5 ACAT-only; discordant sets equal); Brown c 1.01444, df 3.9431 | 1 significant (`open_chromatin_region_ENSR00000763710`); Brown c 1, df 4 (negative covariance -0.0299, so Fisher's method) |

A second, independent Liver-HCC-public run under ACAT gave byte-identical caller and matched
tables. The TERT promoter (`promoter_ENSR00000746745`, two fragments) and the WDR74 CTCF site
(`CTCF_binding_site_ENSR00000431853`), both PCAWG liver elements (Rheinbay et al. 2020), are
significant under both combinations in Liver-HCC-public. Every ACAT-only element has an OncodriveFML
p above 0.05, and every EBM-only element has both caller p-values below 0.05.

**Tolerances.** The baseline ActiveDriverWGS runs used a site R installation; the pipeline image
agrees with it to a maximum relative difference of 1.1e-13 (OpenBLAS 0.3.26 against 0.3.32;
mutation counts identical), so ActiveDriverWGS p-values are compared with a relative tolerance of
1e-9. The baseline's ACAT q-values came from an in-house implementation of the Cauchy combination;
the pipeline uses the ACAT package, so its q-values are compared within the per-cohort bound that
the combination check recorded before any pipeline run (maximum difference in log10: Liver-HCC-public
3.81e-6, LIRI-JP 0.0181, PACA-CA 4.87e-12). OncodriveFML needed no tolerance.

## Annotation layers

Regression of the annotation layers against the annotation-layer baseline
(`tests/regression/PROVENANCE.md`, "Annotation-layer baseline"). Cohort: Liver-HCC-public, run with
`--run_regulatory_context --run_vep --run_motifbreakr`; the caller tasks were reused from the
driver-core runs (`-resume`).

### Continuity run (`--combination_method acat`, 28 significant elements)

| Check                                      | Result                                                                                                                                                                                                                                                                                                                                    |
| ------------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Candidate rows                             | 268: the baseline's 266 sample x mutation x element rows, plus 2 rows in the first fragment of `enhancer_ENSR00000371308` (chr1:144,994,401-144,995,007), which the baseline table lacks because it located candidates by the element's last fragment only                                                                                |
| TERT promoter (`promoter_ENSR00000746745`) | 11 patient rows, 11 donors                                                                                                                                                                                                                                                                                                                |
| Nearest gene                               | equal to the baseline for 266/266 baseline rows (28/28 elements); `distance_to_tss` changes only for the two fragmented elements, measured from the midpoint of the element span (TERT promoter 156 -> -6, `enhancer_ENSR00000371308` 34,775 -> 34,375)                                                                                   |
| Regulatory context                         | 1,323/1,323 comparable element x file rows identical to the annotation-layer baseline on the same 49 ENCODE files (peak count, covered bases, covered fraction, maximum signal); 621 unstaged tissue x layer slots absent, not 0. `CTCF_binding_site_ENSR00001294671` (chrY) has no baseline value: the baseline table holds no chrY rows |
| Ensembl VEP 116.2 -> 116.3                 | 171/171 baseline variants have identical transcript consequences (Feature, Consequence, IMPACT, SYMBOL, BIOTYPE, CANONICAL, DISTANCE); coding and splice flags on protein-coding transcripts identical for every baseline candidate. **No differences to list**                                                                           |
| motifbreakR                                | TERT C228T (chr5:1,295,113 G>A): ELK4 (MA0076.3) gain, reference p 9.99e-4, alternative p 2.29e-5, as in the 10 Oct 2026 pilot; 507 SNV x TF rows at 139 SNVs, 52 with a p below `calculatePvalue()` resolution (flagged, not floored)                                                                                                    |
| MultiQC                                    | report renders with the per-caller run summary and the VEP summary                                                                                                                                                                                                                                                                        |

### Default run (`--combination_method ebm`, 31 significant elements)

| Check                                                                          | Result                                                                                                                                                                                                                                                                   |
| ------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Candidate rows                                                                 | 315 rows from 129 donors on 31 elements (4 fragmented): the candidate-row baseline of the default run                                                                                                                                                                    |
| Nearest gene                                                                   | one row per significant element, equal to the baseline rule (element span) for 31/31                                                                                                                                                                                     |
| Regulatory context                                                             | 1,519/1,519 element x file rows identical to the annotation-layer baseline; 713 unstaged slots absent                                                                                                                                                                    |
| TP53-adjacent `CTCF_binding_site_ENSR00000548580` (significant under EBM only) | all 3 candidate mutations are splice-site variants of protein-coding TP53 transcripts: chr17:7,673,533 A>G and chr17:7,673,700 C>T (splice donor), chr17:7,673,609 C>A (splice acceptor). Its signal is a coding-gene splice effect, not evidence of a non-coding driver |
| motifbreakR                                                                    | TERT C228T ELK4 gain as above; 699 SNV x TF rows at 175 SNVs, 73 with a p below resolution                                                                                                                                                                               |
| MultiQC                                                                        | renders with the per-caller run summary and the VEP summary                                                                                                                                                                                                              |

### A cohort without a significant element

PACA-CA under `--combination_method acat` has no element with `is_consensus == True`. The run logs
`no element has is_consensus == True; candidate mutations, nearest genes and annotation layers are
skipped, no candidate file is written` and publishes no candidate, nearest-gene or annotation file.

## Input tiers and hg19 liftover (Liver-HCC-public, 261 donors)

The same cohort was run from three inputs: the canonical `mutations` table (the validated
baseline), 261 per-sample VCFs projected from that table (`benchmarks/table_to_vcfs.py`, indels
anchored on the FASTA), and the hg19 PCAWG MAF of the 261 aliquots (a `csvtk grep` subset of
`final_consensus_passonly.snv_mnv_indel.icgc.public.maf.gz`, 3,380,218 records) lifted with
`--allow_hg19_liftover`. The criterion is an equal SNV set; every other difference is listed.

| Check                                          | VCF rows                                                                                                                 | MAF row (hg19, Picard LiftoverVcf)                                                                                                                                                                                                           |
| ---------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Assembly check                                 | 261/261 rows pass; 3,379,234/3,379,234 REF alleles match                                                                 | 3,152,665/3,152,665 MAF SNVs match hg19                                                                                                                                                                                                      |
| SNV rows (3,151,863)                           | sorted SHA-256 equal to the `mutations` table                                                                            | sorted SHA-256 equal to the `mutations` table                                                                                                                                                                                                |
| DNP, TNP, ONP (32,377, 361, 19)                | identical                                                                                                                | identical                                                                                                                                                                                                                                    |
| Deletions (129,518 in the table)               | 129,443: 405 left-aligned by `bcftools norm`; 75 were listed twice in the table (below)                                  | 129,439: 402 left-aligned; the same 75 duplicates; 4 rejected by Picard                                                                                                                                                                      |
| Insertions (65,116 in the table)               | 65,095: 5 left-aligned; 1 listed twice; 20 on alt or unplaced contigs                                                    | 65,095: 154 left-aligned; 1 listed twice; 20 on alt or unplaced contigs                                                                                                                                                                      |
| Records outside the analysis contigs           | 20 insertions on 7 alt/unplaced contigs (absent from the 25-contig FASTA)                                                | the same 20 are not in the run: the chain maps them to contigs the FASTA lacks, which Picard rejects as `NoTarget`                                                                                                                           |
| Caller results (ActiveDriverWGS, OncodriveFML) | both native result tables byte-identical to the `mutations`-tier run                                                     | both native result tables byte-identical to the `mutations`-tier run                                                                                                                                                                         |
| Consensus (`--combination_method ebm`)         | 541,936 elements, identical p and q; 31 significant, the same set; Brown c 1.013805, df 3.94553; TERT q 3.07e-05, rank 5 | identical, as for the VCF rows                                                                                                                                                                                                               |
| Liftover rejects                               |                                                                                                                          | 900 (802 SNVs: 797 `NoTarget`, 5 `MismatchedRefAllele`; 11 MNVs; 87 indels incl. 2 `IndelStraddlesMultipleIntervals`); the baseline's pyliftover conversion failed 875 (176 unmapped, 2 split, 697 REF mismatches) and kept the same SNV set |

**Duplicate representations.** The `mutations` table, converted from the PCAWG MAF without
normalisation, lists 75 deletions and 1 insertion twice for the same donor, in two equivalent
positions of a repeat (for example chr1:19,578,665 `AGAC>-` and chr1:19,578,666 `GACA>-` in a
`GAGACA` repeat). `bcftools norm` left-aligns both copies to one record (chr1:19,578,664
`GAGAC>G`), and the cohort merge keeps it once.

**Picard versus pyliftover.** Four deletions that pyliftover lifted (by their first base) are
rejected by Picard: two span a chain gap (`IndelStraddlesMultipleIntervals`) and two have their
anchor base in an unmapped interval (`NoTarget`).

**Wrong declared assembly.** The same MAF declared `GRCh38` stops at the assembly check with
`ASSEMBLY_MISMATCH` (818,924 of 3,152,665 REF alleles, 25.98%, match GRCh38); no caller task ran.

## Variant flags, evidence table and null calibration

Commit `6e71b36` on Karakoram. SYNTH-NC (60 synthetic donors, one GRCh38 `maf` row) and
Liver-HCC-public with 3 null replicates each; LIRI-JP without. The flags, evidence table and null calibration add files only: the
consensus of both PCAWG cohorts is unchanged.

| Gate                            | Criterion                                                                                                   | Result                                                                                                                                                                                                                                              |
| ------------------------------- | ----------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 1. Spike-in recovery (SYNTH-NC) | Planted drivers recovered at q < 0.1, at least as in `baseline.json` (EBM 4/5, ACAT 5/5)                    | EBM 4/5 (E1-E4), ACAT 5/5; ranks E1 1, E2 2, E3 4, E4 5, E5 8 (EBM q 0.31, at the detection limit as in the baseline); the spectrum decoy E6 rank 3                                                                                                 |
| 2. Null calibration (one-sided) | At most one element at q < 0.1 per null replicate and method; no null lambda above 1.2                      | 0 elements in all 24 checks (2 cohorts x 3 replicates x ActiveDriverWGS, OncodriveFML, EBM, ACAT); highest null lambda 1.127 (OncodriveFML, SYNTH-NC)                                                                                               |
| 3. Liver-HCC-public unchanged   | Same consensus as the annotation-layer run: EBM 31 and ACAT 28 significant, TERT promoter rank 5 under both | 541,936 elements, maximum relative difference of every p- and q-value 0; EBM 31, ACAT 28; TERT promoter rank 5 under both (ranked by p; two elements tie at its ACAT q)                                                                             |
| 4. SRGAP2B flagged (LIRI-JP)    | `enhancer_ENSR00000371308` keeps its rank and carries the segdup and/or low-mappability flag                | rank 4, as in the annotation-layer run; 24 mutations in 22 patients, all in GIAB segmental duplications and low-mappability regions (both fractions 1, `flagged`); Umap uniqueness 0.46; callable fraction 0.498; fails the mappability post-filter |

**Variant flags.** No mutation is removed; the flags and fractions are added files.

| Cohort           | Rows      | Blacklist v2  | Segdup         | Low mappability | Umap mean | In gnomAD v4.1 | `germline_pop_common` | Germline NA |
| ---------------- | --------- | ------------- | -------------- | --------------- | --------- | -------------- | --------------------- | ----------- |
| SYNTH-NC         | 231,472   | 8,326 (3.6%)  | 12,231 (5.3%)  | 18,745 (8.1%)   | 0.964     | 21,525         | 1,498 (0.65%)         | 0           |
| Liver-HCC-public | 3,379,254 | 15,565 (0.5%) | 94,215 (2.8%)  | 80,549 (2.4%)   | 0.997     | 490,354        | 27,252 (0.81%)        | 20          |
| LIRI-JP          | 3,605,363 | 17,100 (0.5%) | 103,411 (2.9%) | 90,946 (2.5%)   | 0.997     | 602,363        | 39,568 (1.10%)        | 33          |

The germline NA rows are the insertions on alt and unplaced contigs that the 25-contig FASTA lacks
(see input tiers above). The planted TERT C228T (chr5:1,295,113 G>A) is a PASS gnomAD v4.1 genomes allele at
fafmax_faf95_max 1.1e-4, below `--germline_faf_threshold` (1e-3): it is not flagged. Of the 544,737
catalogue elements, 10,899 have a callable fraction below 0.5 (the same in every cohort: the
fraction depends on the territory only).

**Post-filters of the significant elements** (`--postfilter_mode flag_only`). In both liver cohorts
the same two significant elements fail, both on mappability (no mutation outside GIAB low
mappability), and every other significant element passes all three filters:

| Element                    | Region                                 | Liver-HCC-public rank (q) | LIRI-JP rank (q) | Mutations / patients | Segdup, low mappability | Callable fraction |
| -------------------------- | -------------------------------------- | ------------------------- | ---------------- | -------------------- | ----------------------- | ----------------- |
| `enhancer_ENSR00000371308` | chr1:144,994,401-145,004,400 (SRGAP2B) | 6 (3.2e-5)                | 4 (4.2e-6)       | 22/22 and 24/22      | 1, 1                    | 0.498             |
| `promoter_ENSR00000982577` | chr14:19,511,000-19,511,400            | 25 (0.030)                | 35 (0.041)       | 4/4 and 4/4          | 1, 1                    | 0                 |

SYNTH-NC: all 6 significant elements pass.

**Null calibration** (`--run_null_calibration`, D32, D35). The null moves each SNV to a random
position with the same trinucleotide within +/- 50 kb, in the same patient (PCAWG Sanger design,
Rheinbay et al. 2020, Supplementary Information section 6); lambda = chi-squared quantile of the
median p / 0.4549 over the elements with an ActiveDriverWGS expected count >= 1.

| Cohort (elements)          | Method          | Null lambda, median [min, max] | Null median p | Observed lambda | Observed median p | Null status  |
| -------------------------- | --------------- | ------------------------------ | ------------- | --------------- | ----------------- | ------------ |
| SYNTH-NC (485)             | ActiveDriverWGS | 0.471 [0.372, 0.513]           | 0.643         | 0.435           | 0.657             | conservative |
|                            | OncodriveFML    | 1.096 [0.925, 1.127]           | 0.480         | 0.933           | 0.515             | within_band  |
|                            | EBM             | 0.651 [0.562, 0.689]           | 0.586         | 0.540           | 0.620             | conservative |
|                            | ACAT            | 0.601 [0.582, 0.686]           | 0.601         | 0.521           | 0.626             | conservative |
| Liver-HCC-public (108,429) | ActiveDriverWGS | 0.452 [0.451, 0.453]           | 0.650         | 0.435           | 0.656             | conservative |
|                            | OncodriveFML    | 0.957 [0.949, 0.958]           | 0.509         | 0.847           | 0.535             | within_band  |
|                            | EBM             | 0.522 [0.520, 0.527]           | 0.626         | 0.463           | 0.646             | conservative |
|                            | ACAT            | 0.579 [0.568, 0.579]           | 0.608         | 0.519           | 0.627             | conservative |

The gate is one-sided (D35): only inflation produces false drivers, and power is gate 1. The
median-based lambda of ActiveDriverWGS sits near 0.45 by construction: it is a one-sided count
test that gives p = 1 to elements without mutations and inverts the p-value of depleted elements,
so the null median p is about 0.65 instead of 0.5; EBM and ACAT inherit this. The observed runs
match their own nulls (observed / null lambda 0.83-0.96), and PCAWG likewise removed only inflated
methods and kept deflated ones (Supplementary Information section 7). `--null_lambda_band` keeps
the default `0.8,1.2`; its lower bound sets the `conservative` label only. Because the null
mirrors ActiveDriverWGS's own +/- 50 kb trinucleotide model, which its authors benchmarked, the
check mainly validates the EBM/ACAT consensus, the pipeline's caller set-up (fragmented catalogue,
sharded ActiveDriverWGS with FDR over all shards) and OncodriveFML. It cannot detect effects below
the 50 kb scale (CTCF-site hypermutation, APOBEC hairpins, mapping artefacts such as the two
elements above): those are what the variant flags and post-filters report. Cost: one observed run
per replicate (Liver-HCC-public with 3 replicates: 10 h, 1,773 CPU-h on Karakoram).

**Test profile.** The SYNTH-TERT fixture carries a chrM row, while the test FASTA covers chr5:1-1.4
Mb: the first test-profile run showed that the sites VCF for the gnomAD lookup must leave out
contigs absent from the FASTA index (bcftools sort stops on an undeclared contig). Fixed before
these runs; the germline columns of such rows are NA.

## Recorded differences from the baseline

- **Candidate rows** come from every fragment of a significant element; the baseline located them by
  the last fragment only (2 extra rows under ACAT, in `enhancer_ENSR00000371308`).
- **Nearest gene** distances are measured from the midpoint of the element span (changes only for
  fragmented elements).
- **Ensembl VEP** 116.3 instead of 116.2: no consequence differs.
- **Liftover** with Picard LiftoverVcf instead of pyliftover: the same SNV set; four deletions that
  pyliftover lifted are rejected (chain gap or unmapped anchor).
- **Duplicate indel representations** in a mutation table are left-aligned to one record by
  `bcftools norm` (75 deletions and 1 insertion in Liver-HCC-public).
- **ActiveDriverWGS** in its image instead of a site R install: relative 1.1e-13.
- **ACAT** from the package instead of an in-house implementation: within the recorded bounds.
- **ACAT near p = 1** amplifies the ActiveDriverWGS image difference to at most 1.01e-8 in an ACAT p-value (release runs), with no change of any q-value, set or rank.
- **Null calibration band**: `--null_lambda_band 0.8,1.2`, one-sided; a lambda below 0.8 is
  reported as `conservative`, never as a failure.
