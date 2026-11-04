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

## Variant annotation layers (Step 7A: AlphaGenome AVI, OpenCRAVAT)

SYNTH-NC (`-profile test_full,karakoram`), 2 November 2026, in one launch folder with `-resume`.
Run A is the `dev` commit b56908a with `--run_vep`; runs B to D are this branch. The AVI score file
is synthetic (`tests/data/avi/make_avi_fixture.sh` over the 42 candidate SNVs, invented scores):
the AlphaGenome download had not been staged.

| Run                                                   | Check                                                                                                             | Result                                                                                                                                                                                                     |
| ----------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| B: `--run_vep`, layers off                            | task cache against run A                                                                                          | 404 of 405 tasks cached, `ENSEMBLVEP_VEP` among them (same hash `01/ddaf8d`); only MultiQC re-ran (its methods file is rewritten every run). 143/143 published files byte-identical to run A               |
| D: all defaults                                       | published files against run A                                                                                     | byte-identical (every file of A except the two VEP files)                                                                                                                                                  |
| D: all defaults                                       | published files against the Step 7 run (`runs/SYNTH-NC-test_full`, 28 Oct), `pipeline_info/` and `multiqc/` aside | 141/141 files present; 24 byte-identical, 116 differ only in `##` header lines (bcftools command lines with their dates), `consensus/intermediate/SYNTH-NC.matched.tsv` has the same rows in another order |
| C: `--run_vep --run_alphagenome_avi --run_opencravat` | VEP table                                                                                                         | 1,342 rows as in run A; the 14 Step 7 columns identical, `AVI_PHRED` and `AVI_RAW` added and filled in all 1,342 rows                                                                                      |
| C                                                     | `candidates/avi/SYNTH-NC.avi.tsv`                                                                                 | 42 of 42 candidate variants (all SNVs) `scored`; every value equal to the score file                                                                                                                       |
| C                                                     | `candidates/avi/SYNTH-NC.avi.element_summary.tsv`                                                                 | 6 elements                                                                                                                                                                                                 |
| C                                                     | `opencravat/SYNTH-NC.opencravat.tsv`                                                                              | 42 of 42 variants `annotated`, 71 columns, no hg19 column; SCREEN class for 10 (PLS 2, dELS 4, CTCF-only 4), a repeat for 10, a ClinVar record for 2                                                       |
| C                                                     | `opencravat/SYNTH-NC.opencravat.network.tsv`, `.telemetry.tsv`                                                    | connection attempts only to 127.0.0.1:9 (3 connect, 3 lookup calls), 0 non-local targets; 1 metrics payload kept locally by OpenCRAVAT and deleted by the task                                             |
| C                                                     | other published files against run D                                                                               | byte-identical                                                                                                                                                                                             |

The Liver-HCC-public check with the real AVI file (at most 201 of its 213 candidate variants can be
scored; 12 are indels) is pending the user's download (to-do S7A.15).

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

## Opt-in callers: OncodriveCLUSTL and `--consensus_callers`

Branch `feat/new-callers` on Karakoram, SYNTH-NC `test_full`. Run A: `--run_oncodriveclustl true
--run_null_calibration true --null_replicates 3` (job 390177; 446 tasks, 104.5 CPU-h). Default-identity
and resume runs, all from one launch and work folder: D, dev at the merge-base `abcbe40` with
defaults (jobs 390406 and 390555, resumed after a scheduler stop); E, the branch (`430cb65`) with
defaults from its own checkout (job 390556, effectively fresh: 3 of 403 tasks reused); E2, the branch
checked out in D's worktree with defaults, `-resume` (job 390602); F2, the same with
`--run_oncodriveclustl true`, `-resume` (job 390603). Comparisons cover every published file except
`pipeline_info/` and the MultiQC report, plots, log, JSON, parquet and source lists (146 files).

**Default identity.** Two fresh runs of the same code already differ in bytes. They differ in the
`##` header lines of the bcftools outputs (command and date), in the tabix index of a VCF whose header
changed length, and in the row order of `consensus/intermediate/<cohort>.matched.tsv` (csvtk join). So
the reference is the fresh dev run D, and two comparisons are made:

| Comparison                                        | Byte-identical | Different | Different after normalisation                                           |
| ------------------------------------------------- | -------------- | --------- | ----------------------------------------------------------------------- |
| E2 (branch, defaults, `-resume` D) vs D           | 146 of 146     | 0         | 0                                                                       |
| E (branch, defaults, fresh) vs D                  | 30 of 146      | 116       | 55, all `harmonisation/vcf/*.pass.vcf.gz.tbi`; their VCFs are identical |
| D (dev, fresh) vs the Step 7 results (dev, fresh) | 32 of 146      | 114       | 53, all `*.pass.vcf.gz.tbi`                                             |

Normalisation drops the `##` lines (60 VCFs) and sorts the `matched.tsv` rows (1 file). The consensus
table is byte-identical in D, E, E2 and the Step 7 run (md5 `6f9ef619a74556c51e48be57b393fd2f`).
In E2 the 13 tasks that ran again are COMBINEPVALUES (its script and command changed: `--callers`) and
the 12 tasks downstream of its consensus table. Their outputs are byte-identical to D's.

**Task hashes.** Nextflow 25.10 names every `eval()` output `nxf_out_eval_<n>` from one counter
over the pipeline's process definitions, and the name enters the task hash
(nextflow-io/nextflow#7572). A new module with an `eval()` output therefore renames the outputs of
later modules and `-resume` re-runs them. The first version of this branch had three nf-core
`gawk`/`bedtools intersect` helper processes for the OncodriveCLUSTL prefilter. Its run A matched D's
`eval` output names in 4 of 34 processes, and a default run resumed from A reused 1 of the 387 tasks it started.
ONCODRIVECLUSTL now emits its version with `env()` and prefilters inside its task
(`bin/oncodriveclustl_prefilter.py`). The branch's default run matches D's `eval` names in 33 of 33
processes, and E2 reuses 390 of 403 tasks. Enabling OncodriveCLUSTL on a resumed default run (F2)
reuses 389 of 405 tasks: ActiveDriverWGS, OncodriveFML and every task upstream stay cached. The 16
that run are ONCODRIVECLUSTL, its HARMONIZECALLS, ANCHOR_JOIN, COMBINEPVALUES and the 12 consensus
consumers. A resume also needs the same checkout path: staged `projectDir` assets are hashed by path,
size and modification time, so E, from another worktree of the same commit, reused 3 tasks.

**Prefilter.** OncodriveCLUSTL analyses only elements with at least `-emut` SNVs, and its run time
grows with the size of the regions file, so the task passes only the rows of those ids. On SYNTH-NC
and its three null replicates, `bin/oncodriveclustl_prefilter.py` writes the same bytes as the
earlier `gawk` + `bedtools intersect -c` + `gawk` path. It keeps 3,088, 3,195, 3,227 and 3,219 of
544,737 ids (4,700-4,938 of 595,725 rows) in 2 s with 0.34 GB maximum RSS. F2 (prefilter in the
task) and A (helper processes) give the same OncodriveCLUSTL tables and the same consensus.

**SYNTH-NC: OncodriveCLUSTL reported alongside.** Cutting the two carried columns
(`p_oncodriveclustl`, `q_oncodriveclustl`) from the consensus of A or F2 leaves the default consensus
(md5 `6f9ef619a74556c51e48be57b393fd2f`; EBM parameters identical). OncodriveCLUSTL analysed 3,088
elements and tested 36; an element is tested only if it has an observed cluster. Two of the 36 are at
q < 0.1. Its tables are identical to the design prototype's (md5 `d4cd8ad7c4f4f71386680b0c68e244ee`,
clusters `64dd8852938a1c7339013bd834c53f60`).

| Planted (design)                                   | SNVs | Clusters (clustered SNVs) | OncodriveCLUSTL p | q      | Rank of 36 tested | Consensus rank (EBM q) |
| -------------------------------------------------- | ---- | ------------------------- | ----------------- | ------ | ----------------- | ---------------------- |
| E1 `promoter_ENSR00000746745` (TERT hotspots)      | 18   | 2 (18)                    | 3.75e-5           | 0.0014 | 1                 | 1 (6.5e-32)            |
| E2 `enhancer_ENSR00001140199` (top-CADD positions) | 14   | 3 (13)                    | 0.0055            | 0.099  | 2                 | 2 (1.7e-28)            |
| E3 `enhancer_ENSR00000839370` (random positions)   | 10   | 0                         | not tested        |        |                   | 4 (1.5e-15)            |
| E4 `enhancer_ENSR00001144769` (4 top-CADD SNVs)    | 4    | 1 (4)                     | 0.23              | 0.33   | 23                | 5 (7.8e-9)             |
| E5 `enhancer_ENSR00001173236` (random positions)   | 3    | 0                         | not tested        |        |                   | 8 (0.31)               |
| E6 `enhancer_ENSR00001022568` (spectrum decoy)     | 12   | 1 (2)                     | 0.76              | 0.76   | 36                | 3 (3.3e-24)            |

OncodriveCLUSTL recovers the two planted elements whose SNVs cluster (E1, E2) and ranks the spectrum
decoy E6 last of the 36 tested. The consensus ranks E6 third. Burden at random positions (E3, E5) is
not a clustering signal. OncodriveCLUSTL complements the combined callers here and does not replace
them.

**Null calibration (OncodriveCLUSTL).** The median-p `lambda` of CALIBRATECALLERS is not defined for
OncodriveCLUSTL. Within the ActiveDriverWGS element set (expected count >= 1), OncodriveCLUSTL has a
p-value for 0, 1, 1 and 0 elements (observed run, nulls 1-3). Its null summary (`lambda` 2.32 from
two single elements, `excluded_for_inflation`) is therefore no verdict. The statistic used instead is
the null false-positive rate over the elements OncodriveCLUSTL analysed:

| Run      | Analysed | Tested | Minimum p | p <= 0.05 | p <= 0.01 | q < 0.1 |
| -------- | -------- | ------ | --------- | --------- | --------- | ------- |
| observed | 3,088    | 36     | 3.75e-5   | 3         | 2         | 2       |
| null 1   | 3,195    | 30     | 0.014     | 3         | 0         | 0       |
| null 2   | 3,227    | 30     | 0.053     | 0         | 0         | 0       |
| null 3   | 3,219    | 34     | 0.057     | 0         | 0         | 0       |

No null element reaches q < 0.1; the gate allows at most one per replicate. Over all analysed null
elements, 3 of 9,641 are at p <= 0.05, with untested elements counted as p = 1. Over the 94 tested null
elements it is 3 (3.2%; binomial P(X >= 3) = 0.85 at a 5% rate), and none is at p <= 0.01.
OncodriveCLUSTL is not inflated on SYNTH-NC. This is one synthetic cohort: combining it
(`--consensus_callers`) still prints the opt-in warning until it passes calibration on real cohorts.

**Shards.** A used one shard (3,088 ids; `--oncodriveclustl_shard_size` 5,000). With seven shards of 500
ids (job 390252) the following are unchanged:

- the observed counts and scores (TOTAL_MUT, CLUSTERED_MUT, CLUSTERS, SCORE, LENGTH);
- the 36 tested elements and the 39 cluster coordinates;
- the two elements at q < 0.1, with q-values that are Benjamini-Hochberg over all shards.

The simulated p-values change with the layout, because each shard draws its own simulations from
`--seed`: 35 of 36 differ, by at most 0.48 in log10 (E1: 3.75e-5 with one shard, 1.25e-5 with seven).

**Run time.** ONCODRIVECLUSTL on SYNTH-NC takes 48.5 s, and each null replicate 44-54 s. The largest
single process uses 1.49 GB (`/usr/bin/time`, job 389973). Nextflow's process-tree `peak_rss` for the same
run was 1.3 GB in A and 9.3 GB in F2: it adds up the forked workers, whose shared pages count once per
worker, at its sampling times. `process_medium` allows 36 GB. For comparison, OncodriveFML takes 17.6
min and ActiveDriverWGS 26 min on the same cohort.

**OncodriveFML coordinates (S7A.14).** OncodriveFML 2.5.0 reads the START of its regions as 1-based when
it maps mutations, but queries its CADD background with 0-based half-open coordinates. With the
adapter's START = BED start (unchanged), the scored positions are exactly the BED interval. One extra 5'
base per fragment is mapped, but it has no score and is dropped (43 of 231,472 SYNTH-NC SNVs). START + 1
would drop the first base of every fragment from the background. The adapter is kept.

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
