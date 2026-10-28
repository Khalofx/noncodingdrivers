# khalofx/noncodingdrivers: Validation

Regression of the annotation layers (`EXPLAIN`) against the frozen baseline
(`tests/regression/PROVENANCE.md`, "Annotation-layer baseline"). Cohort: PCAWG Liver-HCC (261
donors, GRCh38), run on Karakoram with `--run_regulatory_context --run_vep --run_motifbreakr`; the
caller tasks were reused from the Step 3 regression runs (`-resume`).

## Continuity run (`--combination_method acat`, 28 significant elements)

| Check                                      | Result                                                                                                                                                                                                                                                                                                                               |
| ------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Candidate rows                             | 268: the release's 266 sample x mutation x element rows, plus 2 rows in the first fragment of `enhancer_ENSR00000371308` (chr1:144,994,401-144,995,007), which the release missed because it located candidates by the element's last fragment                                                                                       |
| TERT promoter (`promoter_ENSR00000746745`) | 11 patient rows, 11 donors                                                                                                                                                                                                                                                                                                           |
| Nearest gene                               | equal to the release for 266/266 release rows (28/28 elements); `distance_to_tss` changes only for the two fragmented elements, measured from the midpoint of the element span (TERT promoter 156 -> -6, `enhancer_ENSR00000371308` 34,775 -> 34,375)                                                                                |
| Regulatory context                         | 1,323/1,323 comparable element x file rows identical to the earlier implementation on the same 49 ENCODE files (peak count, covered bases, covered fraction, maximum signal); 621 unstaged tissue x layer slots absent, not 0. `CTCF_binding_site_ENSR00001294671` (chrY) has no earlier value: the earlier table holds no chrY rows |
| Ensembl VEP 116.2 -> 116.3                 | 171/171 release variants have identical transcript consequences (Feature, Consequence, IMPACT, SYMBOL, BIOTYPE, CANONICAL, DISTANCE); coding and splice flags on protein-coding transcripts identical for every release candidate. **No differences to list**                                                                        |
| motifbreakR                                | TERT C228T (chr5:1,295,113 G>A): ELK4 (MA0076.3) gain, reference p 9.99e-4, alternative p 2.29e-5, as in the 10 Oct 2026 pilot; 507 SNV x TF rows at 139 SNVs, 52 with a p below `calculatePvalue()` resolution (flagged, not floored)                                                                                               |
| MultiQC                                    | report renders with the per-caller run summary and the VEP summary                                                                                                                                                                                                                                                                   |

## Default run (`--combination_method ebm`, 31 significant elements)

| Check                                                                          | Result                                                                                                                                                                                                                                                                   |
| ------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Candidate rows                                                                 | 315 rows from 129 donors on 31 elements (4 fragmented): the candidate-row baseline of the default run                                                                                                                                                                    |
| Nearest gene                                                                   | one row per significant element, equal to the earlier rule on the element span for 31/31                                                                                                                                                                                 |
| Regulatory context                                                             | 1,519/1,519 element x file rows identical to the earlier implementation; 713 unstaged slots absent                                                                                                                                                                       |
| TP53-adjacent `CTCF_binding_site_ENSR00000548580` (significant under EBM only) | all 3 candidate mutations are splice-site variants of protein-coding TP53 transcripts: chr17:7,673,533 A>G and chr17:7,673,700 C>T (splice donor), chr17:7,673,609 C>A (splice acceptor). Its signal is a coding-gene splice effect, not evidence of a non-coding driver |
| motifbreakR                                                                    | TERT C228T ELK4 gain as above; 699 SNV x TF rows at 175 SNVs, 73 with a p below resolution                                                                                                                                                                               |
| MultiQC                                                                        | renders with the per-caller run summary and the VEP summary                                                                                                                                                                                                              |

## A cohort without a significant element

PACA-CA under `--combination_method acat` has no element with `is_consensus == True`. The run logs
`no element has is_consensus == True; candidate mutations, nearest genes and annotation layers are
skipped, no candidate file is written` and publishes no candidate, nearest-gene or annotation file.

## Input tiers and hg19 liftover (Liver-HCC-public, 261 donors)

The same cohort was run from three inputs: the canonical `mutations` table (the validated
baseline), 261 per-sample VCFs projected from that table (`benchmarks/table_to_vcfs.py`, indels
anchored on the FASTA), and the hg19 PCAWG MAF of the 261 aliquots (a `csvtk grep` subset of
`final_consensus_passonly.snv_mnv_indel.icgc.public.maf.gz`, 3,380,218 records) lifted with
`--allow_hg19_liftover`. The criterion is an equal SNV set; every other difference is listed.

| Check                                          | VCF rows                                                                                                                 | MAF row (hg19, Picard LiftoverVcf)                                                                                                                                                                                                        |
| ---------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Assembly check                                 | 261/261 rows pass; 3,379,234/3,379,234 REF alleles match                                                                 | 3,152,665/3,152,665 MAF SNVs match hg19                                                                                                                                                                                                   |
| SNV rows (3,151,863)                           | sorted SHA-256 equal to the `mutations` table                                                                            | sorted SHA-256 equal to the `mutations` table                                                                                                                                                                                             |
| DNP, TNP, ONP (32,377, 361, 19)                | identical                                                                                                                | identical                                                                                                                                                                                                                                 |
| Deletions (129,518 in the table)               | 129,443: 405 left-aligned by `bcftools norm`; 75 were listed twice in the table (below)                                  | 129,439: 402 left-aligned; the same 75 duplicates; 4 rejected by Picard                                                                                                                                                                   |
| Insertions (65,116 in the table)               | 65,095: 5 left-aligned; 1 listed twice; 20 on alt or unplaced contigs                                                    | 65,095: 154 left-aligned; 1 listed twice; 20 on alt or unplaced contigs                                                                                                                                                                   |
| Records outside the analysis contigs           | 20 insertions on 7 alt/unplaced contigs (absent from the 25-contig FASTA)                                                | the same 20 are not in the run: the chain maps them to contigs the FASTA lacks, which Picard rejects as `NoTarget`                                                                                                                        |
| Caller results (ActiveDriverWGS, OncodriveFML) | both native result tables byte-identical to the `mutations`-tier run                                                     | both native result tables byte-identical to the `mutations`-tier run                                                                                                                                                                      |
| Consensus (`--combination_method ebm`)         | 541,936 elements, identical p and q; 31 significant, the same set; Brown c 1.013805, df 3.94553; TERT q 3.07e-05, rank 5 | identical, as for the VCF rows                                                                                                                                                                                                            |
| Liftover rejects                               |                                                                                                                          | 900 (802 SNVs: 797 `NoTarget`, 5 `MismatchedRefAllele`; 11 MNVs; 87 indels incl. 2 `IndelStraddlesMultipleIntervals`); the earlier pyliftover conversion failed 875 (176 unmapped, 2 split, 697 REF mismatches) and kept the same SNV set |

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

## Variant flags, evidence table and null calibration (Step 6)

Commit `6e71b36` on Karakoram. SYNTH-NC (60 synthetic donors, one GRCh38 `maf` row) and
Liver-HCC-public with 3 null replicates each; LIRI-JP without. Step 6 adds files only: the
consensus of both PCAWG cohorts is unchanged.

| Gate                            | Criterion                                                                                         | Result                                                                                                                                                                                                                                    |
| ------------------------------- | ------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 1. Spike-in recovery (SYNTH-NC) | Planted drivers recovered at q < 0.1, at least as in `baseline.json` (EBM 4/5, ACAT 5/5)          | EBM 4/5 (E1-E4), ACAT 5/5; ranks E1 1, E2 2, E3 4, E4 5, E5 8 (EBM q 0.31, at the detection limit as in the baseline); the spectrum decoy E6 rank 3                                                                                       |
| 2. Null calibration (one-sided) | At most one element at q < 0.1 per null replicate and method; no null lambda above 1.2            | 0 elements in all 24 checks (2 cohorts x 3 replicates x ActiveDriverWGS, OncodriveFML, EBM, ACAT); highest null lambda 1.127 (OncodriveFML, SYNTH-NC)                                                                                     |
| 3. Liver-HCC-public unchanged   | Same consensus as the Step 4 run: EBM 31 and ACAT 28 significant, TERT promoter rank 5 under both | 541,936 elements, maximum relative difference of every p- and q-value 0; EBM 31, ACAT 28; TERT promoter rank 5 under both (ranked by p; two elements tie at its ACAT q)                                                                   |
| 4. SRGAP2B flagged (LIRI-JP)    | `enhancer_ENSR00000371308` keeps its rank and carries the segdup and/or low-mappability flag      | rank 4, as in the Step 4 run; 24 mutations in 22 patients, all in GIAB segmental duplications and low-mappability regions (both fractions 1, `flagged`); Umap uniqueness 0.46; callable fraction 0.498; fails the mappability post-filter |

**Variant flags.** No mutation is removed; the flags and fractions are added files.

| Cohort           | Rows      | Blacklist v2  | Segdup         | Low mappability | Umap mean | In gnomAD v4.1 | `germline_pop_common` | Germline NA |
| ---------------- | --------- | ------------- | -------------- | --------------- | --------- | -------------- | --------------------- | ----------- |
| SYNTH-NC         | 231,472   | 8,326 (3.6%)  | 12,231 (5.3%)  | 18,745 (8.1%)   | 0.964     | 21,525         | 1,498 (0.65%)         | 0           |
| Liver-HCC-public | 3,379,254 | 15,565 (0.5%) | 94,215 (2.8%)  | 80,549 (2.4%)   | 0.997     | 490,354        | 27,252 (0.81%)        | 20          |
| LIRI-JP          | 3,605,363 | 17,100 (0.5%) | 103,411 (2.9%) | 90,946 (2.5%)   | 0.997     | 602,363        | 39,568 (1.10%)        | 33          |

The germline NA rows are the insertions on alt and unplaced contigs that the 25-contig FASTA lacks
(Step 5). The planted TERT C228T (chr5:1,295,113 G>A) is a PASS gnomAD v4.1 genomes allele at
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
