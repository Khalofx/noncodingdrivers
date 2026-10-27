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
