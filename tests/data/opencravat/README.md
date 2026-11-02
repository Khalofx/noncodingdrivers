# OpenCRAVAT fixtures

- `syn5.opencravat_input.txt`: five synthetic variants in OpenCRAVAT native input (chrom, pos,
  strand, ref, alt, sample): TERT promoter chr5:1295113 G>A and chr5:1295135 G>A, chr1:1000000
  T>A, a deletion and an insertion. Invented sample names.
- `candidates.tsv`: the same variants as a candidate table, plus chr5:1295160 C>T, which is not in
  the report (`not_in_report`).
- `syn5.opencravat.report.tsv`: the OpenCRAVAT 3.2.4 text report of `syn5.opencravat_input.txt`
  (`oc run -l hg38 -a encode_tfbs ccre_screen repeat clinvar -t text`, modules of
  `containers/opencravat/modules.freeze.json`, 2 November 2026). Its values come from the module
  data: GENCODE gene models, SCREEN cCREs and ENCODE TFBS (ENCODE: no restrictions), UCSC
  RepeatMasker (free for public and commercial use) and ClinVar (NCBI: no restrictions on use or
  distribution).
- `no_original_input.report.tsv`: a hand-written report without the Original Input columns (negative test).
