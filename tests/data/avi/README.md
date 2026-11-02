# AlphaGenome AVI fixtures

All scores here are invented: none is AlphaGenome output, so the files are licence-clean.

- `avi_scores.tsv.gz` (+ `.tbi`): five rows in the layout of the AlphaGenome "AVI SNV scores"
  download (bgzip TSV: chrom, pos, ref, alt, raw, phred; tabix `-s 1 -b 2 -e 2`):
  chr1:1000000 T>A, chr5:1295113 G>A/C/T and chr5:1295135 G>A.
- `no_index/avi_scores.tsv.gz`: the same file without its `.tbi` (negative parameter test).
- `candidates.tsv`: a candidate table with ten rows and eight distinct variants: four SNVs with a
  score, one SNV without (chr5:1295160 C>T), a deletion, an insertion and a GG>AA substitution.
- `vep_avi.tab.gz`, `vep_no_avi.tab.gz`: hand-written Ensembl VEP tab tables for those variants,
  with and without the plugin columns `AVI_PHRED` and `AVI_RAW`.
- `avi.tsv`: the `AVI_TABLE` output for `candidates.tsv` and `vep_avi.tab.gz` (input of the
  element-summary test).
- `make_avi_fixture.sh`: writes a synthetic score file over the SNVs of any candidate table
  (Karakoram integration run, `docs/validation.md`).
