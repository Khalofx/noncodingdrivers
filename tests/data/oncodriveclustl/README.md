# OncodriveCLUSTL fixtures

- `elements.bed`: the synthetic TERT-window elements (`tests/data/activedriverwgs/elements.bed`) plus
  `elem_split` (two fragments with one SNV each of `tests/data/activedriverwgs/mutations.tsv`, so two
  SNVs as one element) and `elem_empty` (no SNV).
- `elements.prefiltered.bed`: the rows of the ids with at least 2 SNVs (`elem_B`, `TERT_promoter_fixture`,
  `elem_D`, `elem_split`), in catalogue order: what `CLUSTL_ELEMENT_COUNTS` + `CLUSTL_REGIONS` keep.
  `elem_A`, `elem_C` and `elem_E` hold one SNV, `elem_empty` none, `elem_unplaced` is off the primary
  chromosomes.
- `expected.*`, `shard*.*` (unit test of `bin/merge_oncodriveclustl_shards.py`): the 36 tested elements and
  39 clusters of the synthetic SYNTH-NC cohort, as OncodriveCLUSTL 1.2.0 wrote them in one run (`--seed 1`),
  and the same rows split into two shuffled shards with their q-values and ranks overwritten.
