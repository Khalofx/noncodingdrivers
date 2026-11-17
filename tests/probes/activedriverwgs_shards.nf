//
// Test probe for tests/activedriverwgs_shards.nf.test, never run by the pipeline: the ACTIVEDRIVERWGS
// module three times on the committed fixture.
//   ADWGS_ONE_SHARD     shard size 5000: every element id in one ActiveDriverWGS call
//   ADWGS_SHARD_PER_ID  shard size 1: one ActiveDriverWGS call per element id
//   ADWGS_BOUNDARY      shard size 5000 on mutations_boundary.tsv: the fixture plus two SNVs from new
//                       patients around elem_A (BED chr5 1220000 1220600), one at 1-based 1,220,000
//                       (the base 5' of the BED interval) and one at 1,220,601 (the base after its end)
//
include { ACTIVEDRIVERWGS as ADWGS_ONE_SHARD    } from '../../modules/local/activedriverwgs/main'
include { ACTIVEDRIVERWGS as ADWGS_SHARD_PER_ID } from '../../modules/local/activedriverwgs/main'
include { ACTIVEDRIVERWGS as ADWGS_BOUNDARY     } from '../../modules/local/activedriverwgs/main'

workflow ADWGS_SHARDS {
    take:
    mutations          // [ meta, mutations.tsv ]
    mutations_boundary // [ meta, mutations_boundary.tsv ]
    elements           // [ meta, elements.bed ]

    main:
    ADWGS_ONE_SHARD(mutations, elements)
    ADWGS_SHARD_PER_ID(mutations, elements)
    ADWGS_BOUNDARY(mutations_boundary, elements)

    emit:
    one_shard    = ADWGS_ONE_SHARD.out.results
    shard_per_id = ADWGS_SHARD_PER_ID.out.results
    boundary     = ADWGS_BOUNDARY.out.results
}
