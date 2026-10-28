//
// NULL_CALIBRATION: re-run the same callers on neutral (SIMULATENULL) copies of each cohort and
// report per-caller calibration against the observed run. Reported, never used to change
// p-values. SIMULATENULL and CALIBRATECALLERS are in-house code (D32), off by default
// (--run_null_calibration false). MutEnricher (QC-only) is not run on the nulls.
//
include { SIMULATENULL                        } from '../../../modules/local/simulatenull'
include { DRIVER_DISCOVERY as DRIVERS_ON_NULL } from '../driver_discovery'
include { CALIBRATECALLERS                    } from '../../../modules/local/calibratecallers'

workflow NULL_CALIBRATION {

    take:
    ch_mutations  // channel: [ val(meta), path(mutations) ]   observed cohort tables
    ch_consensus  // channel: [ val(meta), path(consensus) ]   observed consensus tables
    ch_raw        // channel: [ val(meta), path(caller table) ] observed caller tables
    ch_elements   // channel: [ val(meta2), path(elements) ]    element catalogue (value channel)
    ch_fasta      // channel: [ val(meta3), path(fasta), path(fai) ] (value channel)
    ch_ofml_score // channel: [ val(meta4), path(score), path(tbi) ] (value channel), or empty
    ofml_config   // path, or []
    ofml_genome   // path or value channel, or []
    callers       // list: callers of the observed run
    replicates    // integer: null replicates per cohort
    seed          // integer: replicate r uses seed + r

    main:
    def null_callers = callers - [ 'mutenricher' ]
    def with_adwgs = 'activedriverwgs' in null_callers
    SIMULATENULL(
        ch_mutations.combine(channel.of(1..replicates)).map { meta, tsv, r -> [ meta, tsv, r, seed + r ] },
        ch_fasta
    )
    // Each replicate is its own cohort id, so the reused subworkflow keeps them apart
    def ch_null_tables = SIMULATENULL.out.mutations
        .map { meta, r, tsv -> [ meta + [ id: "${meta.id}.null${r}".toString(), cohort: "${meta.id}.null${r}".toString(), null_of: meta.id, replicate: r ], tsv ] }
    DRIVERS_ON_NULL(ch_null_tables, channel.empty(), ch_elements, ch_ofml_score, ofml_config, ofml_genome, channel.empty(), null_callers)

    def ch_null_cons = DRIVERS_ON_NULL.out.consensus
        .map { meta, consensus -> [ meta.null_of, consensus ] }
        .groupTuple(size: replicates)
    def ch_null_adwgs = with_adwgs
        ? DRIVERS_ON_NULL.out.raw
            .filter { _meta, table -> table.name.endsWith('.activedriverwgs.tsv') }
            .map { meta, table -> [ meta.null_of, table ] }
            .groupTuple(size: replicates)
        : ch_null_cons.map { id, _consensus -> [ id, [] ] }
    def ch_obs_adwgs = with_adwgs
        ? ch_raw
            .filter { _meta, table -> table.name.endsWith('.activedriverwgs.tsv') }
            .map { meta, table -> [ meta.id, table ] }
        : ch_consensus.map { meta, _consensus -> [ meta.id, [] ] }
    CALIBRATECALLERS(
        ch_consensus
            .map { meta, consensus -> [ meta.id, meta, consensus ] }
            .join(ch_obs_adwgs)
            .join(ch_null_cons)
            .join(ch_null_adwgs)
            .map { _id, meta, observed, observed_adwgs, nulls, null_adwgs -> [ meta, observed, observed_adwgs, nulls, null_adwgs ] }
    )

    emit:
    calibration = CALIBRATECALLERS.out.calibration // channel: [ val(meta), path(null_calibration.tsv) ]
    qq          = CALIBRATECALLERS.out.qq          // channel: [ val(meta), path(null_qq.tsv) ]
    summary     = SIMULATENULL.out.summary         // channel: [ val(meta), val(replicate), path(null_summary.tsv) ]
}
