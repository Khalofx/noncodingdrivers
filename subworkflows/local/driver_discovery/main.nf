//
// DRIVER_DISCOVERY: cohort mutation table -> driver callers -> harmonised calls -> anchored on
// element_id -> combined p-values (published R packages; --combination_method sets consensus_q).
// MutEnricher (QC-only) reads per-sample VCFs; its native table is published and never combined.
//
include { ACTIVEDRIVERWGS           } from '../../../modules/local/activedriverwgs'
include { ONCODRIVEFML              } from '../../../modules/local/oncodrivefml'
include { MUTENRICHER               } from '../../../modules/local/mutenricher'
include { HARMONIZECALLS            } from '../../../modules/local/harmonizecalls'
include { GAWK as ELEMENT_COORDS    } from '../../../modules/nf-core/gawk/main'
include { CSVTK_JOIN as ANCHOR_JOIN } from '../../../modules/nf-core/csvtk/join/main'
include { COMBINEPVALUES            } from '../../../modules/local/combinepvalues'

workflow DRIVER_DISCOVERY {

    take:
    ch_mutations   // channel: [ val(meta), path(mutations) ]               one canonical table per cohort
    ch_cohort_vcfs // channel: [ val(meta), [ vcf.gz ], [ tbi ] ]           per-sample PASS VCFs per cohort; used only by MutEnricher
    ch_elements    // channel: [ val(meta2), path(elements) ]               value channel; element catalogue (headerless BED)
    ch_ofml_score  // channel: [ val(meta3), path(score), path(score_tbi) ] value channel; used only by OncodriveFML
    ofml_config    // path: OncodriveFML configuration, or [] when OncodriveFML is not run
    ofml_genome    // path or value channel: offline bgdata hg38 genome directory, or []
    ch_mappable    // channel: [ val(meta4), path(bed.gz), path(tbi) ]      value channel; used only by MutEnricher
    callers        // list: callers to run, e.g. [ 'activedriverwgs', 'oncodrivefml' ]

    main:
    // One coordinate row per element id: chrom, first start, last end and the number of catalogue
    // rows (fragments left by the coding-sequence exclusion). Layout only (gawk, conf/modules.config).
    ELEMENT_COORDS(ch_elements.map { meta2, bed -> [ meta2, bed, 'tsv' ] }, [], false)
    def ch_coords = ELEMENT_COORDS.out.output.first()

    def ch_raw = channel.empty()
    if ('activedriverwgs' in callers) {
        ACTIVEDRIVERWGS(ch_mutations, ch_elements)
        ch_raw = ch_raw.mix(ACTIVEDRIVERWGS.out.results.map { meta, tsv -> [ meta, 'activedriverwgs', tsv ] })
    }
    if ('oncodrivefml' in callers) {
        ONCODRIVEFML(ch_mutations, ch_elements, ch_ofml_score, ofml_config, ofml_genome)
        ch_raw = ch_raw.mix(ONCODRIVEFML.out.results.map { meta, tsv -> [ meta, 'oncodrivefml', tsv ] })
    }
    // QC-only: published unmodified, not harmonised, not anchored, not combined
    def ch_mutenricher = channel.empty()
    if ('mutenricher' in callers) {
        MUTENRICHER(ch_cohort_vcfs, ch_elements, ch_mappable)
        ch_mutenricher = MUTENRICHER.out.results
    }
    def combined_callers = callers.findAll { caller -> caller != 'mutenricher' }

    // Format adapter per caller: native columns -> element_id, p_<caller>, q_<caller>
    HARMONIZECALLS(ch_raw)

    // Outer join of every caller's table with the element coordinates on element_id. Tables are
    // ordered by file name so the joined column order does not depend on which caller finished first.
    def ch_to_join = HARMONIZECALLS.out.harmonized
        .groupTuple(size: combined_callers.size())
        .combine(ch_coords)
        .map { meta, calls, _meta2, coords -> [ meta, calls.sort { f -> f.name } + [ coords ] ] }
    ANCHOR_JOIN(ch_to_join)

    // EBM, ACAT, RRA, votes and BH; --combination_method selects consensus_q
    COMBINEPVALUES(ANCHOR_JOIN.out.out_file)

    emit:
    coords         = ch_coords                          // channel: [ val(meta2), path(coords) ]
    raw            = ch_raw.map { meta, _caller, tsv -> [ meta, tsv ] } // channel: [ val(meta), path(caller table) ]
    mutenricher    = ch_mutenricher                     // channel: [ val(meta), path(native table) ]
    harmonized     = HARMONIZECALLS.out.harmonized      // channel: [ val(meta), path(harmonized) ]
    matched        = ANCHOR_JOIN.out.out_file           // channel: [ val(meta), path(matched) ]
    consensus      = COMBINEPVALUES.out.consensus       // channel: [ val(meta), path(consensus) ]
    ebm_parameters = COMBINEPVALUES.out.ebm_parameters  // channel: [ val(meta), path(tsv) ]
}
