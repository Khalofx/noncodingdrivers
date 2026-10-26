//
// EXPLAIN: annotation layers on top of the consensus. None of them filters or re-ranks an
// element or a mutation. Every analytic step is a published tool (bedtools, Ensembl VEP,
// motifbreakR); the in-house code is gawk format adapters (assets/awk/*.awk).
//
include { GAWK as SIGNIFICANT_BED              } from '../../../modules/nf-core/gawk/main'
include { GAWK as MUTATION_BED                 } from '../../../modules/nf-core/gawk/main'
include { GAWK as CANDIDATE_TABLE              } from '../../../modules/nf-core/gawk/main'
include { GAWK as MIDPOINT_BED                 } from '../../../modules/nf-core/gawk/main'
include { GAWK as NEAREST_TABLE                } from '../../../modules/nf-core/gawk/main'
include { GAWK as REGCTX_TABLE                 } from '../../../modules/nf-core/gawk/main'
include { GAWK as VEP_INPUT                    } from '../../../modules/nf-core/gawk/main'
include { BEDTOOLS_SORT as SORT_SIGNIFICANT    } from '../../../modules/nf-core/bedtools/sort/main'
include { BEDTOOLS_SORT as SORT_MIDPOINTS      } from '../../../modules/nf-core/bedtools/sort/main'
include { BEDTOOLS_SORT as SORT_PEAKS          } from '../../../modules/nf-core/bedtools/sort/main'
include { BEDTOOLS_INTERSECT as CANDIDATES     } from '../../../modules/nf-core/bedtools/intersect/main'
include { BEDTOOLS_CLOSEST as NEAREST_TSS      } from '../../../modules/nf-core/bedtools/closest/main'
include { BEDTOOLS_COVERAGE as REGCTX_OVERLAP  } from '../../../modules/nf-core/bedtools/coverage/main'
include { BEDTOOLS_MAP as REGCTX_SIGNAL        } from '../../../modules/nf-core/bedtools/map/main'
include { ENSEMBLVEP_VEP                       } from '../../../modules/nf-core/ensemblvep/vep/main'
include { MOTIFBREAKR                          } from '../../../modules/local/motifbreakr'

workflow EXPLAIN {

    take:
    ch_consensus   // channel: [ val(meta), path(consensus) ]
    ch_mutations   // channel: [ val(meta), path(mutations) ]   cohorts without a mutation table are skipped
    ch_elements    // channel: [ val(meta2), path(elements_bed) ] element catalogue (value channel)
    ch_gene_tss    // channel: [ val(meta3), path(tss_bed) ]     sorted gene TSS BED (value channel), or empty
    ch_peaks       // channel: [ val(meta4), path(peaks) ]       one item per ENCODE peak file, or empty
    ch_vep_cache   // channel: [ val(meta5), path(cache) ]       offline VEP cache (value channel), or empty
    awk_dir        // path: directory with the gawk adapter programs (assets/awk)
    motifs         // path: JASPAR PFM file, or []
    options        // map: [ regctx: bool, vep: bool, motifs: bool, vep_cache_version: int ]

    main:
    def ch_regctx      = channel.empty()
    def ch_vep         = channel.empty()
    def ch_motif_tf    = channel.empty()
    def ch_mqc         = channel.empty()

    // Calls-only cohorts have a consensus but no mutation table: log the skip, publish nothing.
    ch_consensus
        .join(ch_mutations, remainder: true)
        .filter { _meta, consensus, mutations -> consensus != null && mutations == null }
        .subscribe { meta, _consensus, _mutations ->
            log.info "[noncodingdrivers] ${meta.id}: no mutation table (calls-only input); candidate mutations and their annotation layers are skipped, no candidate file is written"
        }

    //
    // 1. Fragments of the significant elements (BED4). A cohort with no significant element
    //    is logged and skipped: no candidate, nearest-gene or annotation file is written.
    //
    SIGNIFICANT_BED(
        ch_consensus.combine(ch_elements).map { meta, consensus, _meta2, elements -> [ meta, [ consensus, elements ], 'bed' ] },
        awk_dir.resolve('significant_bed.awk'),
        false
    )
    SORT_SIGNIFICANT(SIGNIFICANT_BED.out.output, [], 'sorted.bed')
    def ch_significant = SORT_SIGNIFICANT.out.sorted.branch { _meta, bed ->
        none: bed.size() == 0
        some: true
    }
    ch_significant.none.subscribe { meta, _bed ->
        log.info "[noncodingdrivers] ${meta.id}: no element has is_consensus == True; candidate mutations, nearest genes and annotation layers are skipped, no candidate file is written"
    }

    //
    // 2. Patient-resolved candidate mutations: one row per sample x mutation x significant
    //    element (bedtools intersect -wa -wb on the element fragments).
    //
    MUTATION_BED(ch_mutations.map { meta, tsv -> [ meta, tsv, 'bed' ] }, awk_dir.resolve('mutation_bed.awk'), false)
    CANDIDATES(MUTATION_BED.out.output.join(ch_significant.some), [ [], [] ], 'bed')
    CANDIDATE_TABLE(CANDIDATES.out.intersect.map { meta, bed -> [ meta, bed, 'tsv' ] }, awk_dir.resolve('candidates_table.awk'), false)
    def ch_candidates = CANDIDATE_TABLE.out.output

    //
    // 3. Nearest gene (proximity only): bedtools closest from the element midpoint to the
    //    sorted GENCODE TSS BED; signed distance in the direction of transcription.
    //
    MIDPOINT_BED(
        ch_consensus.join(ch_significant.some).map { meta, consensus, _bed -> [ meta, consensus, 'bed' ] },
        awk_dir.resolve('midpoint_bed.awk'),
        false
    )
    SORT_MIDPOINTS(MIDPOINT_BED.out.output, [], 'sorted.bed')
    NEAREST_TSS(
        SORT_MIDPOINTS.out.sorted.combine(ch_gene_tss).map { meta, midpoints, _meta3, tss -> [ meta, midpoints, tss ] },
        []
    )
    NEAREST_TABLE(NEAREST_TSS.out.output.map { meta, bed -> [ meta, bed, 'tsv' ] }, awk_dir.resolve('nearest_table.awk'), false)

    //
    // 4. Regulatory context (opt-in): per staged ENCODE file, peak count and covered bases
    //    (bedtools coverage) and maximum signal (bedtools map), collected into one long table.
    //
    if (options.regctx) {
        SORT_PEAKS(ch_peaks, [], 'sorted.bed')
        def ch_sig_x_peaks = ch_significant.some
            .combine(SORT_PEAKS.out.sorted)
            .map { meta, bed, meta4, peaks ->
                [ meta + [ tissue: meta4.tissue, layer: meta4.layer, accession: meta4.accession, signal_col: meta4.signal_col ], bed, peaks ]
            }
        REGCTX_OVERLAP(ch_sig_x_peaks, [])
        REGCTX_SIGNAL(ch_sig_x_peaks, [ [], [] ])
        def layer_keys = [ 'tissue', 'layer', 'accession', 'signal_col' ]
        REGCTX_TABLE(
            REGCTX_OVERLAP.out.bed
                .mix(REGCTX_SIGNAL.out.mapped)
                .map { meta, bed -> [ meta.findAll { key, _value -> !(key in layer_keys) }, bed ] }
                .groupTuple()
                .map { meta, beds -> [ meta, beds.sort { bed -> bed.name }, 'tsv' ] },
            awk_dir.resolve('regctx_table.awk'),
            false
        )
        ch_regctx = REGCTX_TABLE.out.output
    }

    //
    // 5. Ensembl VEP on candidate mutations only (opt-in): default-format input of the distinct
    //    variants, tab output with the v0.1 field list. --regulatory is never used: the element
    //    catalogue (Regulatory Build 110) is the only regulatory build.
    //
    if (options.vep) {
        VEP_INPUT(ch_candidates.map { meta, tsv -> [ meta, tsv, 'txt' ] }, awk_dir.resolve('vep_input.awk'), false)
        ENSEMBLVEP_VEP(
            VEP_INPUT.out.output.map { meta, txt -> [ meta, txt, [] ] },
            'GRCh38',
            'homo_sapiens',
            options.vep_cache_version,
            ch_vep_cache,
            [ [], [] ],
            [],
            [ [], [] ]
        )
        ch_vep = ENSEMBLVEP_VEP.out.tab
        ch_mqc = ch_mqc.mix(ENSEMBLVEP_VEP.out.report.map { _meta, _process, _tool, report -> report })
    }

    //
    // 6. motifbreakR on candidate SNVs (opt-in; annotation and hypothesis only)
    //
    if (options.motifs && motifs) {
        MOTIFBREAKR(ch_candidates, motifs)
        ch_motif_tf = MOTIFBREAKR.out.per_tf
    }

    emit:
    candidates   = ch_candidates              // channel: [ val(meta), path(candidates.tsv) ]
    nearest      = NEAREST_TABLE.out.output   // channel: [ val(meta), path(nearest_gene.tsv) ]
    regctx       = ch_regctx                  // channel: [ val(meta), path(regulatory_context.tsv) ]
    vep          = ch_vep                     // channel: [ val(meta), path(vep.tab.gz) ]
    motif_per_tf = ch_motif_tf                // channel: [ val(meta), path(motifbreakr.per_tf.tsv) ]
    mqc          = ch_mqc                     // channel: path(report) for MultiQC
}
