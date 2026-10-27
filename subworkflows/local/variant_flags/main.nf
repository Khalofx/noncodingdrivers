//
// VARIANT_FLAGS: per-row mask and germline flags of each cohort mutation table, and their
// element-level aggregation. Flag, never filter: the mutation table, the callers and the
// consensus are unchanged. Every analytic step is a published tool (bedtools, echtvar,
// bcftools); the in-house code is gawk format adapters (assets/awk/*.awk).
//
include { GAWK as VARIANT_SITES                         } from '../../../modules/nf-core/gawk/main'
include { GAWK as INDEL_ANCHORS                         } from '../../../modules/nf-core/gawk/main'
include { GAWK as SITES_VCF                             } from '../../../modules/nf-core/gawk/main'
include { GAWK as VARIANT_FLAGS_TABLE                   } from '../../../modules/nf-core/gawk/main'
include { GAWK as FLAGS_BED                             } from '../../../modules/nf-core/gawk/main'
include { GAWK as ELEMENT_FLAGS_TABLE                   } from '../../../modules/nf-core/gawk/main'
include { BEDTOOLS_INTERSECT as MASK_BLACKLIST          } from '../../../modules/nf-core/bedtools/intersect/main'
include { BEDTOOLS_INTERSECT as MASK_SEGDUP             } from '../../../modules/nf-core/bedtools/intersect/main'
include { BEDTOOLS_INTERSECT as MASK_LOW_MAPPABILITY    } from '../../../modules/nf-core/bedtools/intersect/main'
include { BEDTOOLS_COVERAGE  as MASK_UMAP               } from '../../../modules/nf-core/bedtools/coverage/main'
include { BEDTOOLS_GETFASTA  as ANCHOR_BASES            } from '../../../modules/nf-core/bedtools/getfasta/main'
include { BCFTOOLS_SORT      as SORT_SITES_VCF          } from '../../../modules/nf-core/bcftools/sort/main'
include { ECHTVAR_ANNO       as GNOMAD_ANNO             } from '../../../modules/nf-core/echtvar/anno/main'
include { BCFTOOLS_QUERY     as GERMLINE_QUERY          } from '../../../modules/nf-core/bcftools/query/main'
include { BEDTOOLS_INTERSECT as ELEMENT_HITS            } from '../../../modules/nf-core/bedtools/intersect/main'
include { BEDTOOLS_SUBTRACT  as CALLABLE_TERRITORY      } from '../../../modules/nf-core/bedtools/subtract/main'
include { BEDTOOLS_COVERAGE  as ELEMENT_CALLABLE        } from '../../../modules/nf-core/bedtools/coverage/main'

workflow VARIANT_FLAGS {

    take:
    ch_mutations // channel: [ val(meta), path(mutations) ]  one canonical table per cohort
    ch_elements  // channel: [ val(meta2), path(elements) ]  element catalogue (value channel)
    ch_masks     // channel: [ val(meta3), path(blacklist), path(segdup), path(low_mappability), path(umap) ] (value channel)
    ch_gnomad    // channel: [ val(meta4), path(echtvar_zip) ] (value channel); used only when germline is true
    ch_fasta     // channel: [ val(meta5), path(fasta), path(fai) ] (value channel); used only when germline is true
    awk_dir      // path: directory with the gawk adapter programs (assets/awk)
    germline     // boolean: annotate gnomAD germline frequencies (echtvar anno)

    main:
    //
    // 1. Distinct variants with their span (BED4, key chrom:pos:ref:alt), and the mask flags:
    //    overlap counts with the ENCODE Blacklist v2 and the GIAB segdup and low-mappability
    //    stratifications; covered fraction of the span by Umap k100 unique sequence.
    //
    VARIANT_SITES(ch_mutations.map { meta, tsv -> [ meta, tsv, 'bed' ] }, awk_dir.resolve('variant_sites.awk'), false)
    def ch_sites_masks = VARIANT_SITES.out.output.combine(ch_masks)
    MASK_BLACKLIST(
        ch_sites_masks.map { meta, sites, _meta3, blacklist, _segdup, _lowmap, _umap -> [ meta, sites, blacklist ] },
        [ [], [] ],
        'bed'
    )
    MASK_SEGDUP(
        ch_sites_masks.map { meta, sites, _meta3, _blacklist, segdup, _lowmap, _umap -> [ meta, sites, segdup ] },
        [ [], [] ],
        'bed'
    )
    MASK_LOW_MAPPABILITY(
        ch_sites_masks.map { meta, sites, _meta3, _blacklist, _segdup, lowmap, _umap -> [ meta, sites, lowmap ] },
        [ [], [] ],
        'bed'
    )
    MASK_UMAP(
        ch_sites_masks.map { meta, sites, _meta3, _blacklist, _segdup, _lowmap, umap -> [ meta, sites, umap ] },
        []
    )
    def ch_flag_inputs = MASK_BLACKLIST.out.intersect
        .join(MASK_SEGDUP.out.intersect)
        .join(MASK_LOW_MAPPABILITY.out.intersect)
        .join(MASK_UMAP.out.bed)
        .map { meta, blacklist, segdup, lowmap, umap -> [ meta, [ blacklist, segdup, lowmap, umap ] ] }

    //
    // 2. Germline frequencies (optional): indel anchor bases from the FASTA, a sites-only VCF,
    //    echtvar anno against gnomAD v4.1 genomes, and the annotation as a key table.
    //
    if (germline) {
        INDEL_ANCHORS(ch_mutations.map { meta, tsv -> [ meta, tsv, 'bed' ] }, awk_dir.resolve('variant_sites.awk'), false)
        ANCHOR_BASES(INDEL_ANCHORS.out.output, ch_fasta.map { _meta5, fasta, _fai -> fasta })
        SITES_VCF(
            ANCHOR_BASES.out.fasta
                .join(VARIANT_SITES.out.output)
                .combine(ch_fasta)
                .map { meta, anchors, sites, _meta5, _fasta, fai -> [ meta, [ fai, anchors, sites ], 'vcf' ] },
            awk_dir.resolve('sites_vcf.awk'),
            false
        )
        SORT_SITES_VCF(SITES_VCF.out.output)
        GNOMAD_ANNO(SORT_SITES_VCF.out.vcf, ch_gnomad, 'vcf.gz')
        GERMLINE_QUERY(GNOMAD_ANNO.out.vcf.map { meta, vcf -> [ meta, vcf, [] ] }, [], [], [])
        ch_flag_inputs = ch_flag_inputs
            .join(GERMLINE_QUERY.out.output)
            .map { meta, masks, query -> [ meta, masks + [ query ] ] }
    }

    //
    // 3. Per-row flags table (one row per mutation-table row, same order)
    //
    VARIANT_FLAGS_TABLE(
        ch_flag_inputs
            .join(ch_mutations)
            .map { meta, inputs, tsv -> [ meta, inputs + [ tsv ], 'tsv' ] },
        awk_dir.resolve('variant_flags.awk'),
        false
    )

    //
    // 4. Element level: rows assigned to elements by the point rule, aggregated per element id,
    //    with the callable fraction of each element (Umap k100 unique minus Blacklist v2).
    //
    FLAGS_BED(VARIANT_FLAGS_TABLE.out.output.map { meta, tsv -> [ meta, tsv, 'bed' ] }, awk_dir.resolve('flags_bed.awk'), false)
    ELEMENT_HITS(
        FLAGS_BED.out.output.combine(ch_elements).map { meta, bed, _meta2, elements -> [ meta, bed, elements ] },
        [ [], [] ],
        'bed'
    )
    CALLABLE_TERRITORY(ch_masks.map { meta3, blacklist, _segdup, _lowmap, umap -> [ meta3, umap, blacklist ] })
    ELEMENT_CALLABLE(
        ch_elements.combine(CALLABLE_TERRITORY.out.bed).map { meta2, elements, _meta3, callable -> [ meta2, elements, callable ] },
        []
    )
    ELEMENT_FLAGS_TABLE(
        ELEMENT_HITS.out.intersect
            .combine(ELEMENT_CALLABLE.out.bed)
            .map { meta, hits, _meta2, callable -> [ meta, [ callable, hits ], 'tsv' ] },
        awk_dir.resolve('element_flags.awk'),
        false
    )

    emit:
    variant_flags = VARIANT_FLAGS_TABLE.out.output // channel: [ val(meta), path(variant_flags.tsv) ]
    element_flags = ELEMENT_FLAGS_TABLE.out.output // channel: [ val(meta), path(element_flags.tsv) ]
}
