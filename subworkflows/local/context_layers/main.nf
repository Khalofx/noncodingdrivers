//
// CONTEXT_LAYERS: known-cancer and regulatory-target context of the significant elements, on
// top of EXPLAIN. Every layer is off by default and none of them filters, re-ranks or weights
// an element or a mutation; each publishes its own tables.
//   COSMIC         exact matches of the candidate mutations in the user's COSMIC NCV file
//   cancer genes   list labels of the nearest gene, nearest protein-coding gene and rE2G targets
//   ENCODE-rE2G    enhancer-gene links of the element fragments for the requested tissues
//   SCREEN V4      registry cCRE classes (cell-type-agnostic and, optionally, per biosample)
// The analytic steps are bedtools; the in-house code is gawk format adapters (assets/awk) and
// two standard-library Python scripts (bin/cosmic_context.py, bin/cancer_gene_context.py).
//
include { COSMICCONTEXT                               } from '../../../modules/local/cosmiccontext'
include { CANCERGENECONTEXT                           } from '../../../modules/local/cancergenecontext'
include { GAWK as PC_TSS                              } from '../../../modules/nf-core/gawk/main'
include { GAWK as RE2G_LINKS                          } from '../../../modules/nf-core/gawk/main'
include { GAWK as RE2G_SUMMARY                        } from '../../../modules/nf-core/gawk/main'
include { GAWK as CCRE_TABLE                          } from '../../../modules/nf-core/gawk/main'
include { GAWK as CCRE_SUMMARY                        } from '../../../modules/nf-core/gawk/main'
include { GAWK as CCRE_CELLTYPE_TABLE                 } from '../../../modules/nf-core/gawk/main'
include { BEDTOOLS_INTERSECT as RE2G_OVERLAP          } from '../../../modules/nf-core/bedtools/intersect/main'
include { BEDTOOLS_INTERSECT as CCRE_OVERLAP          } from '../../../modules/nf-core/bedtools/intersect/main'
include { BEDTOOLS_INTERSECT as CCRE_CELLTYPE_OVERLAP } from '../../../modules/nf-core/bedtools/intersect/main'
include { BEDTOOLS_CLOSEST as NEAREST_PC_TSS          } from '../../../modules/nf-core/bedtools/closest/main'

workflow CONTEXT_LAYERS {

    take:
    ch_consensus     // channel: [ val(meta), path(consensus) ]
    ch_significant   // channel: [ val(meta), path(significant.sorted.bed) ]  fragments, cohorts with >= 1 significant element
    ch_candidates    // channel: [ val(meta), path(candidates.tsv) ]
    ch_nearest       // channel: [ val(meta), path(nearest_gene.tsv) ], or empty
    ch_midpoints     // channel: [ val(meta), path(midpoints.sorted.bed) ]
    ch_gene_tss      // channel: [ val(meta2), path(tss.bed) ] sorted gene TSS BED (value channel), or empty
    ch_cosmic_ncv    // channel: path(Cosmic_NonCodingVariants_v<NNN>_GRCh38.tsv.gz) (value channel), or empty
    ch_re2g          // channel: [ val(meta3), path(links.bed.gz) ] one per staged rE2G file (tissue, biosample, accession), or empty
    ch_ccre_v4       // channel: [ val(meta4), path(cCREs.bed.gz) ] SCREEN V4 agnostic cCREs (value channel), or empty
    ch_ccre_celltype // channel: [ val(meta5), path(cCREs.bed.gz) ] one per biosample cCRE file (tissue, biosample, accession), or empty
    gene_list        // path: cancer gene list, or []
    awk_dir          // path: directory with the gawk adapter programs (assets/awk)
    options          // map: [ cosmic: bool, cosmic_release: str, cosmic_md5: str, re2g: bool, re2g_tissues: 'a,b',
                     //        re2g_n_files: 'a:N,b:M', ccre: bool, ccre_celltype: bool, gene_list_name: str ]

    main:
    def ch_cosmic_mutations = channel.empty()
    def ch_cosmic_elements  = channel.empty()
    def ch_cosmic_prov      = channel.empty()
    def ch_cosmic_exact     = channel.empty()
    def ch_links            = channel.empty()
    def ch_re2g_summary     = channel.empty()
    def ch_ccre             = channel.empty()
    def ch_ccre_summary     = channel.empty()
    def ch_ccre_ct          = channel.empty()
    def ch_cancer_genes     = channel.empty()
    def file_keys           = [ 'tissue', 'biosample', 'accession' ]

    //
    // 1. COSMIC (user-supplied NCV file): candidate table + consensus (rank, span, q)
    //
    if (options.cosmic) {
        COSMICCONTEXT(
            ch_candidates.join(ch_consensus),
            ch_cosmic_ncv,
            options.cosmic_release,
            options.cosmic_md5 ?: ''
        )
        ch_cosmic_mutations = COSMICCONTEXT.out.mutation_summary
        ch_cosmic_elements  = COSMICCONTEXT.out.element_summary
        ch_cosmic_prov      = COSMICCONTEXT.out.provenance
        ch_cosmic_exact     = COSMICCONTEXT.out.exact_matches
    }

    //
    // 2. ENCODE-rE2G: element fragments x each staged link file (bedtools intersect -wo), then
    //    the long links table and the per element x tissue summary (self-promoter split, nearest gene)
    //
    if (options.re2g) {
        RE2G_OVERLAP(
            ch_significant
                .combine(ch_re2g)
                .map { meta, bed, meta3, links -> [ meta + meta3.subMap(file_keys), bed, links ] },
            [ [], [] ],
            'bed'
        )
        RE2G_LINKS(
            RE2G_OVERLAP.out.intersect
                .map { meta, bed -> [ meta.findAll { key, _value -> !(key in file_keys) }, bed ] }
                .groupTuple()
                .map { meta, beds -> [ meta, beds.sort { bed -> bed.name }, 'tsv' ] },
            awk_dir.resolve('re2g_links.awk'),
            false
        )
        ch_links = RE2G_LINKS.out.output
        RE2G_SUMMARY(
            ch_significant
                .join(ch_links)
                .join(ch_nearest, remainder: true)
                .filter { _meta, bed, _links, _nearest -> bed != null }
                .map { meta, bed, links, nearest ->
                    [ meta + [ re2g_tissues: options.re2g_tissues, re2g_n_files: options.re2g_n_files ], [ bed, links ] + (nearest ? [ nearest ] : []), 'tsv' ]
                },
            awk_dir.resolve('re2g_summary.awk'),
            false
        )
        ch_re2g_summary = RE2G_SUMMARY.out.output
            .map { meta, tsv -> [ meta.findAll { key, _value -> !(key in [ 're2g_tissues', 're2g_n_files' ]) }, tsv ] }
    }

    //
    // 3. SCREEN Registry V4: agnostic cCRE classes (bedtools intersect -wao) and, optionally, the
    //    biosample-specific classes (bedtools intersect -wo per staged cell-type file)
    //
    if (options.ccre) {
        CCRE_OVERLAP(
            ch_significant.combine(ch_ccre_v4).map { meta, bed, _meta4, ccre -> [ meta, bed, ccre ] },
            [ [], [] ],
            'bed'
        )
        CCRE_TABLE(CCRE_OVERLAP.out.intersect.map { meta, bed -> [ meta, bed, 'tsv' ] }, awk_dir.resolve('ccre_v4_table.awk'), false)
        CCRE_SUMMARY(CCRE_OVERLAP.out.intersect.map { meta, bed -> [ meta, bed, 'tsv' ] }, awk_dir.resolve('ccre_v4_table.awk'), false)
        ch_ccre         = CCRE_TABLE.out.output
        ch_ccre_summary = CCRE_SUMMARY.out.output
        if (options.ccre_celltype) {
            CCRE_CELLTYPE_OVERLAP(
                ch_significant
                    .combine(ch_ccre_celltype)
                    .map { meta, bed, meta5, ccre -> [ meta + meta5.subMap(file_keys), bed, ccre ] },
                [ [], [] ],
                'bed'
            )
            CCRE_CELLTYPE_TABLE(
                CCRE_CELLTYPE_OVERLAP.out.intersect
                    .map { meta, bed -> [ meta.findAll { key, _value -> !(key in file_keys) }, bed ] }
                    .groupTuple()
                    .map { meta, beds -> [ meta, beds.sort { bed -> bed.name }, 'tsv' ] },
                awk_dir.resolve('ccre_celltype_table.awk'),
                false
            )
            ch_ccre_ct = CCRE_CELLTYPE_TABLE.out.output
        }
    }

    //
    // 4. Cancer gene list: nearest protein-coding gene (bedtools closest -d -t first on the
    //    protein_coding TSSs), then the list labels of nearest, protein-coding and rE2G genes
    //
    if (gene_list) {
        PC_TSS(ch_gene_tss.map { meta2, tss -> [ meta2, tss, 'bed' ] }, awk_dir.resolve('protein_coding_tss.awk'), false)
        NEAREST_PC_TSS(
            ch_midpoints.combine(PC_TSS.out.output).map { meta, midpoints, _meta2, tss -> [ meta, midpoints, tss ] },
            []
        )
        def ch_gene_inputs = ch_nearest.join(NEAREST_PC_TSS.out.output)
        if (options.re2g) {
            ch_gene_inputs = ch_gene_inputs
                .join(ch_links, remainder: true)
                .filter { _meta, nearest, _pc, _links -> nearest != null }
                .map { meta, nearest, pc, links -> [ meta, nearest, pc, links ?: [] ] }
        } else {
            ch_gene_inputs = ch_gene_inputs.map { meta, nearest, pc -> [ meta, nearest, pc, [] ] }
        }
        CANCERGENECONTEXT(ch_gene_inputs, gene_list, options.gene_list_name ?: '')
        ch_cancer_genes = CANCERGENECONTEXT.out.context
    }

    emit:
    cosmic_mutations = ch_cosmic_mutations   // channel: [ val(meta), path(cosmic.mutation_summary.tsv) ]
    cosmic_elements  = ch_cosmic_elements    // channel: [ val(meta), path(cosmic.element_summary.tsv) ]
    cosmic_prov      = ch_cosmic_prov        // channel: [ val(meta), path(cosmic.provenance.json) ]
    cosmic_exact     = ch_cosmic_exact       // channel: [ val(meta), path(cosmic.exact_matches.tsv) ]  restricted
    re2g_links       = ch_links              // channel: [ val(meta), path(enhancer_gene_links.tsv) ]
    re2g_summary     = ch_re2g_summary       // channel: [ val(meta), path(enhancer_gene_summary.tsv) ]
    ccre             = ch_ccre               // channel: [ val(meta), path(ccre_v4.tsv) ]
    ccre_summary     = ch_ccre_summary       // channel: [ val(meta), path(ccre_v4_summary.tsv) ]
    ccre_celltype    = ch_ccre_ct            // channel: [ val(meta), path(ccre_v4_celltype.tsv) ]
    cancer_genes     = ch_cancer_genes       // channel: [ val(meta), path(cancer_gene_context.tsv) ]
}
