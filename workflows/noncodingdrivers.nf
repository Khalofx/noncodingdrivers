/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_noncodingdrivers_pipeline'
include { CSVTK_CONCAT as COHORT_CONCAT     } from '../modules/nf-core/csvtk/concat/main'
include { GAWK as COHORT_DEDUP              } from '../modules/nf-core/gawk/main'
include { UNTAR as UNTAR_ONCODRIVEFML_GENOME } from '../modules/nf-core/untar/main'
include { DRIVER_DISCOVERY                  } from '../subworkflows/local/driver_discovery'
include { RUNSUMMARYMQC                     } from '../modules/local/runsummarymqc'
include { GAWK as GENE_TSS                  } from '../modules/nf-core/gawk/main'
include { GAWK as SIGPROFILER_INPUT         } from '../modules/nf-core/gawk/main'
include { SIGPROFILER                       } from '../modules/nf-core/sigprofiler/main'
include { EXPLAIN                           } from '../subworkflows/local/explain'
include { peakFileMeta                      } from '../subworkflows/local/utils_nfcore_noncodingdrivers_pipeline'
include { resolveReferences                 } from '../subworkflows/local/utils_nfcore_noncodingdrivers_pipeline'
include { INPUT_HARMONISATION               } from '../subworkflows/local/input_harmonisation'
include { VARIANT_FLAGS                     } from '../subworkflows/local/variant_flags'
include { NULL_CALIBRATION                  } from '../subworkflows/local/null_calibration'
include { POSTFILTERFLAGS                   } from '../modules/local/postfilterflags'
include { WRITEDRIVEREVIDENCE               } from '../modules/local/writedriverevidence'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NONCODINGDRIVERS {

    take:
    ch_samplesheet // channel: [ val(meta), mutations, vcf, index, maf ] read in from --input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    //
    // Reference set of the run: datasheet assets/references/<genome>.yml with explicit
    // parameters taking precedence (PIPELINE_INITIALISATION has validated it)
    //
    def refs = resolveReferences()

    def ch_rows = ch_samplesheet.branch { meta, _mutations, _vcf, _index, _maf ->
        mutations: meta.entry == 'mutations'
        variants: true
    }

    //
    // SUBWORKFLOW: vcf, sarek and maf rows -> per-sample canonical tables and PASS VCFs (GRCh38).
    // Each row first passes the assembly check against the FASTA of its declared build.
    //
    def ch_harmonised  = channel.empty()
    def ch_cohort_vcfs = channel.empty()
    if (refs.fasta) {
        def ch_fasta = channel.value([ [ id: 'GRCh38' ], file(refs.fasta, checkIfExists: true), file(refs.fasta_fai, checkIfExists: true) ])
        def ch_dict  = channel.value(params.allow_hg19_liftover
            ? [ [ id: 'GRCh38' ], file(refs.fasta_dict, checkIfExists: true) ]
            : [ [:], [] ])
        def ch_source = channel.value(params.allow_hg19_liftover && refs.liftover_source_fasta
            ? [ [ id: 'hg19' ], file(refs.liftover_source_fasta, checkIfExists: true), file(refs.liftover_source_fai, checkIfExists: true) ]
            : [ [:], [], [] ])
        def ch_chain = channel.value(params.allow_hg19_liftover && refs.liftover_chain
            ? [ [ id: 'hg19ToHg38' ], file(refs.liftover_chain, checkIfExists: true) ]
            : [ [:], [] ])
        INPUT_HARMONISATION(
            ch_rows.variants.map { meta, _mutations, vcf, index, maf -> [ meta, vcf, index, maf ] },
            ch_fasta,
            ch_dict,
            ch_source,
            ch_chain,
            file(refs.rename_chrs, checkIfExists: true),
            params.assembly_min_match,
            params.assembly_min_records
        )
        ch_harmonised = INPUT_HARMONISATION.out.mutations
        // MutEnricher (QC-only) reads the per-sample PASS VCFs of each cohort
        ch_cohort_vcfs = INPUT_HARMONISATION.out.vcf
            .map { meta, vcf, tbi -> [ meta.cohort, vcf, tbi ] }
            .groupTuple()
            .map { cohort, vcfs, tbis -> [ [ id: cohort, cohort: cohort, assembly: 'GRCh38' ], vcfs, tbis ] }
    }

    //
    // One canonical mutation table per cohort. csvtk concat aligns the tables by column name;
    // a gawk program (conf/modules.config) then drops exact duplicate rows. Every row keeps its
    // sample column, so donors stay distinct after the merge. mutations rows enter as given;
    // vcf, sarek and maf rows enter as the per-sample tables of INPUT_HARMONISATION.
    //
    def ch_cohort_tables = ch_rows.mutations
        .map { meta, mutations, _vcf, _index, _maf -> [ meta.cohort, mutations ] }
        .mix(ch_harmonised.map { meta, table -> [ meta.cohort, table ] })
        .groupTuple()
        .map { cohort, tables -> [ [ id: cohort, cohort: cohort, assembly: 'GRCh38' ], tables.sort { table -> [ table.name, table.toString() ].join('\t') } ] }
    COHORT_CONCAT(ch_cohort_tables, 'tsv', 'tsv')
    COHORT_DEDUP(COHORT_CONCAT.out.csv.map { meta, tsv -> [ meta, tsv, 'mutations.tsv' ] }, [], false)
    def ch_mutations = COHORT_DEDUP.out.output

    //
    // References for the callers
    //
    def ch_elements = channel.value([ [ id: 'elements', assembly: 'GRCh38' ], file(refs.elements_bed, checkIfExists: true) ])


    //
    // Variant flags: mask and germline flags per mutation row, element-level fractions and the
    // callable fraction of each element. Added files only; nothing downstream reads them for
    // ranking.
    //
    def ch_variant_flags = channel.empty()
    def ch_element_flags = channel.empty()
    if (params.run_variant_flags) {
        def germline = !params.skip_germline_flags
        VARIANT_FLAGS(
            ch_mutations,
            ch_elements,
            channel.value([
                [ id: 'masks' ],
                file(refs.blacklist_bed, checkIfExists: true),
                file(refs.segdup_bed, checkIfExists: true),
                file(refs.low_mappability_bed, checkIfExists: true),
                file(refs.umap_bed, checkIfExists: true)
            ]),
            germline ? channel.value([ [ id: 'gnomad' ], file(refs.gnomad_echtvar, checkIfExists: true) ]) : channel.empty(),
            germline ? channel.value([ [ id: 'GRCh38' ], file(refs.fasta, checkIfExists: true), file(refs.fasta_fai, checkIfExists: true) ]) : channel.empty(),
            file("${projectDir}/assets/awk", checkIfExists: true),
            germline
        )
        ch_variant_flags = VARIANT_FLAGS.out.variant_flags
        ch_element_flags = VARIANT_FLAGS.out.element_flags
    }
    def ch_ofml_score = channel.empty()
    def ofml_config   = []
    def ofml_genome   = []
    if (params.run_oncodrivefml) {
        ch_ofml_score = channel.value([
            [ id: 'oncodrivefml_score' ],
            file(params.oncodrivefml_score, checkIfExists: true),
            file("${params.oncodrivefml_score}.tbi", checkIfExists: true)
        ])
        ofml_config = file(params.oncodrivefml_conf ?: "${projectDir}/assets/oncodrivefml.CADD_GRCh38.conf", checkIfExists: true)
        // Offline bgdata genome: a directory, or a .tar.gz of it unpacked by the nf-core untar module
        if (refs.oncodrivefml_genome && refs.oncodrivefml_genome.toString().endsWith('.tar.gz')) {
            UNTAR_ONCODRIVEFML_GENOME(channel.of([ [ id: 'oncodrivefml_genome' ], file(refs.oncodrivefml_genome, checkIfExists: true) ]))
            ofml_genome = UNTAR_ONCODRIVEFML_GENOME.out.untar.map { _meta, genome -> genome }.first()
        } else if (refs.oncodrivefml_genome) {
            ofml_genome = file(refs.oncodrivefml_genome, checkIfExists: true)
        }
    }

    // MutEnricher (QC-only) reads the per-sample PASS VCFs of the vcf, sarek and maf rows
    // (ch_cohort_vcfs above); mutations-tier cohorts have no VCF and are not run through MutEnricher.
    def ch_mappable = params.run_mutenricher
        ? channel.value([
            [ id: 'mappable' ],
            file(params.mutenricher_mappable_regions, checkIfExists: true),
            file("${params.mutenricher_mappable_regions}.tbi", checkIfExists: true)
        ])
        : channel.empty()

    def callers = [
        params.run_activedriverwgs ? 'activedriverwgs' : null,
        params.run_oncodrivefml    ? 'oncodrivefml'    : null,
        params.run_mutenricher     ? 'mutenricher'     : null
    ].findAll { caller -> caller != null }
    def combined_callers = callers.findAll { caller -> caller != 'mutenricher' }

    //
    // SUBWORKFLOW: callers -> harmonised calls -> anchored on element_id -> combined p-values
    //
    DRIVER_DISCOVERY(
        ch_mutations,
        ch_cohort_vcfs,
        ch_elements,
        ch_ofml_score,
        ofml_config,
        ofml_genome,
        ch_mappable,
        callers
    )

    //
    // MODULE: run summary for MultiQC (counts describing the run; no candidate ranking)
    //
    //
    // Null calibration (in-house, off by default): the same callers on SIMULATENULL copies of
    // each cohort; reported only
    //
    def ch_null_calibration = channel.empty()
    if (params.run_null_calibration) {
        NULL_CALIBRATION(
            ch_mutations,
            DRIVER_DISCOVERY.out.consensus,
            DRIVER_DISCOVERY.out.raw,
            ch_elements,
            channel.value([ [ id: 'GRCh38' ], file(refs.fasta, checkIfExists: true), file(refs.fasta_fai, checkIfExists: true) ]),
            ch_ofml_score,
            ofml_config,
            ofml_genome,
            callers,
            params.null_replicates,
            params.seed
        )
        ch_null_calibration = NULL_CALIBRATION.out.calibration
    }
    //
    // Post-filter flags and the long-form evidence table (D33). Added files only: the consensus
    // is not filtered or re-ranked; tiers are assigned downstream.
    //
    def ch_evidence = channel.empty()
    if (params.run_driver_evidence) {
        // coordinate provenance of each cohort's mutations: GRCh37/hg19 rows are lifted
        def ch_coordinates = ch_samplesheet
            .map { meta, _mutations, _vcf, _index, _maf -> [ meta.cohort, (meta.assembly_declared ?: 'GRCh38') in [ 'GRCh37', 'hg19' ] ? 'lifted' : 'native' ] }
            .unique()
            .groupTuple()
            .map { cohort, labels -> [ cohort, labels.unique().size() > 1 ? 'mixed' : labels[0] ] }
        def ch_cons_by_id = DRIVER_DISCOVERY.out.consensus.map { meta, consensus -> [ meta.id, meta, consensus ] }
        def ch_flags_by_id = ch_element_flags.map { meta, flags -> [ meta.id, flags ] }
        POSTFILTERFLAGS(
            ch_cons_by_id
                .join(ch_flags_by_id, remainder: true)
                .filter { _id, meta, _consensus, _flags -> meta != null }
                .map { _id, meta, consensus, flags -> [ meta, consensus, flags ?: [] ] }
        )
        WRITEDRIVEREVIDENCE(
            ch_cons_by_id
                .join(ch_flags_by_id, remainder: true)
                .filter { _id, meta, _consensus, _flags -> meta != null }
                .join(POSTFILTERFLAGS.out.postfilter.map { meta, postfilter -> [ meta.id, postfilter ] })
                .join(ch_null_calibration.map { meta, calibration -> [ meta.id, calibration ] }, remainder: true)
                .filter { _id, meta, _consensus, _flags, _postfilter, _calibration -> meta != null }
                .join(ch_coordinates)
                .map { _id, meta, consensus, flags, postfilter, calibration, coordinates ->
                    [ meta, consensus, flags ?: [], postfilter, calibration ?: [], coordinates ]
                },
            file("${projectDir}/assets/evidence_sources.tsv", checkIfExists: true),
            workflow.manifest.version
        )
        ch_evidence = WRITEDRIVEREVIDENCE.out.evidence
    }
    RUNSUMMARYMQC(
        DRIVER_DISCOVERY.out.consensus
            .join(DRIVER_DISCOVERY.out.raw.groupTuple(size: combined_callers.size()))
    )
    ch_multiqc_files = ch_multiqc_files.mix(RUNSUMMARYMQC.out.sections.map { _meta, sections -> sections })

    //
    // References for the annotation layers. The nearest-gene table uses a sorted gene TSS BED:
    // --gene_tss_bed, or one built from --gencode_gtf by assets/awk/gene_tss.awk (TSS = gene
    // start on +, gene end on -). Without either, no nearest-gene table is written.
    //
    def ch_gene_tss = channel.empty()
    if (params.gene_tss_bed) {
        ch_gene_tss = channel.value([ [ id: 'gene_tss' ], file(params.gene_tss_bed, checkIfExists: true) ])
    } else if (refs.gencode_gtf) {
        def gtf = file(refs.gencode_gtf, checkIfExists: true)
        GENE_TSS(
            channel.of([ [ id: gtf.name.replaceAll(/\.gtf(\.gz)?$/, '') ], gtf, 'tss.bed' ]),
            file("${projectDir}/assets/awk/gene_tss.awk", checkIfExists: true),
            false
        )
        ch_gene_tss = GENE_TSS.out.output
    } else {
        log.warn "Neither --gene_tss_bed nor --gencode_gtf is set: the nearest-gene table is not written"
    }

    // ENCODE peak files, one per tissue and layer: <regulatory_peaks_dir>/<tissue>/<layer>.<accession>.bed.gz
    def ch_peaks = channel.empty()
    if (params.run_regulatory_context) {
        ch_peaks = channel.fromPath("${refs.regulatory_peaks_dir}/*/*.{bed,bed.gz}", checkIfExists: true)
            .map { peak -> [ peakFileMeta(peak), peak ] }
    }

    def ch_vep_cache = channel.empty()
    if (params.run_vep) {
        ch_vep_cache = channel.value([ [ id: "vep_${params.vep_cache_version}" ], file(refs.vep_cache, checkIfExists: true) ])
    }
    def motifs = params.run_motifbreakr
        ? file(refs.jaspar_motifs, checkIfExists: true)
        : []

    //
    // SUBWORKFLOW: annotation layers (none of them filters or re-ranks a candidate)
    //
    EXPLAIN(
        DRIVER_DISCOVERY.out.consensus,
        ch_mutations,
        ch_elements,
        ch_gene_tss,
        ch_peaks,
        ch_vep_cache,
        file("${projectDir}/assets/awk", checkIfExists: true),
        motifs,
        [ regctx: params.run_regulatory_context, vep: params.run_vep, motifs: params.run_motifbreakr,
          vep_cache_version: params.vep_cache_version ]
    )
    ch_multiqc_files = ch_multiqc_files.mix(EXPLAIN.out.mqc)

    //
    // Optional cohort-level QC: SigProfiler signature exposures per sample (nf-core sigprofiler).
    // No element filter and no scoring use these exposures.
    //
    if (params.run_sigprofiler) {
        SIGPROFILER_INPUT(ch_mutations.map { meta, tsv -> [ meta, tsv, 'tsv' ] }, file("${projectDir}/assets/awk/sigprofiler_input.awk", checkIfExists: true), false)
        SIGPROFILER(
            SIGPROFILER_INPUT.out.output,
            'GRCh38',
            params.sigprofiler_genome ? file(params.sigprofiler_genome, checkIfExists: true) : [],
            []
        )
        // per-sample COSMIC SBS96 activities reach MultiQC as custom content (assets/multiqc_config.yml)
        ch_multiqc_files = ch_multiqc_files.mix(SIGPROFILER.out.results_sigprofiler.map { _meta, results -> results })
    }

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'noncodingdrivers_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'noncodingdrivers'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    emit:
    mutations      = ch_mutations                      // channel: [ val(meta), path(mutations) ] one table per cohort
    consensus      = DRIVER_DISCOVERY.out.consensus    // channel: [ val(meta), path(consensus) ]
    candidates     = EXPLAIN.out.candidates            // channel: [ val(meta), path(candidates) ]
    variant_flags  = ch_variant_flags                  // channel: [ val(meta), path(variant_flags.tsv) ]
    element_flags  = ch_element_flags                  // channel: [ val(meta), path(element_flags.tsv) ]
    null_calibration = ch_null_calibration             // channel: [ val(meta), path(null_calibration.tsv) ]
    evidence       = ch_evidence                       // channel: [ val(meta), path(driver_evidence.tsv) ]
    multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                       // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
