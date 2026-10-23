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
    // One canonical mutation table per cohort. csvtk concat aligns the tables by column name;
    // a gawk program (conf/modules.config) then drops exact duplicate rows. Every row keeps its
    // sample column, so donors stay distinct after the merge. Only the mutations entry tier
    // reaches this point (PIPELINE_INITIALISATION rejects vcf, sarek and maf rows for now).
    //
    def ch_cohort_tables = ch_samplesheet
        .map { meta, mutations, _vcf, _index, _maf -> [ meta.cohort, mutations ] }
        .groupTuple()
        .map { cohort, tables -> [ [ id: cohort, cohort: cohort, assembly: 'GRCh38' ], tables.sort { table -> table.toString() } ] }
    COHORT_CONCAT(ch_cohort_tables, 'tsv', 'tsv')
    COHORT_DEDUP(COHORT_CONCAT.out.csv.map { meta, tsv -> [ meta, tsv, 'mutations.tsv' ] }, [], false)
    def ch_mutations = COHORT_DEDUP.out.output

    //
    // References for the callers
    //
    def ch_elements = channel.value([ [ id: 'elements', assembly: 'GRCh38' ], file(params.elements_bed, checkIfExists: true) ])

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
        if (params.oncodrivefml_genome && params.oncodrivefml_genome.toString().endsWith('.tar.gz')) {
            UNTAR_ONCODRIVEFML_GENOME(channel.of([ [ id: 'oncodrivefml_genome' ], file(params.oncodrivefml_genome, checkIfExists: true) ]))
            ofml_genome = UNTAR_ONCODRIVEFML_GENOME.out.untar.map { _meta, genome -> genome }.first()
        } else if (params.oncodrivefml_genome) {
            ofml_genome = file(params.oncodrivefml_genome, checkIfExists: true)
        }
    }

    // MutEnricher (QC-only) reads per-sample PASS VCFs, which only the vcf, sarek and maf entry
    // tiers provide; mutations-tier cohorts have no VCF and are not run through MutEnricher.
    // Until those tiers are added, the per-sample VCF channel is empty and MutEnricher has no task.
    def ch_cohort_vcfs = channel.empty()
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
    RUNSUMMARYMQC(
        DRIVER_DISCOVERY.out.consensus
            .join(DRIVER_DISCOVERY.out.raw.groupTuple(size: combined_callers.size()))
    )
    ch_multiqc_files = ch_multiqc_files.mix(RUNSUMMARYMQC.out.sections.map { _meta, sections -> sections })

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
    multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                       // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
