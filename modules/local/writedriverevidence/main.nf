process WRITEDRIVEREVIDENCE {
    tag "${meta.id}"
    label 'process_single'

    // Long-form evidence table (D33; bin/write_driver_evidence.py), one row per cohort x element x
    // evidence item, in the layout of assets/schema_driver_evidence.json. No tiers: those are
    // assigned downstream (reports repository). `tables` holds the optional per-cohort tables of the
    // layers that are on (and the OncodriveCLUSTL null FPR table), recognised by file name; with
    // every layer off it is empty and the table is unchanged.
    // ghcr.io/khalofx/noncodingdrivers-utils:1.0.1 (containers/utils), pinned by digest
    container 'ghcr.io/khalofx/noncodingdrivers-utils@sha256:9c8a5f128d8c64e581c977d58950ac700b6aabf9903e40af77e2e25db6599c6f'

    input:
    tuple val(meta), path(consensus), path(element_flags), path(postfilter), path(calibration), val(coordinates), path(oncodriveclustl), path(tables, stageAs: 'tables/*')
    path(sources)
    val(pipeline_version)

    output:
    tuple val(meta), path("*.driver_evidence.tsv"), emit: evidence
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def optional = [
        element_flags ? "--element-flags ${element_flags}" : '',
        postfilter    ? "--postfilter ${postfilter}" : '',
        calibration   ? "--calibration ${calibration}" : '',
        oncodriveclustl ? "--oncodriveclustl ${oncodriveclustl}" : ''
    ].findAll { arg -> arg }.join(' ')
    // layer summary tables (S7A.33): file-name suffix -> option of bin/write_driver_evidence.py
    def table_options = [
        '.avi.element_summary.tsv'   : '--avi-summary',
        '.opencravat.tsv'            : '--opencravat',
        '.cosmic.element_summary.tsv': '--cosmic-elements',
        '.cancer_gene_context.tsv'   : '--cancer-genes',
        '.enhancer_gene_summary.tsv' : '--re2g-summary',
        '.ccre_v4_summary.tsv'       : '--ccre-summary',
        '.null_fpr.tsv'              : '--null-fpr'
    ]
    def layer_args = [ tables ].flatten().sort { t -> t.name }.collect { t ->
        def suffix = table_options.keySet().find { s -> t.name.endsWith(s) }
        suffix ? "${table_options[suffix]} ${t}" : ''
    }.findAll { arg -> arg }.join(' ')
    """
    write_driver_evidence.py \\
        --cohort ${meta.id} \\
        --consensus ${consensus} \\
        ${optional} \\
        ${layer_args} \\
        --sources ${sources} \\
        --coordinates ${coordinates} \\
        --pipeline-version ${pipeline_version} \\
        --output ${prefix}.driver_evidence.tsv \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'cohort\\telement_id\\tevidence_item\\tvalue\\tunit\\tstatus\\tthreshold\\tsource\\tsource_version\\tnative_or_lifted\\n' > ${prefix}.driver_evidence.tsv
    """
}
