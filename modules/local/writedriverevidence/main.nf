process WRITEDRIVEREVIDENCE {
    tag "${meta.id}"
    label 'process_single'

    // Long-form evidence table (D33; bin/write_driver_evidence.py), one row per cohort x element x
    // evidence item, in the layout of assets/schema_driver_evidence.json. No tiers: those are
    // assigned downstream (reports repository).
    container 'ghcr.io/khalofx/noncodingdrivers-utils:1.0.0'

    input:
    tuple val(meta), path(consensus), path(element_flags), path(postfilter), path(calibration), val(coordinates)
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
        calibration   ? "--calibration ${calibration}" : ''
    ].findAll { arg -> arg }.join(' ')
    """
    write_driver_evidence.py \\
        --cohort ${meta.id} \\
        --consensus ${consensus} \\
        ${optional} \\
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
