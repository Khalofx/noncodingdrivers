process RUNSUMMARYMQC {
    tag "${meta.id}"
    label 'process_single'

    // MultiQC custom content describing the run (counts only; no candidate ranking), written by
    // bin/run_summary_mqc.py from the consensus table and the callers' own output tables.
    // No conda directive: the script is pipeline-owned; the utils image carries its python.
    container 'ghcr.io/khalofx/noncodingdrivers-utils:1.0.0'

    input:
    tuple val(meta), path(consensus), path(caller_tables, stageAs: 'callers/*')

    output:
    tuple val(meta), path("*_mqc.yaml"), emit: sections
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Paths only: the caller name is derived from each file name inside the (unit-tested) script.
    // A single table arrives as a Path rather than a list, hence the flatten.
    // Sorted by name, so the command line does not depend on which caller finished first.
    def tables = [ caller_tables ].flatten().sort { table -> table.name }.collect { table -> "--caller-table ${table}" }.join(' ')
    """
    run_summary_mqc.py \\
        --consensus ${consensus} \\
        ${tables} \\
        --sample ${meta.id} \\
        --outdir . \\
        ${args}
    """

    stub:
    """
    printf 'id: "noncodingdrivers_overview"\\nsection_name: "stub"\\nplot_type: "table"\\ndata:\\n    "%s":\\n        elements_tested: 0\\n' "${meta.id}" \\
        > ${meta.id}_noncodingdrivers_overview_mqc.yaml
    """
}
