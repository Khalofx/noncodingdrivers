process POSTFILTERFLAGS {
    tag "${meta.id}"
    label 'process_single'

    // PCAWG-style post-filter flags per element (D33; bin/postfilter_flags.py): minimum support,
    // mappability, maximum sample fraction. Flag only by default; --postfilter_mode pcawg adds
    // q_postfilter and never changes consensus_q.
    // ghcr.io/khalofx/noncodingdrivers-utils:1.0.1 (containers/utils), pinned by digest
    container 'ghcr.io/khalofx/noncodingdrivers-utils@sha256:9c8a5f128d8c64e581c977d58950ac700b6aabf9903e40af77e2e25db6599c6f'

    input:
    tuple val(meta), path(consensus), path(element_flags)

    output:
    tuple val(meta), path("*.postfilter.tsv"), emit: postfilter
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def flags = element_flags ? "--element-flags ${element_flags}" : ''
    """
    postfilter_flags.py \\
        --consensus ${consensus} \\
        ${flags} \\
        --output ${prefix}.postfilter.tsv \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'element_id\\tpostfilter_status\\n' > ${prefix}.postfilter.tsv
    """
}
