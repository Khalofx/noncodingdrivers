process COSMICCONTEXT {
    tag "${meta.id}"
    label 'process_single'

    // Exact GRCh38 matching of the candidate mutations against the user's own COSMIC Non-Coding
    // Variants file (bin/cosmic_context.py; Python standard library only). Annotation only.
    // Count tables (no COSMIC identifiers) and the provenance JSON go to cosmic/; the COSMIC
    // records (exact_matches) and NOTICE.txt go to cosmic/restricted/ (conf/modules.config).
    // ghcr.io/khalofx/noncodingdrivers-utils:1.0.1 (containers/utils), pinned by digest
    container 'ghcr.io/khalofx/noncodingdrivers-utils@sha256:9c8a5f128d8c64e581c977d58950ac700b6aabf9903e40af77e2e25db6599c6f'

    input:
    tuple val(meta), path(candidates), path(consensus)
    path cosmic_ncv
    val release
    val expected_md5

    output:
    tuple val(meta), path("*.cosmic.mutation_summary.tsv"), emit: mutation_summary
    tuple val(meta), path("*.cosmic.element_summary.tsv") , emit: element_summary
    tuple val(meta), path("*.cosmic.provenance.json")     , emit: provenance
    tuple val(meta), path("*.cosmic.exact_matches.tsv")   , emit: exact_matches
    tuple val(meta), path("NOTICE.txt")                   , emit: notice
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def md5    = expected_md5 ? "--expected-md5 ${expected_md5}" : ''
    def cons   = consensus ? "--consensus ${consensus}" : ''
    """
    cosmic_context.py \\
        --cosmic ${cosmic_ncv} \\
        --release ${release} \\
        ${md5} \\
        --candidates ${candidates} \\
        ${cons} \\
        --cohort ${prefix} \\
        --outdir . \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'cohort\\tcosmic_release\\n' > ${prefix}.cosmic.mutation_summary.tsv
    printf 'cohort\\tcosmic_release\\n' > ${prefix}.cosmic.element_summary.tsv
    printf 'cohort\\tcosmic_release\\n' > ${prefix}.cosmic.exact_matches.tsv
    echo '{}' > ${prefix}.cosmic.provenance.json
    touch NOTICE.txt
    """
}
