process COHORT_FILTER {
    tag "${meta.id}"
    label 'process_single'

    // One sample per donor (with donor information) and the hypermutator filter, applied once to the cohort table, so every caller, the null calibration
    // and the flags read the same samples (bin/filter_hypermutators.py; --hypermutator_filter).
    // ghcr.io/khalofx/noncodingdrivers-utils:1.0.1 (containers/utils), pinned by digest
    container 'ghcr.io/khalofx/noncodingdrivers-utils@sha256:9c8a5f128d8c64e581c977d58950ac700b6aabf9903e40af77e2e25db6599c6f'

    input:
    // staged under input/: the upstream table is also named <cohort>.mutations.tsv
    tuple val(meta), path(mutations, stageAs: 'input/*'), path(sheet_donors, stageAs: 'donors/samplesheet/*')
    path(donor_map, stageAs: 'donors/map/*')

    output:
    tuple val(meta), path("${prefix}.mutations.tsv")        , emit: mutations
    tuple val(meta), path("${prefix}.sample_filter.tsv")    , emit: sample_filter
    tuple val(meta), path("${prefix}.removed_mutations.tsv"), emit: removed
    tuple val(meta), path("${prefix}.sample_filter_mqc.tsv"), emit: mqc
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    prefix    = task.ext.prefix ?: "${meta.id}"
    def maps  = [ sheet_donors, donor_map ].flatten().findAll { f -> f }.collect { f -> "--donor-map ${f}" }.join(' ')
    """
    filter_hypermutators.py \\
        --input ${mutations} \\
        --prefix ${prefix} \\
        --cohort ${meta.cohort ?: meta.id} \\
        ${maps} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    head -n 1 ${mutations} > ${prefix}.removed_mutations.tsv
    cp ${mutations} ${prefix}.mutations.tsv
    printf 'sample\\tdonor\\tn_mutations\\tkept\\tremoved_by\\treason\\trule\\tmin_mutations\\tiqr_factor\\tq1\\tq3\\tiqr_threshold\\n' > ${prefix}.sample_filter.tsv
    printf 'cohort\\tsamples\\tsamples_removed\\n' > ${prefix}.sample_filter_mqc.tsv
    """
}
