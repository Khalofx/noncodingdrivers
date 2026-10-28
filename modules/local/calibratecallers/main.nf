process CALIBRATECALLERS {
    tag "${meta.id}"
    label 'process_single'

    // In-house calibration report (D32; off by default, --run_null_calibration): per-method
    // inflation factor and q < FDR counts on the SIMULATENULL cohorts, and QQ points
    // (bin/calibrate_callers.py). Reported only; no p- or q-value is changed.
    // ghcr.io/khalofx/noncodingdrivers-utils:1.0.1 (containers/utils), pinned by digest
    container 'ghcr.io/khalofx/noncodingdrivers-utils@sha256:9c8a5f128d8c64e581c977d58950ac700b6aabf9903e40af77e2e25db6599c6f'

    input:
    tuple val(meta), path(observed_consensus), path(observed_adwgs), path(null_consensus, stageAs: 'null/*'), path(null_adwgs, stageAs: 'null/*')

    output:
    tuple val(meta), path("*.null_calibration.tsv"), emit: calibration
    tuple val(meta), path("*.null_qq.tsv"), emit: qq
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // A single file arrives as a Path rather than a list, hence the flatten; sorted by name so the
    // command line does not depend on which replicate finished first.
    def nulls = [ null_consensus ].flatten().sort { f -> f.name }.join(' ')
    def null_adwgs_arg = null_adwgs ? "--null-adwgs " + [ null_adwgs ].flatten().sort { f -> f.name }.join(' ') : ''
    def observed_adwgs_arg = observed_adwgs ? "--observed-adwgs ${observed_adwgs}" : ''
    """
    calibrate_callers.py \\
        --cohort ${meta.id} \\
        --observed-consensus ${observed_consensus} \\
        ${observed_adwgs_arg} \\
        --null-consensus ${nulls} \\
        ${null_adwgs_arg} \\
        --calibration ${prefix}.null_calibration.tsv \\
        --qq ${prefix}.null_qq.tsv \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'cohort\\trun\\treplicate\\tmethod\\tlambda\\tstatus\\n' > ${prefix}.null_calibration.tsv
    printf 'cohort\\trun\\tmethod\\texpected_mlog10p\\tobserved_mlog10p\\n' > ${prefix}.null_qq.tsv
    """
}
