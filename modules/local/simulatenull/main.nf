process SIMULATENULL {
    tag "${meta.id}.null${replicate}"
    label 'process_single'

    // In-house null simulation (D32; off by default, --run_null_calibration): each SNV moves to a
    // random position with the same trinucleotide within +/- 50 kb, in the same patient
    // (bin/simulate_null.py). No conda directive: the utils image carries python, numpy and pandas.
    // ghcr.io/khalofx/noncodingdrivers-utils:1.0.1 (containers/utils), pinned by digest
    container 'ghcr.io/khalofx/noncodingdrivers-utils@sha256:9c8a5f128d8c64e581c977d58950ac700b6aabf9903e40af77e2e25db6599c6f'

    input:
    tuple val(meta), path(mutations), val(replicate), val(seed)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), val(replicate), path("*.mutations.tsv"), emit: mutations
    tuple val(meta), val(replicate), path("*.null_summary.tsv"), emit: summary
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.null${replicate}"
    """
    simulate_null.py \\
        --mutations ${mutations} \\
        --fasta ${fasta} \\
        --fai ${fai} \\
        --seed ${seed} \\
        --output ${prefix}.mutations.tsv \\
        --summary ${prefix}.null_summary.tsv \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.null${replicate}"
    """
    printf 'sample\\tchrom\\tpos\\tref\\talt\\tvtype\\n' > ${prefix}.mutations.tsv
    printf 'metric\\tvalue\\nsnvs_out\\t0\\n' > ${prefix}.null_summary.tsv
    """
}
