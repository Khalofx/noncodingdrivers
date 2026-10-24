process COMBINEPVALUES {
    tag "${meta.id}"
    label 'process_single'

    // R wrapper (bin/combine_pvalues.R) that only calls package functions:
    // EmpiricalBrownsMethod, ACAT (authors' package), RobustRankAggreg::rhoScores, p.adjust.
    // ACAT 0.91 is on neither CRAN nor conda: the image installs it from
    // github.com/yaowuliu/ACAT @ e737d7d (GPL-3) until an r-acat conda-forge recipe exists.
    // With -profile conda, ACAT is absent, so run with --skip_acat true.
    conda "${moduleDir}/environment.yml"
    container 'ghcr.io/khalofx/noncodingdrivers-combinepvalues:1.0.0'   // built from containers/combinepvalues/Dockerfile

    input:
    tuple val(meta), path(matched)                // element_id + p_<caller> (+ q_<caller>) per element

    output:
    tuple val(meta), path("${prefix}.consensus.tsv"),                emit: consensus
    tuple val(meta), path("${prefix}.consensus.ebm_parameters.tsv"), emit: ebm_parameters, optional: true
    tuple val("${task.process}"), val('empiricalbrownsmethod'), eval("grep '^EmpiricalBrownsMethod' combinepvalues.versions.tsv | cut -f2"), topic: versions, emit: versions_empiricalbrownsmethod
    tuple val("${task.process}"), val('robustrankaggreg'), eval("grep '^RobustRankAggreg' combinepvalues.versions.tsv | cut -f2"), topic: versions, emit: versions_robustrankaggreg
    tuple val("${task.process}"), val('acat'), eval("grep '^ACAT' combinepvalues.versions.tsv | cut -f2"), topic: versions, emit: versions_acat
    tuple val("${task.process}"), val('r-base'), eval("grep '^r-base' combinepvalues.versions.tsv | cut -f2"), topic: versions, emit: versions_r

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''                // --method ebm|acat --fdr 0.1 [--skip-acat]
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    combine_pvalues.R \\
        --input ${matched} \\
        --out ${prefix}.consensus.tsv \\
        --versions combinepvalues.versions.tsv \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.consensus.tsv ${prefix}.consensus.ebm_parameters.tsv
    printf 'EmpiricalBrownsMethod\\t1.34.0\\nRobustRankAggreg\\t1.2.1\\nACAT\\t0.91\\nr-base\\t4.4.3\\n' > combinepvalues.versions.tsv
    """
}
