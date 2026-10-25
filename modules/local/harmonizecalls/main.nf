process HARMONIZECALLS {
    tag "${meta.id}:${tool}"
    label 'process_single'

    // Format adapter only: selects the element id, p and q columns of one caller's native table
    // and renames them to the shared schema (element_id, p_<caller>, q_<caller>). No statistics
    // are computed. Same csvtk package and image as the installed nf-core csvtk modules.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/91/917edb71b915f07fa2838c20e3c731181d3d315cbf8a9bfead41412d2b4ae062/data':
    'community.wave.seqera.io/library/csvtk:0.37.0--113625988dd3285d' }"

    input:
    tuple val(meta), val(tool), path(raw)            // one native result table per (cohort, caller)

    output:
    tuple val(meta), path("${prefix}.${tool}.harmonized.tsv"), emit: harmonized
    tuple val("${task.process}"), val('csvtk'), eval("csvtk version | sed -e 's/csvtk v//g'"), topic: versions, emit: versions_csvtk

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''                   // the caller's native columns: --fields <id>,<p>,<q> (conf/modules.config)
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    csvtk cut --tabs ${args} ${raw} \\
        | csvtk rename --tabs --fields 1,2,3 --names element_id,p_${tool},q_${tool} \\
        > ${prefix}.${tool}.harmonized.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'element_id\\tp_${tool}\\tq_${tool}\\n' > ${prefix}.${tool}.harmonized.tsv
    """
}
