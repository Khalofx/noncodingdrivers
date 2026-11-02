process CANCERGENECONTEXT {
    tag "${meta.id}"
    label 'process_single'

    // Cancer gene-list labels of the nearest gene, the nearest protein-coding gene and the
    // ENCODE-rE2G targets of every significant element (bin/cancer_gene_context.py; Python
    // standard library only). Annotation only; nearest_gene.tsv itself is not changed.
    // ghcr.io/khalofx/noncodingdrivers-utils:1.0.1 (containers/utils), pinned by digest
    container 'ghcr.io/khalofx/noncodingdrivers-utils@sha256:9c8a5f128d8c64e581c977d58950ac700b6aabf9903e40af77e2e25db6599c6f'

    input:
    tuple val(meta), path(nearest), path(nearest_pc), path(links)
    path gene_list
    val list_name

    output:
    tuple val(meta), path("*.cancer_gene_context.tsv"), emit: context
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def pc     = nearest_pc ? "--nearest-pc ${nearest_pc}" : ''
    def re2g   = links ? "--links ${links}" : ''
    def name   = (list_name ?: gene_list.name).toString().replace("'", '')
    """
    cancer_gene_context.py \\
        --nearest ${nearest} \\
        ${pc} \\
        ${re2g} \\
        --gene-list ${gene_list} \\
        --list-name '${name}' \\
        --output ${prefix}.cancer_gene_context.tsv \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'element_id\\tnearest_in_list\\n' > ${prefix}.cancer_gene_context.tsv
    """
}
