process MOTIFBREAKR {
    tag "${meta.id}"
    label 'process_high'

    // R adapter (bin/run_motifbreakr.R) that only calls package functions: universalmotif reads
    // the JASPAR PFMs, motifbreakR scores both alleles and calculatePvalue() gives allele p-values.
    // Annotation only: nothing downstream filters or ranks on the output.
    conda "${moduleDir}/environment.yml"
    // ghcr.io/khalofx/noncodingdrivers-motifbreakr:2.24.0, built from containers/motifbreakr (explicit conda lock); pinned by digest
    container 'ghcr.io/khalofx/noncodingdrivers-motifbreakr@sha256:8f6d759729e82d7d759ac0d6d2a19c00c00242efd29a643742348d2585125116'

    input:
    tuple val(meta), path(candidates)   // candidate table: chrom, pos, ref, alt, vtype (+ sample, element_id)
    path motifs                         // JASPAR-format PFM file

    output:
    tuple val(meta), path("${prefix}.motifbreakr.tsv"),        emit: motifbreakr
    tuple val(meta), path("${prefix}.motifbreakr.per_tf.tsv"), emit: per_tf
    tuple val("${task.process}"), val('motifbreakr'), eval("grep '^motifbreakR' motifbreakr.versions.tsv | cut -f2"), topic: versions, emit: versions_motifbreakr
    tuple val("${task.process}"), val('universalmotif'), eval("grep '^universalmotif' motifbreakr.versions.tsv | cut -f2"), topic: versions, emit: versions_universalmotif
    tuple val("${task.process}"), val('bsgenome.hsapiens.ucsc.hg38'), eval("grep '^BSgenome' motifbreakr.versions.tsv | cut -f2"), topic: versions, emit: versions_bsgenome
    tuple val("${task.process}"), val('r-base'), eval("grep '^r-base' motifbreakr.versions.tsv | cut -f2"), topic: versions, emit: versions_r

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''      // --threshold 1e-4 --pseudocount 0.8 (pilot settings; the defaults)
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    run_motifbreakr.R \\
        --candidates ${candidates} \\
        --motifs ${motifs} \\
        --prefix ${prefix} \\
        --cores ${task.cpus} \\
        --versions motifbreakr.versions.tsv \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'seqnames\\tstart\\tend\\tSNP_id\\tgeneSymbol\\tproviderId\\tRefpvalue\\tAltpvalue\\teffect\\n' > ${prefix}.motifbreakr.tsv
    printf 'chrom\\tpos\\tref\\talt\\tsnp_id\\ttf\\tmotif_id\\tmotif_length\\tref_pvalue\\talt_pvalue\\tbest_pvalue\\tdirection\\tfold_change\\tp_below_resolution\\tn_motifs_for_tf\\tmotifbreakr_effect\\n' > ${prefix}.motifbreakr.per_tf.tsv
    printf 'motifbreakR\\t2.24.0\\nuniversalmotif\\t1.28.0\\nBSgenome.Hsapiens.UCSC.hg38\\t1.4.5\\nr-base\\t4.5.3\\n' > motifbreakr.versions.tsv
    """
}
