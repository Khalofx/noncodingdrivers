process MAF2VCF {
    tag "${meta.id}"
    label 'process_single'

    // maf2vcf.pl (vcf2maf 1.6.22, Bioconda) writes one tumour's MAF rows as VCF, taking indel
    // anchor bases from the FASTA of the declared build. It drops rows whose Reference_Allele
    // disagrees with that FASTA and lists them in <prefix>.skipped.tsv (published). The build
    // itself is checked before this step, on the MAF rows (ASSEMBLYCHECK), because these drops
    // would otherwise hide a wrong build. maf2vcf.pl does not report its version.
    //
    // After maf2vcf.pl the VCF keeps the tumour column only (column 10; maf2vcf.pl adds the
    // matched normal as column 11), and a missing FILTER ('.') becomes PASS: a MAF lists called
    // variants, and FILTER-aware readers (MutEnricher) need PASS. A FILTER column of the MAF is
    // kept as maf2vcf.pl writes it. maf2vcf.pl writes no ##contig lines, which bcftools needs:
    // they are added from the FASTA index. Records are not sorted (BCFTOOLS_SORT follows).
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vcf2maf:1.6.22--hdfd78af_2' :
        'biocontainers/vcf2maf:1.6.22--hdfd78af_2' }"

    input:
    tuple val(meta), path(maf)                // one tumour (assets/awk/maf_split.awk), UCSC contig names
    tuple val(meta2), path(fasta), path(fai)  // FASTA of the declared build (GRCh38, or hg19 for liftover rows)

    output:
    tuple val(meta), path("${prefix}.vcf"),         emit: vcf
    tuple val(meta), path("${prefix}.skipped.tsv"), emit: skipped, optional: true
    tuple val("${task.process}"), val('vcf2maf'), eval("echo 1.6.22"), topic: versions, emit: versions_vcf2maf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    // WARN: maf2vcf.pl has no version option; update the version string with the container.
    """
    maf2vcf.pl \\
        --input-maf ${maf} \\
        --output-dir maf2vcf \\
        --output-vcf maf2vcf/${prefix}.maf2vcf.vcf \\
        --ref-fasta ${fasta} \\
        ${args}

    awk 'BEGIN { FS = OFS = "\\t" }
        NR == FNR { contig[++n] = \$1; length_of[\$1] = \$2; next }
        /^##FORMAT=<ID=AD,/ { has_ad = 1 }
        /^##FORMAT=<ID=DP,/ { has_dp = 1 }
        /^##contig=/ { has_contig = 1 }
        /^##/ { print; next }
        /^#CHROM/ {
            if (!has_ad) print "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\\"Allelic depths of REF and ALT\\">"
            if (!has_dp) print "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\\"Read depth\\">"
            if (!has_contig) for (i = 1; i <= n; i++) print "##contig=<ID=" contig[i] ",length=" length_of[contig[i]] ">"
            print "##FILTER=<ID=PASS,Description=\\"All filters passed (MAF row without a FILTER value)\\">"
            print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10; next
        }
        { if (\$7 == "." || \$7 == "") \$7 = "PASS"; print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10 }' \\
        ${fai} maf2vcf/${prefix}.maf2vcf.vcf > ${prefix}.vcf

    for f in maf2vcf/*.skipped.tsv; do
        if [ -s "\$f" ]; then mv "\$f" ${prefix}.skipped.tsv; fi
    done
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '##fileformat=VCFv4.2\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t${meta.id}\\n' > ${prefix}.vcf
    """
}
