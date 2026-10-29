process MUTENRICHER {
    tag "${meta.id}"
    label 'process_medium'

    // MutEnricher 1.3.3 (MIT), QC-only. Its element p-values correlate with ActiveDriverWGS
    // (Spearman rho +0.67 in two cohorts) and not with OncodriveFML: it replicates the
    // mutation-burden family rather than adding independent evidence, so its native table is
    // published as is and never enters the p-value combination.
    // No conda directive: MutEnricher has no Bioconda or conda-forge package. The image
    // (containers/mutenricher) pins python 3.7, builds cyvcf2 from source and carries bgzip and
    // tabix from htslib 1.21.
    // ghcr.io/khalofx/noncodingdrivers-mutenricher:1.3.3, pinned by digest
    container 'ghcr.io/khalofx/noncodingdrivers-mutenricher@sha256:896eafb86136dff87a0fb6e2232d81ccdeb0737f0766f86c958d6716bae2deb0'

    input:
    tuple val(meta), path(vcfs, stageAs: 'vcfs/*'), path(tbis, stageAs: 'vcfs/*') // per-sample PASS VCFs (bgzip) and their .tbi
    tuple val(meta2), path(elements)                                              // element catalogue, headerless BED (>= 4 columns)
    tuple val(meta3), path(mappable), path(mappable_tbi)                          // BGZF + tabix callable/mappable regions BED

    output:
    tuple val(meta), path("${prefix}.mutenricher.tsv"), emit: results
    tuple val("${task.process}"), val('mutenricher'), eval('cat mutenricher.version.txt'), topic: versions, emit: versions_mutenricher

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    python /opt/MutEnricher/mutEnricher.py --version 2>&1 | grep -o -E '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1 > mutenricher.version.txt

    # Manifest: <absolute path> TAB <sample>, one per-sample VCF per line (sample = file name
    # without .vcf.gz). FILTER=PASS on every record is load-bearing: cyvcf2 reports PASS as None
    # and MutEnricher skips any other FILTER value, so a '.' FILTER yields zero variants silently.
    for vcf in vcfs/*.vcf.gz; do
        printf '%s\\t%s\\n' "\$PWD/\$vcf" "\$(basename "\$vcf" .vcf.gz)"
    done > vcf_files.txt

    # BED4 projection of the catalogue on the 24 primary chromosomes, the element universe of the
    # other callers (fragments keep their id). MutEnricher queries every region in every VCF, and a
    # contig the VCFs lack (e.g. an unplaced scaffold) aborts it: "could not create iterator".
    awk 'BEGIN{FS=OFS="\\t"; for(i=1;i<=22;i++) p["chr" i]=1; p["chrX"]=p["chrY"]=1}
         (\$1 in p) && NF >= 4 {print \$1, \$2, \$3, \$4}' ${elements} > me_elements.bed

    # 'noncoding' takes no --anno-type. -o names an OUTPUT DIRECTORY; the files inside carry --prefix.
    python /opt/MutEnricher/mutEnricher.py noncoding \\
        me_elements.bed vcf_files.txt \\
        -o ${prefix} \\
        --prefix ${prefix} \\
        --processors ${task.cpus} \\
        --mappable-regions ${mappable} \\
        ${args}

    # The native region table, published unmodified
    cp ${prefix}/${prefix}_region_WAP_hotspot_Fisher_enrichments.txt ${prefix}.mutenricher.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mutenricher.tsv
    echo "1.3.3" > mutenricher.version.txt
    """
}
