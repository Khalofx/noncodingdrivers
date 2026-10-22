process ACTIVEDRIVERWGS {
    tag "${meta.id}"
    label 'process_high'

    // No Conda package exists for ActiveDriverWGS (CRAN only). Checked 2026-10-20: no recipe
    // in bioconda-recipes (b95d6dd); absent from conda-forge noarch and linux-64 repodata.
    // So this module has no conda directive.
    // The image is built from containers/activedriverwgs/Dockerfile (R 4.5.3 + CRAN
    // ActiveDriverWGS 1.2.1 + BSgenome.Hsapiens.UCSC.hg38 1.4.5). Pin by digest in Step 7.
    container 'ghcr.io/khalofx/noncodingdrivers-activedriverwgs:1.2.1'

    input:
    tuple val(meta), path(mutations)     // canonical cohort table: sample chrom pos ref alt vtype (GRCh38, chr-prefixed)
    tuple val(meta2), path(elements)     // CDS-excluded element catalogue, headerless BED (>= 4 columns)

    output:
    tuple val(meta), path("${prefix}.activedriverwgs.tsv"),          emit: results
    tuple val(meta), path("${prefix}.activedriverwgs.untested.txt"), emit: untested, optional: true
    tuple val("${task.process}"), val('activedriverwgs'), eval('cat activedriverwgs.version.txt'), topic: versions, emit: versions_activedriverwgs
    tuple val("${task.process}"), val('r-base'), eval('cat r.version.txt'), topic: versions, emit: versions_r

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''              // --ref-genome hg38 --shard-size N (conf/modules.config)
    prefix   = task.ext.prefix ?: "${meta.id}"
    // Transitional escape hatch for a site R install (Karakoram): ext.rbin / ext.rlib point at the
    // site Rscript and its prepared library. Container runs use the image's Rscript and library.
    // Removed in Step 7 once the image is the only supported runtime.
    def rbin = task.ext.rbin ?: 'Rscript'
    def rlib = task.ext.rlib ?: ''
    def libs = rlib ? ".libPaths(c('${rlib}', .libPaths()));" : ''
    def rlib_arg = rlib ? "--rlib ${rlib}" : ''
    """
    ${rbin} -e "${libs} cat(as.character(packageVersion('ActiveDriverWGS')))" > activedriverwgs.version.txt
    ${rbin} -e "cat(paste(R.version\\\$major, R.version\\\$minor, sep='.'))" > r.version.txt

    # Canonical table -> ActiveDriverWGS input (chr, pos1, pos2, ref, alt, patient), pos1 == pos2.
    # format_muts accepts only single-base A/C/G/T or '-' alleles and one bad row kills a whole
    # shard, so only SNVs (single-base ACGT on both sides) on the 24 standard chromosomes are
    # passed. Indels and MNVs are still scored by OncodriveFML.
    awk 'BEGIN{FS=OFS="\\t"; print "chr","pos1","pos2","ref","alt","patient";
              for(i=1;i<=22;i++) ok["chr"i]=1; ok["chrX"]=1; ok["chrY"]=1}
         NR==1{for(i=1;i<=NF;i++)c[\$i]=i; next}
         {
           chr=\$c["chrom"]; ref=\$c["ref"]; alt=\$c["alt"];
           if(!(chr in ok)) next;
           if(ref ~ /^[ACGTacgt]\$/ && alt ~ /^[ACGTacgt]\$/)
             print chr,\$c["pos"],\$c["pos"],ref,alt,\$c["sample"];
         }' \\
        ${mutations} > adwgs_input.tsv

    # Element catalogue -> BED4 (chrom, start, end, id) on the 24 primary chromosomes; rows with
    # non-numeric or empty coordinates or ids are dropped, and an empty result fails the task.
    # Same output as the headerless-BED branch of the earlier adapter script (md5 checked on the
    # GRCh38 catalogue, tests/regression/PROVENANCE.md).
    awk 'BEGIN{FS=OFS="\\t"; for(i=1;i<=22;i++) p[i ""]=1; p["X"]=p["Y"]=1}
         {
           chrom=\$1; sub(/^chr/, "", chrom);
           if(!(chrom in p)) next;
           if(\$2 !~ /^[0-9]+\$/ || \$3 !~ /^[0-9]+\$/ || \$3+0 <= \$2+0 || \$4 == "") next;
           print "chr" chrom, \$2, \$3, \$4; n++
         }
         END{ if(n == 0){ print "ERROR: no valid primary-chromosome elements" > "/dev/stderr"; exit 2 } }' \\
        ${elements} > adwgs_elements.bed

    ${rbin} \$(command -v run_activedriverwgs.R) \\
        --mutations adwgs_input.tsv \\
        --elements adwgs_elements.bed \\
        --out ${prefix}.activedriverwgs.tsv \\
        --cores ${task.cpus} \\
        ${rlib_arg} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    // Stub runs need no R: the version files hold fixed strings here; real runs write them
    // from the R installation actually used (script block above).
    """
    touch ${prefix}.activedriverwgs.tsv
    touch ${prefix}.activedriverwgs.untested.txt
    echo "1.2.1" > activedriverwgs.version.txt
    echo "4.5.3" > r.version.txt
    """
}
