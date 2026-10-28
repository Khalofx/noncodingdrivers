process ONCODRIVEFML {
    tag "${meta.id}"
    label 'process_high'

    // Licence: OncodriveFML ships an AGPL LICENSE plus an academic/non-commercial notice. The
    // module is never enabled by default; DRIVER_DISCOVERY runs it only when the user sets
    // --run_oncodrivefml true and --accept_oncodrivefml_licence true.
    // No conda directive: no oncodrivefml/bgdata package in Bioconda or conda-forge (2026-10-20).
    // Image: the baseline OncodriveFML image rebuilt with procps (containers/oncodrivefml).
    // The rights-holder image docker.io/bbglab/oncodrivefml:2.5.0 has no `ps` (checked 2026-10-21).
    // Tag 2.5.0, pinned by digest. The image is not public (OncodriveFML may not be redistributed
    // without a licence): build containers/oncodrivefml/ and set the container in a config (docs/usage.md).
    container 'ghcr.io/khalofx/noncodingdrivers-oncodrivefml@sha256:5a113bf33516c42a07c3c557732e157807f2a95508d80b77245140d43f64f2d2'
    // The image's entrypoint starts OncodriveFML; Docker must run Nextflow's task shell instead.
    containerOptions { workflow.containerEngine == 'docker' ? '--entrypoint ""' : '' }

    input:
    tuple val(meta), path(mutations)                 // canonical cohort table: sample chrom pos ref alt vtype
    tuple val(meta2), path(elements)                 // CDS-excluded element catalogue, headerless BED (>= 4 columns)
    tuple val(meta3), path(score), path(score_tbi)   // CADD GRCh38 v1.7 SNV track (bgzip) + .tbi, user-supplied
    path config                                      // OncodriveFML config naming the score file by basename (seed, cores, sampling)
    path bgdata_genome                               // offline bgdata hg38 genome directory, or [] for bgdata's own cache

    output:
    tuple val(meta), path("${prefix}.oncodrivefml.tsv"),        emit: results
    tuple val(meta), path("${prefix}.oncodrivefml.status.txt"), emit: status
    tuple val("${task.process}"), val('oncodrivefml'), eval('cat oncodrivefml.version.txt'), topic: versions, emit: versions_oncodrivefml

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def genome_cache = bgdata_genome ? """
    # Offline genome: build the bgdata mapping inside the task, so it never holds a host-only path.
    export BGDATA_LOCAL="\$PWD/.bgdata"
    export BGDATA_BUILDS="\$PWD/bgdata-builds.conf"
    export BGDATA_OFFLINE=TRUE
    mkdir -p "\$BGDATA_LOCAL"
    printf '[paths]\\ndatasets/genomereference/hg38 = %s\\n' "\$PWD/${bgdata_genome}" > "\$BGDATA_BUILDS"
    """ : ''
    """
    # bgconfig writes runtime configuration on first import; keep it task-local and writable.
    export XDG_CONFIG_HOME="\$PWD/.config"
    export MPLCONFIGDIR="\$PWD/.matplotlib"
    mkdir -p "\$XDG_CONFIG_HOME" "\$MPLCONFIGDIR"
    ${genome_cache}
    python -c 'from oncodrivefml import __version__; print(__version__)' > oncodrivefml.version.txt

    # Canonical table -> OncodriveFML input (CHROMOSOME, POSITION, REF, ALT, SAMPLE), SNVs only.
    # CADD GRCh38 uses Ensembl contig names (5, not chr5) and OncodriveFML queries the score
    # with the chromosome string verbatim, so 'chr' is stripped here and from the elements.
    awk 'BEGIN{FS=OFS="\\t"; print "CHROMOSOME","POSITION","REF","ALT","SAMPLE"}
         NR==1{for(i=1;i<=NF;i++)c[\$i]=i; next}
         (\$c["vtype"]=="SNP"){ch=\$c["chrom"]; sub(/^chr/,"",ch); print ch,\$c["pos"],\$c["ref"],\$c["alt"],\$c["sample"]}' \\
        ${mutations} > ofml_input.tsv

    # Element catalogue -> BED4 on the 24 primary chromosomes (same adapter as ACTIVEDRIVERWGS),
    # then the headered table OncodriveFML reads, with 'chr' stripped to match CADD.
    awk 'BEGIN{FS=OFS="\\t"; for(i=1;i<=22;i++) p[i ""]=1; p["X"]=p["Y"]=1}
         {
           chrom=\$1; sub(/^chr/, "", chrom);
           if(!(chrom in p)) next;
           if(\$2 !~ /^[0-9]+\$/ || \$3 !~ /^[0-9]+\$/ || \$3+0 <= \$2+0 || \$4 == "") next;
           print "chr" chrom, \$2, \$3, \$4; n++
         }
         END{ if(n == 0){ print "ERROR: no valid primary-chromosome elements" > "/dev/stderr"; exit 2 } }' \\
        ${elements} > ofml_elements.bed
    awk 'BEGIN{FS=OFS="\\t"; print "CHROMOSOME","START","END","ELEMENT","SYMBOL"}
         {ch=\$1; sub(/^chr/,"",ch); print ch,\$2,\$3,\$4,\$4}' \\
        ofml_elements.bed > ofml_elements.tsv

    # OncodriveFML can exit non-zero on a clean completion (for example when no element passes
    # correction), so success is judged by the result table; the exit code is kept as evidence.
    set +e
    oncodrivefml \\
        -i ofml_input.tsv \\
        -e ofml_elements.tsv \\
        --sequencing wgs \\
        -c ${config} \\
        -o oncodrivefml_out \\
        ${args}
    echo "oncodrivefml_exit_code=\$?" > ${prefix}.oncodrivefml.status.txt
    set -e

    if ls oncodrivefml_out/*.tsv.gz > /dev/null 2>&1; then
        gunzip -c oncodrivefml_out/*.tsv.gz > ${prefix}.oncodrivefml.tsv
    elif ls oncodrivefml_out/*.tsv > /dev/null 2>&1; then
        cp oncodrivefml_out/*.tsv ${prefix}.oncodrivefml.tsv
    else
        echo "ERROR: OncodriveFML produced no result table" >&2
        exit 1
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.oncodrivefml.tsv
    echo "oncodrivefml_exit_code=0" > ${prefix}.oncodrivefml.status.txt
    echo "2.5.0" > oncodrivefml.version.txt
    """
}
