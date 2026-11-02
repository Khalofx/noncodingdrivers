process ONCODRIVECLUSTL {
    tag "${meta.id}"
    label 'process_medium'

    // Licence: OncodriveCLUSTL is AGPL-3.0-or-later (IRB Barcelona) with no non-commercial clause, so
    // it needs no licence gate and its image is public; the pipeline runs it as a separate program.
    // Reported alongside the consensus by default (--consensus_callers).
    // No conda directive: no oncodriveclustl/bgdata package in Bioconda or conda-forge (2026-10-29);
    // environment.yml documents the runtime. Image: containers/oncodriveclustl (explicit conda lock,
    // hashed pip lock, AGPL text, bundled COSMIC CGC list reduced to its header).
    // ghcr.io/khalofx/noncodingdrivers-oncodriveclustl:1.2.0, pinned by digest
    container 'ghcr.io/khalofx/noncodingdrivers-oncodriveclustl@sha256:8925475d0045b13dac79e4ac7ee1852db1b123843c9c31c9cfed2141f5aa5c7f'

    input:
    tuple val(meta), path(mutations), path(elements) // canonical cohort table (sample chrom pos ref alt vtype); element rows (headerless BED >= 4 columns), usually the ids with >= -emut SNVs
    path bgdata_genome                               // offline bgdata hg38 genome directory, or [] for bgdata's own cache

    output:
    tuple val(meta), path("${prefix}.oncodriveclustl.tsv"),          emit: results
    tuple val(meta), path("${prefix}.oncodriveclustl.clusters.tsv"), emit: clusters
    tuple val(meta), path("${prefix}.oncodriveclustl.log"),          emit: log
    tuple val("${task.process}"), val('oncodriveclustl'), eval('cat oncodriveclustl.version.txt'), topic: versions, emit: versions_oncodriveclustl

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''             // OncodriveCLUSTL options, with --seed (conf/modules.config)
    def shard_size = task.ext.shard_size ?: 5000     // element ids per OncodriveCLUSTL run
    prefix         = task.ext.prefix ?: "${meta.id}"
    if (!(args =~ /(^|\s)--seed\s+\d+/)) {
        error("ONCODRIVECLUSTL needs --seed <n> in ext.args: without it the simulations are not reproducible.")
    }
    if (args.contains('region_normalized')) {
        error("ONCODRIVECLUSTL: -sigcalc region_normalized is not supported: the signature is computed once from all SNVs " +
              "(frequencies) and shared by every shard, and the regions are prefiltered to the analysed ids.")
    }
    def genome_cache = bgdata_genome ? """
    # Offline genome: build the bgdata mapping inside the task, so it never holds a host-only path.
    export BGDATA_LOCAL="\$PWD/.bgdata"
    export BGDATA_BUILDS="\$PWD/bgdata-builds.conf"
    export BGDATA_OFFLINE=TRUE
    mkdir -p "\$BGDATA_LOCAL"
    printf '[paths]\\ndatasets/genomereference/hg38 = %s\\n' "\$PWD/${bgdata_genome}" > "\$BGDATA_BUILDS"
    """ : ''
    """
    # bgconfig and matplotlib write runtime configuration on first import; keep it task-local.
    export XDG_CONFIG_HOME="\$PWD/.config"
    export MPLCONFIGDIR="\$PWD/.matplotlib"
    mkdir -p "\$XDG_CONFIG_HOME" "\$MPLCONFIGDIR"
    ${genome_cache}
    python -c 'import importlib.metadata as m; print(m.version("oncodriveclustl"))' > oncodriveclustl.version.txt

    # Canonical table -> OncodriveCLUSTL input (CHROMOSOME, POSITION, REF, ALT, SAMPLE): SNVs on the
    # 24 primary chromosomes (the contigs of the regions; the signature is computed from these SNVs),
    # 'chr' stripped for the bgdata hg38 contigs. The file name is the signature key (clustl_input).
    awk 'BEGIN{FS=OFS="\\t"; for(i=1;i<=22;i++) p[i ""]=1; p["X"]=p["Y"]=1; print "CHROMOSOME","POSITION","REF","ALT","SAMPLE"}
         NR==1{for(i=1;i<=NF;i++)c[\$i]=i; next}
         (\$c["vtype"]=="SNP"){ch=\$c["chrom"]; sub(/^chr/,"",ch); if(ch in p) print ch,\$c["pos"],\$c["ref"],\$c["alt"],\$c["sample"]}' \\
        ${mutations} > clustl_input.tsv

    # Element rows -> OncodriveCLUSTL regions (CHROMOSOME, START, END, ELEMENT, SYMBOL) on the 24
    # primary chromosomes, one file per shard of ${shard_size} element ids; every fragment of an id
    # stays in its shard, so OncodriveCLUSTL joins them into one element. OncodriveCLUSTL reads
    # [START, END] as 1-based inclusive: START = BED start + 1. ELEMENT = SYMBOL = element id.
    awk -v n=${shard_size} 'BEGIN{FS=OFS="\\t"; for(i=1;i<=22;i++) p[i ""]=1; p["X"]=p["Y"]=1}
         {
           ch=\$1; sub(/^chr/, "", ch);
           if(!(ch in p)) next;
           if(\$2 !~ /^[0-9]+\$/ || \$3 !~ /^[0-9]+\$/ || \$3+0 <= \$2+0 || \$4 == "") next;
           if(!(\$4 in shard)){ shard[\$4] = int(k / n) + 1; k++ }
           f = sprintf("regions.%05d.tsv", shard[\$4]);
           if(!(f in seen)){ print "CHROMOSOME","START","END","ELEMENT","SYMBOL" > f; seen[f] = 1 }
           print ch, \$2 + 1, \$3, \$4, \$4 > f
         }' ${elements}

    # One OncodriveCLUSTL run per shard. The first run computes the mutational signature from all
    # SNVs of the cohort (OncodriveCLUSTL's frequencies method); every later shard reads that
    # signature (-sig), so all shards share it. With one shard the result is OncodriveCLUSTL's own.
    # A shard without any element above the mutation cutoff writes no table (OncodriveCLUSTL stops
    # with "No element found with enough mutations"); its elements are untested.
    for shard in regions.*.tsv; do
        [ -e "\$shard" ] || continue
        out="out.\${shard#regions.}"; out="\${out%.tsv}"
        sig=''
        if [ -s signature.json ]; then sig='-sig signature.json'; fi
        set +e
        oncodriveclustl -i clustl_input.tsv -r "\$shard" -o "\$out" -g hg38 -c ${task.cpus} \$sig ${args} > "\$out.stderr.log" 2>&1
        rc=\$?
        set -e
        if [ \$rc -ne 0 ]; then
            if grep -q 'No element found with enough mutations' "\$out.stderr.log"; then
                echo "\$shard: no element with enough mutations; no table" >&2
                continue
            fi
            cat "\$out.stderr.log" >&2
            exit \$rc
        fi
        if [ ! -s signature.json ]; then
            python -c 'import glob, json, pickle, sys; json.dump(pickle.load(open(glob.glob(sys.argv[1] + "/cache/clustl_input_kmer_*.pickle")[0], "rb")), open("signature.json", "w"))' "\$out"
        fi
    done

    # Shard tables -> one cohort table. Two or more: the q-values are recomputed over all shards
    # (Benjamini-Hochberg with OncodriveCLUSTL's own statsmodels call) and the cluster ranks follow.
    tables=\$(ls out.*/elements_results.txt 2> /dev/null || true)
    n_tables=\$(echo "\$tables" | grep -c . || true)
    if [ "\$n_tables" -eq 0 ]; then
        printf 'SYMBOL\\tENSID\\tCGC\\tCHROMOSOME\\tSTRAND\\tLENGTH\\tTOTAL_MUT\\tCLUSTERED_MUT\\tCLUSTERS\\tSIM_CLUSTERS\\tSCORE\\tP_EMPIRICAL\\tQ_EMPIRICAL\\tP_ANALYTICAL\\tQ_ANALYTICAL\\tP_TOPCLUSTER\\tQ_TOPCLUSTER\\n' > ${prefix}.oncodriveclustl.tsv
        printf 'RANK\\tSYMBOL\\tENSID\\tCGC\\tCHROMOSOME\\tSTRAND\\tCOORDINATES\\tMAX_COORD\\tWIDTH\\tN_MUT\\tN_SAMPLES\\tFRA_UNIQ_SAMPLES\\tSCORE\\tP\\n' > ${prefix}.oncodriveclustl.clusters.tsv
        echo "OncodriveCLUSTL: no element with enough mutations in any shard; header-only tables" | tee -a ${prefix}.oncodriveclustl.log >&2
    elif [ "\$n_tables" -eq 1 ]; then
        cp \$tables ${prefix}.oncodriveclustl.tsv
        cp \$(dirname \$tables)/clusters_results.tsv ${prefix}.oncodriveclustl.clusters.tsv
    else
        merge_oncodriveclustl_shards.py \\
            --elements \$tables \\
            --clusters \$(for t in \$tables; do echo \$(dirname \$t)/clusters_results.tsv; done) \\
            --out-elements ${prefix}.oncodriveclustl.tsv \\
            --out-clusters ${prefix}.oncodriveclustl.clusters.tsv
    fi
    for log in out.*/results.log; do
        [ -e "\$log" ] || continue
        echo "== \$(dirname \$log)" >> ${prefix}.oncodriveclustl.log
        cat "\$log" >> ${prefix}.oncodriveclustl.log
    done
    touch ${prefix}.oncodriveclustl.log
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'SYMBOL\\tENSID\\tCGC\\tCHROMOSOME\\tSTRAND\\tLENGTH\\tTOTAL_MUT\\tCLUSTERED_MUT\\tCLUSTERS\\tSIM_CLUSTERS\\tSCORE\\tP_EMPIRICAL\\tQ_EMPIRICAL\\tP_ANALYTICAL\\tQ_ANALYTICAL\\tP_TOPCLUSTER\\tQ_TOPCLUSTER\\n' > ${prefix}.oncodriveclustl.tsv
    printf 'RANK\\tSYMBOL\\tENSID\\tCGC\\tCHROMOSOME\\tSTRAND\\tCOORDINATES\\tMAX_COORD\\tWIDTH\\tN_MUT\\tN_SAMPLES\\tFRA_UNIQ_SAMPLES\\tSCORE\\tP\\n' > ${prefix}.oncodriveclustl.clusters.tsv
    touch ${prefix}.oncodriveclustl.log
    echo "1.2.0" > oncodriveclustl.version.txt
    """
}
