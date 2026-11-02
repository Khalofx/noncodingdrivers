process OPENCRAVAT {
    tag "${meta.id}"
    label 'process_low'

    // OpenCRAVAT 3.2.4 in the pipeline image (containers/opencravat): its base modules and the
    // annotators ccre_screen, repeat and clinvar at the versions frozen in
    // containers/opencravat/modules.freeze.json. Annotation only: nothing downstream filters or
    // ranks on it.
    // GRCh38 only, no liftover anywhere: -l hg38 is fixed (any other -l in ext.args stops the task),
    // and the image holds neither the hg19 annotator nor a UCSC chain file (non-commercial licence),
    // so OpenCRAVAT cannot lift coordinates in either direction.
    // Telemetry: OpenCRAVAT 3.2.4 posts usage metrics after every run with no switch to stop it. The
    // image points every remote address at a closed local port; this task also routes all HTTP(S)
    // requests to that port, records every connection attempt (netaudit), fails if one is not local,
    // and deletes the copy of the metrics that OpenCRAVAT saves after the failed POST (it holds a
    // MAC-derived machine id).
    // ghcr.io/khalofx/noncodingdrivers-opencravat:3.2.4, pinned by digest
    container 'ghcr.io/khalofx/noncodingdrivers-opencravat@sha256:TODO_PIN_AFTER_PUBLICATION'

    input:
    tuple val(meta), path(variants)     // OpenCRAVAT native input: chrom, pos, strand, ref, alt, sample (GRCh38)

    output:
    tuple val(meta), path("${prefix}.opencravat.report.tsv"),    emit: report
    tuple val(meta), path("${prefix}.opencravat.sqlite"),        emit: sqlite
    tuple val(meta), path("${prefix}.opencravat.log"),           emit: log
    tuple val(meta), path("${prefix}.opencravat.network.tsv"),   emit: network
    tuple val(meta), path("${prefix}.opencravat.telemetry.tsv"), emit: telemetry
    tuple val("${task.process}"), val('opencravat'), eval("grep '^open-cravat' opencravat.versions.tsv | cut -f2"), topic: versions, emit: versions_opencravat
    tuple val("${task.process}"), val('opencravat-hg38'), eval("grep '^hg38' opencravat.versions.tsv | cut -f2"), topic: versions, emit: versions_hg38
    tuple val("${task.process}"), val('opencravat-ccre_screen'), eval("grep '^ccre_screen' opencravat.versions.tsv | cut -f2"), topic: versions, emit: versions_ccre_screen
    tuple val("${task.process}"), val('opencravat-repeat'), eval("grep '^repeat' opencravat.versions.tsv | cut -f2"), topic: versions, emit: versions_repeat
    tuple val("${task.process}"), val('opencravat-clinvar'), eval("grep '^clinvar' opencravat.versions.tsv | cut -f2"), topic: versions, emit: versions_clinvar

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''    // the annotators: -a ccre_screen repeat clinvar (conf/modules.config)
    prefix     = task.ext.prefix ?: "${meta.id}"
    def closed = 'http://127.0.0.1:9'
    // all HTTP(S) requests to the closed port; ext.closed_proxy = false (tests only) leaves the
    // network open, to show that the image configuration alone keeps every request local
    def proxy  = task.ext.closed_proxy == false ? '' : "http_proxy=${closed} https_proxy=${closed} HTTP_PROXY=${closed} HTTPS_PROXY=${closed} no_proxy= NO_PROXY="
    if (args =~ /(^|\s)(-l|--liftover)(\s|=|$)/) {
        error("OPENCRAVAT: the input is GRCh38 and -l hg38 is fixed; remove -l/--liftover from ext.args.")
    }
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("OPENCRAVAT runs in its image only (containers/opencravat: the OpenCRAVAT modules are installed when the image is built); conda is not supported.")
    }
    """
    # 1. Task-local OpenCRAVAT configuration: the image's, with the metrics folder in this task.
    #    Every remote address in it must be the closed local port.
    conf_image=\${OPENCRAVAT_IMAGE_CONF:?OPENCRAVAT needs the noncodingdrivers OpenCRAVAT image (containers/opencravat)}
    mkdir -p oc_metrics
    sed "s|^metrics_dir:.*|metrics_dir: \$PWD/oc_metrics|" "\$conf_image" > cravat-system.yml
    if grep -E '^[a-z_]+_url:' cravat-system.yml | grep -v -F ': ${closed}'; then
        echo "OPENCRAVAT_TELEMETRY: an address in the OpenCRAVAT configuration is not the closed local port ${closed}" >&2
        exit 1
    fi
    export OPENCRAVAT_SYSTEM_CONF_PATH=\$PWD/cravat-system.yml

    # 2. oc run on GRCh38 (-l hg38: no liftover of the input); every HTTP(S) request goes to the closed
    #    port; every connection attempt of OpenCRAVAT and its child processes is recorded.
    : > netaudit.raw
    env ${proxy} \\
        PYTHONPATH=\${OPENCRAVAT_NETAUDIT_DIR}\${PYTHONPATH:+:\$PYTHONPATH} OC_NETAUDIT_LOG=\$PWD/netaudit.raw \\
        oc run ${variants} \\
            -l hg38 \\
            -n ${prefix}.opencravat \\
            -t text \\
            -d . \\
            --mp ${task.cpus} \\
            ${args} \\
        > oc_run.out 2>&1 || { cat oc_run.out >&2; exit 1; }
    cat oc_run.out
    if ! grep -q '^Finished normally' oc_run.out; then
        echo "OPENCRAVAT: oc run did not finish normally (see ${prefix}.opencravat.log)" >&2
        exit 1
    fi
    mv ${prefix}.opencravat.tsv ${prefix}.opencravat.report.tsv

    # 3. Network audit: one row per call x target; any target that is not local fails the task
    awk -F '\\t' -v OFS='\\t' '{ n[\$2 OFS \$3 OFS \$4]++ }
        END { print "call", "host", "port", "n_calls", "is_local"
              for (k in n) { split(k, f, OFS); print k, n[k], (f[2] ~ /^(127\\.0\\.0\\.1|::1|localhost)\$/ ? "yes" : "no") } }' netaudit.raw \\
        | (read -r h; echo "\$h"; LC_ALL=C sort) > ${prefix}.opencravat.network.tsv
    non_local=\$(awk -F '\\t' 'NR > 1 && \$5 != "yes"' ${prefix}.opencravat.network.tsv | wc -l)
    if [ "\$non_local" -gt 0 ]; then
        echo "OPENCRAVAT_NETWORK: OpenCRAVAT tried to reach \$non_local non-local target(s); see ${prefix}.opencravat.network.tsv" >&2
        exit 1
    fi

    # 4. The metrics copy OpenCRAVAT saves after the failed POST leaves no output: count it, delete it
    n_metrics=\$(find oc_metrics -name 'ocmetric_*.json' | wc -l)
    rm -rf oc_metrics netaudit.raw
    {
        printf 'item\\tvalue\\n'
        grep -E '^(metrics|store|publish)_url:' cravat-system.yml | sed 's/: /\\t/'
        printf 'http_proxy\\t%s\\n' '${task.ext.closed_proxy == false ? 'none (test: open network)' : closed}'
        printf 'non_local_connection_targets\\t%s\\n' "\$non_local"
        printf 'metrics_payloads_saved_locally_and_deleted\\t%s\\n' "\$n_metrics"
    } > ${prefix}.opencravat.telemetry.tsv

    python - ${prefix}.opencravat.sqlite > opencravat.versions.tsv <<'PY'
    import sqlite3, sys
    info = dict(sqlite3.connect(sys.argv[1]).execute("select colkey, colval from info"))
    print("open-cravat\\t" + info["open-cravat"])
    print("hg38\\t" + info["_mapper"].split(":")[1])
    for entry in info["_annotators"].split(","):
        name, _, version = entry.partition(":")
        if version:
            print(f"{name}\\t{version}")
    PY
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '#CRAVAT Report\\n#Report level: variant\\n#\\nVariant Annotation\\t\\t\\t\\t\\tOriginal Input\\t\\t\\t\\nUID\\tChrom\\tPosition\\tRef Base\\tAlt Base\\tChrom\\tPos\\tReference allele\\tAlternate allele\\n' > ${prefix}.opencravat.report.tsv
    touch ${prefix}.opencravat.sqlite ${prefix}.opencravat.log
    printf 'call\\thost\\tport\\tn_calls\\tis_local\\n' > ${prefix}.opencravat.network.tsv
    printf 'item\\tvalue\\nmetrics_url\\thttp://127.0.0.1:9\\n' > ${prefix}.opencravat.telemetry.tsv
    printf 'open-cravat\\t3.2.4\\nhg38\\t1.47.2\\nccre_screen\\t1.0.2\\nclinvar\\t2026.10.01\\nrepeat\\t2020.10.16\\n' > opencravat.versions.tsv
    """
}
