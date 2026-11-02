#!/usr/bin/env bash
# Run inside the selected image with the repository mounted read-only at /pipeline:
#   docker run --rm --network none --entrypoint bash -v "$PWD:/pipeline:ro" <image> \
#       /pipeline/containers/smoke-test.sh <name>
set -euo pipefail

name="${1:?Specify utils, activedriverwgs, mutenricher, oncodrivefml, oncodriveclustl, combinepvalues or motifbreakr}"

# Nextflow's task wrapper runs ps inside every task container to collect CPU and memory metrics.
command -v ps > /dev/null || { echo "ps (procps) is missing from the ${name} image" >&2; exit 1; }

case "${name}" in
    utils)
        python -c 'import pandas, numpy, scipy.stats'
        ;;
    activedriverwgs)
        Rscript --vanilla -e "options(error = function() { traceback(2); quit(status = 1) }); library(ActiveDriverWGS); library(BSgenome.Hsapiens.UCSC.hg38); stopifnot(as.character(packageVersion('ActiveDriverWGS')) == '1.2.1'); stopifnot(as.character(packageVersion('BSgenome.Hsapiens.UCSC.hg38')) == '1.4.5'); stopifnot(length(BSgenome.Hsapiens.UCSC.hg38::Hsapiens[['chr1']]) == 248956422)"
        # the adapter must load and reject a call without arguments with its usage line
        if adapter_output=$(Rscript /pipeline/bin/run_activedriverwgs.R 2>&1); then
            echo 'Adapter unexpectedly accepted a call without arguments' >&2
            exit 1
        fi
        if ! grep -Fq 'usage: run_activedriverwgs.R' <<< "$adapter_output"; then
            printf '%s\n' "$adapter_output" >&2
            echo 'Adapter did not report the expected usage message' >&2
            exit 1
        fi
        ;;
    mutenricher)
        python /opt/MutEnricher/mutEnricher.py --version
        python /opt/MutEnricher/mutEnricher.py noncoding --help > /dev/null
        bgzip --version > /dev/null
        tabix --version > /dev/null
        ;;
    oncodrivefml)
        oncodrivefml --help > /dev/null
        python -c "import importlib.metadata as m; v = m.version('oncodrivefml'); assert v == '2.5.0', v"
        ;;
    oncodriveclustl)
        oncodriveclustl --help > /dev/null
        python -c "import importlib.metadata as m; v = m.version('oncodriveclustl'); assert v == '1.2.0', v"
        python -c "import oncodriveclustl.main, bgreference, bgsignature, bgdata, pycurl"
        # full AGPL-3.0 text and the modification notice ship with the image
        grep -q 'GNU AFFERO GENERAL PUBLIC LICENSE' /usr/share/licenses/oncodriveclustl/AGPL-3.0.txt
        grep -q 'CGC_all_Oct15_10_29_09_2018.tsv' /usr/share/licenses/oncodriveclustl/NOTICE
        # the COSMIC Cancer Gene Census list bundled in the wheel is reduced to its header line
        python -c "import os, oncodriveclustl; p = os.path.join(os.path.dirname(oncodriveclustl.__file__), 'data', 'CGC_all_Oct15_10_29_09_2018.tsv'); t = open(p).read(); assert t == 'Gene Symbol\\tTier\\n', repr(t)"
        ;;
    combinepvalues)
        Rscript -e "library(EmpiricalBrownsMethod); library(RobustRankAggreg); library(ACAT)"
        ;;
    motifbreakr)
        Rscript -e "library(motifbreakR); library(universalmotif); library(BSgenome.Hsapiens.UCSC.hg38); stopifnot(as.character(packageVersion('motifbreakR')) == '2.24.0', as.character(packageVersion('universalmotif')) == '1.28.0', as.character(packageVersion('BSgenome.Hsapiens.UCSC.hg38')) == '1.4.5'); stopifnot(as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38, 'chr5', 1295113, 1295113)) == 'G')"
        # The adapter end to end: TERT C228T (chr5:1295113 G>A) against JASPAR MA0076.3 (ELK4) must
        # give the ETS gain of the 10 Oct 2026 pilot (ref p 9.994507e-04, alt p 2.288818e-05).
        tmp=$(mktemp -d)
        printf 'sample\tchrom\tpos\tref\talt\tvtype\telement_id\nS1\tchr5\t1295113\tG\tA\tSNP\tpromoter_ENSR00000746745\n' > "$tmp/cand.tsv"
        awk '/^>/ { p = ($1 == ">MA0076.3") } p' /pipeline/assets/motifs/jaspar_core_vertebrates_2026.pfm > "$tmp/elk4.pfm"
        (cd "$tmp" && Rscript /pipeline/bin/run_motifbreakr.R --candidates cand.tsv --motifs elk4.pfm --prefix smoke --cores 1)
        if ! awk -F '\t' 'NR > 1 && $7 == "MA0076.3" && $12 == "gain" && sprintf("%.4g", $9) == "0.0009995" && sprintf("%.4g", $10) == "2.289e-05" { ok = 1 } END { exit ok ? 0 : 1 }' "$tmp/smoke.motifbreakr.per_tf.tsv"; then
            cat "$tmp/smoke.motifbreakr.per_tf.tsv" >&2
            echo 'TERT C228T ELK4 gain not reproduced' >&2
            exit 1
        fi
        ;;
    *)
        echo "Unknown image: ${name}" >&2
        exit 2
        ;;
esac
echo "smoke test passed: ${name}"
