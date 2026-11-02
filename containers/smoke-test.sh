#!/usr/bin/env bash
# Run inside the selected image with the repository mounted read-only at /pipeline:
#   docker run --rm --network none --entrypoint bash -v "$PWD:/pipeline:ro" <image> \
#       /pipeline/containers/smoke-test.sh <name>
set -euo pipefail

name="${1:?Specify utils, activedriverwgs, mutenricher, oncodrivefml, oncodriveclustl, combinepvalues, motifbreakr or opencravat}"

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
    opencravat)
        python -c "import importlib.metadata as m; v = m.version('open-cravat'); assert v == '3.2.4', v"
        # telemetry neutralised in the image: every remote address is the closed local port, the
        # metrics folder is empty
        python - "$OPENCRAVAT_IMAGE_CONF" <<'PY'
import sys, yaml
conf = yaml.safe_load(open(sys.argv[1]))
urls = {k: v for k, v in conf.items() if k.endswith("_url")}
assert urls and all(v == "http://127.0.0.1:9" for v in urls.values()), urls
assert conf["user_email_opt_out"] is True
PY
        test -z "$(ls -A /opt/opencravat/metrics)"
        # oc run on three synthetic SNVs (one, chr1:1000000 T>A, is in no ClinVar record), audited:
        # with no network the run must still finish, every connection attempt must be local and the
        # metrics payload must stay on disk (its POST to the closed port failed)
        tmp=$(mktemp -d)
        printf 'sample\tchrom\tpos\tref\talt\tvtype\telement_id\nS1\tchr5\t1295113\tG\tA\tSNP\tpromoter_ENSR00000746745\nS2\tchr5\t1295135\tG\tA\tSNP\tpromoter_ENSR00000746745\nS1\tchr1\t1000000\tT\tA\tSNP\telem_B\n' > "$tmp/cand.tsv"
        awk -F '\t' -v OFS='\t' -f /pipeline/assets/awk/opencravat_input.awk "$tmp/cand.tsv" > "$tmp/in.txt"
        mkdir "$tmp/metrics" "$tmp/work"
        sed -e "s|^metrics_dir:.*|metrics_dir: $tmp/metrics|" -e "s|^work_dir:.*|work_dir: $tmp/work|" "$OPENCRAVAT_IMAGE_CONF" > "$tmp/cravat-system.yml"
        (cd "$tmp" && env OPENCRAVAT_SYSTEM_CONF_PATH="$tmp/cravat-system.yml" PYTHONPATH="$OPENCRAVAT_NETAUDIT_DIR" OC_NETAUDIT_LOG="$tmp/net.tsv" \
            http_proxy=http://127.0.0.1:9 https_proxy=http://127.0.0.1:9 HTTP_PROXY=http://127.0.0.1:9 HTTPS_PROXY=http://127.0.0.1:9 \
            oc run in.txt -l hg38 -n smoke -a ccre_screen repeat clinvar -t text -d out > run.log 2>&1) || { cat "$tmp/run.log" >&2; exit 1; }
        grep -q '^Finished normally' "$tmp/run.log" || { cat "$tmp/run.log" >&2; echo 'oc run did not finish normally' >&2; exit 1; }
        if awk -F '\t' '$3 !~ /^(127\.0\.0\.1|::1|localhost)$/ { bad = 1 } END { exit bad ? 0 : 1 }' "$tmp/net.tsv"; then
            cat "$tmp/net.tsv" >&2; echo 'OpenCRAVAT tried to reach a non-local address' >&2; exit 1
        fi
        grep -q $'connect\t127.0.0.1\t9' "$tmp/net.tsv" || { echo 'no metrics POST to the closed local port was recorded' >&2; exit 1; }
        test "$(find "$tmp/metrics" -name 'ocmetric_*.json' | wc -l)" -ge 1
        test -z "$(ls -A /opt/opencravat/metrics)"
        # no liftover anywhere: no UCSC chain file in the image and no hg19 column in the report
        if find / -xdev \( -name '*.over.chain' -o -name '*.over.chain.gz' \) 2>/dev/null | grep -q .; then echo 'UCSC chain file found in the image' >&2; exit 1; fi
        if grep -v '^#' "$tmp/out/smoke.tsv" | head -2 | grep -q 'hg19'; then echo 'hg19 columns in the OpenCRAVAT report' >&2; exit 1; fi
        # the adapter joins the report back: TERT C228T in a SCREEN promoter-like cCRE
        awk -F '\t' -v OFS='\t' -f /pipeline/assets/awk/opencravat_table.awk "$tmp/cand.tsv" "$tmp/out/smoke.tsv" > "$tmp/table.tsv"
        awk -F '\t' 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; if (index($0, "hg19")) exit 1; next }
            $c["variant_id"] == "chr5:1295113:G:A" && $c["oc_status"] == "annotated" && $c["ccre_screen__classification"] == "PLS" { ok = 1 }
            END { exit ok ? 0 : 1 }' "$tmp/table.tsv" || { cat "$tmp/table.tsv" >&2; echo 'OpenCRAVAT table check failed' >&2; exit 1; }
        ;;
    *)
        echo "Unknown image: ${name}" >&2
        exit 2
        ;;
esac
echo "smoke test passed: ${name}"
