#!/usr/bin/env bash
# Run inside the selected image with the repository mounted read-only at /pipeline:
#   docker run --rm --network none --entrypoint bash -v "$PWD:/pipeline:ro" <image> \
#       /pipeline/containers/smoke-test.sh <name>
set -euo pipefail

name="${1:?Specify utils, activedriverwgs, mutenricher, oncodrivefml or combinepvalues}"

# Nextflow's task wrapper runs ps inside every task container to collect CPU and memory metrics.
command -v ps > /dev/null || { echo "ps (procps) is missing from the ${name} image" >&2; exit 1; }

case "${name}" in
    utils)
        python -c 'import pandas, numpy, scipy.stats'
        ;;
    activedriverwgs)
        Rscript -e "library(ActiveDriverWGS); library(BSgenome.Hsapiens.UCSC.hg38); \
            stopifnot(as.character(packageVersion('ActiveDriverWGS')) == '1.2.1'); \
            stopifnot(GenomeInfoDb::seqlengths(BSgenome.Hsapiens.UCSC.hg38::Hsapiens)[['chr1']] == 248956422)"
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
    combinepvalues)
        Rscript -e "library(EmpiricalBrownsMethod); library(RobustRankAggreg); library(ACAT)"
        ;;
    *)
        echo "Unknown image: ${name}" >&2
        exit 2
        ;;
esac
echo "smoke test passed: ${name}"
