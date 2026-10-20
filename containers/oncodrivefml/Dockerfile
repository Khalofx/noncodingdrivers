# FALLBACK ONLY (build guide Step 2, R5): the pipeline uses the rights-holder image
# docker.io/bbglab/oncodrivefml:2.5.0 unless the S2.3 check fails. This recipe rebuilds the
# image that produced the regression baseline: the upstream image below plus procps.
#
# Nextflow's task wrapper runs `ps` inside the task container to record CPU and
# memory use. The upstream image has no `ps`, so every task fails before its
# script starts: "Command 'ps' required by nextflow to collect task metrics
# cannot be found" (Karakoram, Apptainer, 22 Sep 2026).
#
# The base is pinned by digest, not tag, so OncodriveFML, its Python and its
# dependencies are exactly those of ghcr.io/khan-lab/oncodrivefml:sha-0661d53,
# which produced the regression baseline (tests/regression/images/oncodrivefml.packages.txt). The upstream entrypoint
# (run_oncodrivefml.sh) and its non-root user (biouser) are kept unchanged.
FROM ghcr.io/khan-lab/oncodrivefml@sha256:81f51bd0f8f08f12270221aa9a89405a47281a9fd01bf52e075a3705ee1b3f25

LABEL org.opencontainers.image.title="noncodingdrivers-oncodrivefml" \
      org.opencontainers.image.version="2.5.0" \
      org.opencontainers.image.source="https://github.com/Khalofx/noncodingdrivers" \
      org.opencontainers.image.base.name="ghcr.io/khan-lab/oncodrivefml:sha-0661d53" \
      org.opencontainers.image.base.digest="sha256:81f51bd0f8f08f12270221aa9a89405a47281a9fd01bf52e075a3705ee1b3f25"

USER root
RUN apt-get update && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*
USER biouser

RUN ps --version && oncodrivefml --help > /dev/null
