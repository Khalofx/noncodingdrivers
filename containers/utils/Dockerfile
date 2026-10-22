# Runtime for the pipeline-owned Python code: RUNSUMMARYMQC and the in-house exceptions
# listed in the provenance section of the build guide (bin/*.py).
#
# bin/ is NOT copied in: Nextflow puts the pipeline's own bin/ on the task PATH, so a copy
# baked into the image would silently diverge from the checked-out code. The image supplies
# the interpreter and libraries only.
#
# Same base and pins as the image that produced the regression baseline (utils 0.2.2,
# tests/regression/images/utils.packages.txt), minus pysam and pyliftover: their only users
# were scripts that v1.0 replaces with established tools. That difference is recorded in
# tests/regression/PROVENANCE.md.
FROM python:3.11.14-slim-bookworm

LABEL org.opencontainers.image.title="noncodingdrivers-utils" \
      org.opencontainers.image.version="1.0.1" \
      org.opencontainers.image.source="https://github.com/Khalofx/noncodingdrivers" \
      org.opencontainers.image.description="Python runtime for noncodingdrivers pipeline-owned scripts (pandas, numpy, scipy)." \
      org.opencontainers.image.licenses="MIT"

# procps supplies `ps`, which Nextflow's task wrapper runs inside the container to record
# CPU and memory use; without it every task fails before its script starts.
RUN apt-get update && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir \
      "pandas==2.3.3" \
      "numpy==1.26.4" \
      "scipy==1.17.1" \
      "python-dateutil==2.9.0.post0" \
      "pytz==2026.3.post1" \
      "tzdata==2026.4"

# Fail the build rather than ship an image the pipeline cannot run in.
RUN python -c "import pandas, numpy, scipy.stats; print('utils deps OK')" && ps --version
