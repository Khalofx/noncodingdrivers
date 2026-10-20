# MutEnricher 1.3.3 (QC-only caller in noncodingdrivers v1.0; mutation-burden family).
#
# Licence: MIT. Verified by reading LICENSE.txt in the upstream repo verbatim; the
# GitHub API reports NOASSERTION only because a copyright line precedes the licence
# text. The image is therefore redistributable.
#
# TWO BUILD CONSTRAINTS, both learned by hitting them:
#
# 1. python:3.7 is REQUIRED, not conservative pinning. MutEnricher's setup.py builds a
#    Cython extension (math_funcs) against the interpreter's C API, and the upstream
#    pins (numpy 1.17.3 / scipy 1.3.1) predate the numpy 2.x ABI. Newer Python or newer
#    numpy breaks the extension at import, not at build.
#
# 2. cyvcf2 MUST be built from source, and the flag must be scoped to cyvcf2 alone.
#    Its prebuilt wheel is compiled against a different numpy ABI and segfaults on
#    import next to numpy 1.17.3. `--no-binary :all:` recurses into build dependencies
#    and deadlocks pip's build tracker, so scope it: `--no-binary cyvcf2`.
#
# VERIFYING THIS IMAGE: math_funcs ships without an __init__.py, so
# `python -c "from math_funcs import WAP"` fails from any cwd other than the module
# directory. That is not a broken build — smoke-test through the entry point instead:
#     python /opt/MutEnricher/mutEnricher.py --version
#
# The full buildpack base supplies the legacy compiler/development libraries.
# This avoids installing from the retired Bullseye security package repository.
FROM python:3.7.17-bullseye

ARG MUTENRICHER_COMMIT=ddb1863d32c248aeec8797004d13193d836573a1

LABEL org.opencontainers.image.title="noncodingdrivers-mutenricher" \
      org.opencontainers.image.version="1.3.3" \
      org.opencontainers.image.licenses="MIT" \
      org.opencontainers.image.source="https://github.com/Khalofx/noncodingdrivers" \
      org.opencontainers.image.revision="${MUTENRICHER_COMMIT}"

RUN curl --fail --location --retry 5 \
        https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2 \
        -o /tmp/htslib.tar.bz2 \
    && echo '84b510e735f4963641f26fd88c8abdee81ff4cb62168310ae716636aac0f1823  /tmp/htslib.tar.bz2' | sha256sum --check - \
    && tar -xjf /tmp/htslib.tar.bz2 -C /tmp \
    && cd /tmp/htslib-1.21 && ./configure --disable-libcurl \
    && make -j2 bgzip tabix \
    && install bgzip tabix /usr/local/bin/ \
    && rm -rf /tmp/htslib-1.21 /tmp/htslib.tar.bz2

# Keep the legacy extension build tools compatible with Python 3.7 and pysam.
# bgzip/tabix above are standalone executables, not Python bindings, and do not
# change the numpy ABI. MutEnricher reads bgzipped, tabix-indexed per-sample VCFs.
RUN pip install --no-cache-dir pip==23.1.2 setuptools==57.5.0 wheel==0.38.4 \
    && pip install --no-cache-dir numpy==1.17.3 scipy==1.3.1 cython==0.29.13 pysam==0.15.2 \
    && pip install --no-cache-dir --no-build-isolation --no-binary cyvcf2 cyvcf2==0.20.0

RUN cd /opt \
    && curl --fail --location --retry 5 https://github.com/asoltis/MutEnricher/archive/${MUTENRICHER_COMMIT}.tar.gz -o me.tar.gz \
    && tar -xzf me.tar.gz && rm me.tar.gz \
    && mv MutEnricher-${MUTENRICHER_COMMIT} MutEnricher \
    && cd /opt/MutEnricher/math_funcs && python setup.py build_ext --inplace

# Fail the build rather than ship an image whose extension does not load.
RUN apt-get update && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*

RUN cd /opt/MutEnricher && python mutEnricher.py --version \
    && bgzip --version && tabix --version && ps --version

ENV MUTENRICHER=/opt/MutEnricher
WORKDIR /data
