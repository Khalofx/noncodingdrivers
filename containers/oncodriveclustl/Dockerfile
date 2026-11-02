# Image for ONCODRIVECLUSTL: OncodriveCLUSTL 1.2.0 (PyPI; no Bioconda or conda-forge package) on
# python 3.12 (conda-forge), installed from the two explicit locks in this folder:
#   conda-lock.linux-64.txt  python, pip, setuptools, wheel, pycurl, procps-ng (every package, URL, MD5)
#   requirements.lock        oncodriveclustl==1.2.0 and its dependencies, sha256-hashed (pip --require-hashes)
# Licence: OncodriveCLUSTL is AGPL-3.0-or-later (IRB Barcelona); no non-commercial clause (checked
# 2026-10-29). The upstream LICENSE has no licence text, so the full AGPL-3.0 text ships in
# /usr/share/licenses/oncodriveclustl with NOTICE (the one modification and the corresponding source).
# Modification: the COSMIC Cancer Gene Census list bundled in the wheel is replaced by its header
# line (no COSMIC content is redistributed; with --genome hg38 the CGC column is "Non Available"
# anyway). The recipe files in this folder are MIT, like the pipeline. The same locks build the
# Apptainer image (oncodriveclustl.def).
FROM mambaorg/micromamba:2.0.5

# Tools on PATH without the entrypoint (Apptainer exec and `--entrypoint bash` skip activation).
ENV PATH=/opt/conda/bin:$PATH

LABEL org.opencontainers.image.title="noncodingdrivers-oncodriveclustl" \
      org.opencontainers.image.version="1.2.0" \
      org.opencontainers.image.source="https://github.com/Khalofx/noncodingdrivers" \
      org.opencontainers.image.licenses="AGPL-3.0-or-later AND MIT"

COPY --chown=$MAMBA_USER:$MAMBA_USER conda-lock.linux-64.txt requirements.lock AGPL-3.0.txt NOTICE /tmp/build/
RUN micromamba install -y -n base -f /tmp/build/conda-lock.linux-64.txt \
    && micromamba clean -a -y
ARG MAMBA_DOCKERFILE_ACTIVATE=1
# --no-deps: requirements.lock lists every package; --no-build-isolation: the sdist-only bg*
# packages build with the locked setuptools/wheel, not with whatever an isolated build would fetch
RUN pip install --no-cache-dir --no-deps --no-build-isolation --require-hashes -r /tmp/build/requirements.lock \
    && pip check

USER root
RUN cgc=$(python -c "import os, oncodriveclustl; print(os.path.join(os.path.dirname(oncodriveclustl.__file__), 'data', 'CGC_all_Oct15_10_29_09_2018.tsv'))") \
    && printf 'Gene Symbol\tTier\n' > "$cgc" \
    && mkdir -p /usr/share/licenses/oncodriveclustl \
    && cp /tmp/build/AGPL-3.0.txt /tmp/build/NOTICE /usr/share/licenses/oncodriveclustl/ \
    && rm -rf /tmp/build
USER $MAMBA_USER

RUN oncodriveclustl --help > /dev/null \
    && python -c "import importlib.metadata as m; v = m.version('oncodriveclustl'); assert v == '1.2.0', v" \
    && ps --version
