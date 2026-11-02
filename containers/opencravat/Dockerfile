# Image for OPENCRAVAT: OpenCRAVAT 3.2.4 (MIT) on Python 3.11 with gawk and procps, the OpenCRAVAT
# base modules and the annotators encode_tfbs, ccre_screen, repeat and clinvar, installed at build
# time at the versions of modules.freeze.json. Telemetry is neutralised: the run-time system
# configuration (cravat-system.yml) points the metrics, store and publish addresses at a closed
# local port, the metrics folder is empty and the email prompt is off; netaudit/ records every
# connection attempt when the OPENCRAVAT module enables it. Module data licences: README.md.
# The same locks build the Apptainer image (opencravat.def).
FROM mambaorg/micromamba:2.0.5

ENV PATH=/opt/conda/bin:$PATH \
    OPENCRAVAT_IMAGE_CONF=/opt/conda/lib/python3.11/site-packages/cravat/conf/cravat-system.yml \
    OPENCRAVAT_NETAUDIT_DIR=/opt/opencravat/netaudit \
    PYTHONNOUSERSITE=1

LABEL org.opencontainers.image.title="noncodingdrivers-opencravat" \
      org.opencontainers.image.version="3.2.4" \
      org.opencontainers.image.source="https://github.com/Khalofx/noncodingdrivers" \
      org.opencontainers.image.description="OpenCRAVAT 3.2.4 with base modules, encode_tfbs, ccre_screen, repeat and clinvar; telemetry to a closed local port" \
      org.opencontainers.image.licenses="MIT"

USER root
COPY conda-lock.linux-64.txt pip-lock.txt modules.freeze.json cravat-system.yml setup-opencravat.sh /opt/opencravat/build/
COPY netaudit/sitecustomize.py /opt/opencravat/netaudit/sitecustomize.py
RUN micromamba install -y -p /opt/conda -f /opt/opencravat/build/conda-lock.linux-64.txt \
    && micromamba clean -a -y
RUN pip install --no-cache-dir --no-deps --require-hashes -r /opt/opencravat/build/pip-lock.txt \
    && pip check
RUN bash /opt/opencravat/build/setup-opencravat.sh
USER $MAMBA_USER
