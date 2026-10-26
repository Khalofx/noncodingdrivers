# Image for MOTIFBREAKR: motifbreakR 2.24.0 + universalmotif 1.28.0 + BSgenome.Hsapiens.UCSC.hg38
# 1.4.5 (Bioconda, R 4.5), installed from the explicit conda lock in this folder.
# Component licences: motifbreakR GPL-2, universalmotif GPL-3, BSgenome packages Artistic-2.0;
# the UCSC hg38 sequence is freely redistributable. The recipe itself is MIT.
# The same lock builds the Apptainer image (motifbreakr.def).
FROM mambaorg/micromamba:2.0.5

# Tools on PATH without the entrypoint (Apptainer exec and `--entrypoint bash` skip activation).
ENV PATH=/opt/conda/bin:$PATH

LABEL org.opencontainers.image.title="noncodingdrivers-motifbreakr" \
      org.opencontainers.image.version="2.24.0" \
      org.opencontainers.image.source="https://github.com/Khalofx/noncodingdrivers" \
      org.opencontainers.image.licenses="MIT"

COPY --chown=$MAMBA_USER:$MAMBA_USER conda-lock.linux-64.txt /tmp/conda-lock.txt
RUN micromamba install -y -n base -f /tmp/conda-lock.txt \
    && micromamba clean -a -y \
    && rm /tmp/conda-lock.txt
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN Rscript -e "library(motifbreakR); library(universalmotif); library(BSgenome.Hsapiens.UCSC.hg38); \
      stopifnot(as.character(packageVersion('motifbreakR')) == '2.24.0', \
                as.character(packageVersion('universalmotif')) == '1.28.0', \
                as.character(packageVersion('BSgenome.Hsapiens.UCSC.hg38')) == '1.4.5')" \
    && ps --version
