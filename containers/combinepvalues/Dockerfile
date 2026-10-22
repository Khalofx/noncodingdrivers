# Image for COMBINEPVALUES: R 4.4.3 + EmpiricalBrownsMethod 1.34.0 (Bioconda) +
# RobustRankAggreg 1.2.1 (conda-forge) + ACAT 0.91 from the authors' repository at a pinned commit.
# ACAT is GPL-3: it is installed into the image, never copied into this MIT repository.
FROM mambaorg/micromamba:2.0.5

# Tools on PATH without the entrypoint (Apptainer exec and `--entrypoint bash` skip activation).
ENV PATH=/opt/conda/bin:$PATH

LABEL org.opencontainers.image.title="noncodingdrivers-combinepvalues" \
      org.opencontainers.image.version="1.0.0" \
      org.opencontainers.image.source="https://github.com/Khalofx/noncodingdrivers" \
      org.opencontainers.image.licenses="MIT"

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/env.yml
RUN micromamba install -y -n base -f /tmp/env.yml \
    && micromamba install -y -n base -c conda-forge r-matrix git procps-ng \
    && micromamba clean -a -y
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN git clone https://github.com/yaowuliu/ACAT /tmp/ACAT \
    && git -C /tmp/ACAT checkout e737d7d \
    && R CMD INSTALL /tmp/ACAT && rm -rf /tmp/ACAT \
    && Rscript -e "library(EmpiricalBrownsMethod); library(RobustRankAggreg); library(ACAT); \
         stopifnot(as.character(packageVersion('EmpiricalBrownsMethod')) == '1.34.0')" \
    && ps --version
