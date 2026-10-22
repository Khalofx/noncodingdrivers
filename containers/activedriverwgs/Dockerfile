# ActiveDriverWGS 1.2.1 (CRAN) with the matching BSgenome, so no host R library is needed.
# R 4.5.3 and BSgenome.Hsapiens.UCSC.hg38 1.4.5 match the runtime that produced the
# regression baseline (site R 4.5.3, tests/regression/images/activedriverwgs.packages.txt).
# Bioconductor 3.22 is the release for R 4.5.
FROM rocker/r-ver:4.5.3

LABEL org.opencontainers.image.title="noncodingdrivers-activedriverwgs" \
      org.opencontainers.image.version="1.2.1" \
      org.opencontainers.image.source="https://github.com/Khalofx/noncodingdrivers" \
      org.opencontainers.image.licenses="GPL-3.0-only"

# procps: `ps` for Nextflow task metrics.
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates curl procps libcurl4-openssl-dev libssl-dev libxml2-dev \
        libbz2-dev liblzma-dev zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Bioconductor dependencies first (CRAN cannot resolve them), then the exact caller source
# tarball, checked by SHA-256, rather than whatever CRAN serves today.
RUN Rscript -e "options(timeout=3600, Ncpus=2); \
        install.packages('BiocManager', repos='https://cloud.r-project.org'); \
        BiocManager::install(version='3.22', ask=FALSE, update=FALSE); \
        BiocManager::install(c('BSgenome', 'Biostrings', 'GenomeInfoDb', \
          'GenomicRanges', 'IRanges', 'S4Vectors', 'BSgenome.Hsapiens.UCSC.hg38'), \
          ask=FALSE, update=FALSE)" \
    && (curl --fail --location --retry 5 \
        https://cran.r-project.org/src/contrib/ActiveDriverWGS_1.2.1.tar.gz \
        --output /tmp/ActiveDriverWGS.tar.gz \
        || curl --fail --location --retry 5 \
        https://cran.r-project.org/src/contrib/Archive/ActiveDriverWGS/ActiveDriverWGS_1.2.1.tar.gz \
        --output /tmp/ActiveDriverWGS.tar.gz) \
    && echo '1e01413214990b2846c761f212215e2e672323afe4ae81c76017eba9d6b7d4ec  /tmp/ActiveDriverWGS.tar.gz' | sha256sum --check - \
    && R CMD INSTALL /tmp/ActiveDriverWGS.tar.gz \
    && rm /tmp/ActiveDriverWGS.tar.gz

RUN Rscript -e "library(ActiveDriverWGS); library(BSgenome.Hsapiens.UCSC.hg38); \
        stopifnot(as.character(packageVersion('ActiveDriverWGS')) == '1.2.1'); \
        stopifnot(as.character(packageVersion('BSgenome.Hsapiens.UCSC.hg38')) == '1.4.5'); \
        stopifnot(length(BSgenome.Hsapiens.UCSC.hg38::Hsapiens[['chr1']]) == 248956422)" \
    && ps --version

WORKDIR /data
