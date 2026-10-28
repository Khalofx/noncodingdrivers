# khalofx/noncodingdrivers

[![GitHub Actions CI Status](https://github.com/khalofx/noncodingdrivers/actions/workflows/nf-test.yml/badge.svg)](https://github.com/khalofx/noncodingdrivers/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/khalofx/noncodingdrivers/actions/workflows/linting.yml/badge.svg)](https://github.com/khalofx/noncodingdrivers/actions/workflows/linting.yml)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)
[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.10.4-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.1.0-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.1.0)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/khalofx/noncodingdrivers)

## Introduction

**khalofx/noncodingdrivers** finds candidate non-coding driver elements in cohorts of cancer
genomes and reports the evidence for and against each one. It starts from pre-called somatic
variants on GRCh38, either VCF (including nf-core/sarek output) or MAF. It does not call variants.

For every element of a non-coding catalogue, the pipeline asks whether the cohort carries more
mutations, or more high-impact mutations, than expected. It reports per-element p- and q-values.
Beside them it reports the flags that can explain a signal away, such as low mappability, germline
leakage or mutations concentrated in one patient, and the annotation of each candidate's
mutations. No step filters a mutation or re-ranks an element; the flags and annotations are
added files.

1. Input harmonisation: an assembly check of every row against its declared build (EBI
   vcf_assembly_checker), contig renaming, MAF conversion (maf2vcf.pl), an opt-in GRCh37/hg19 to
   GRCh38 liftover (Picard LiftoverVcf), normalisation and PASS filtering (bcftools), and one
   mutation table per cohort.
2. Driver discovery on a CDS-excluded catalogue of Ensembl Regulatory Build elements:
   ActiveDriverWGS (default) and OncodriveFML (opt-in, licence-gated, with user-supplied CADD
   scores); MutEnricher as a QC-only caller.
3. Combination of the caller p-values with the empirical Brown's method (default) or ACAT, with
   Benjamini–Hochberg q-values; ACAT, robust rank aggregation and caller votes are kept as
   sensitivity columns.
4. Variant flags (on by default): ENCODE Blacklist v2, GIAB v3.6 segmental duplications and low
   mappability, Umap k100 uniqueness and gnomAD v4.1 genomes, per mutation and per element, and
   each element's callable fraction.
5. PCAWG-style post-filter flags and a long-form evidence table per cohort.
6. Optional null calibration: the callers re-run on neutral copies of the cohort.
7. Annotation of the significant elements: their mutations per patient and the nearest gene, plus
   the optional ENCODE regulatory context, Ensembl VEP and motifbreakR. SigProfiler signature QC
   is optional.
8. A MultiQC run summary.

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/get_started/environment_setup/overview) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/get_started/run-your-first-pipeline) with `-profile test` before running the workflow on actual data.

Prepare a samplesheet with one row per input file. Rows with the same `cohort` are merged into
one cohort:

```csv title="samplesheet.csv"
sample,entry,cohort,cancer_type,assembly,mutations,vcf,index,maf,vcf_tumour_id
D01,vcf,LIVER,HCC,GRCh38,,D01.vcf.gz,D01.vcf.gz.tbi,,
D02,sarek,LIVER,HCC,GRCh38,,D02.mutect2.filtered.vcf.gz,D02.mutect2.filtered.vcf.gz.tbi,,D02_T
pcawg_liver,maf,LIVER,HCC,GRCh38,,,,liver.maf.gz,
```

Then run the pipeline with the GRCh38 reference set of the datasheet
`assets/references/GRCh38.yml` (see [references](docs/references.md)):

```bash
nextflow run khalofx/noncodingdrivers -r 1.0.0 \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --reference_base /path/to/references \
   --outdir <OUTDIR>
```

OncodriveFML runs only with `--run_oncodrivefml true --accept_oncodrivefml_licence true` and a
CADD v1.7 score file (`--oncodrivefml_score`). Its licence and CADD's restrict commercial use, and
its image must be built locally ([usage](docs/usage.md#callers-and-licences)).

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).

Two test profiles are provided: `test` runs a small synthetic TERT-window cohort and needs no
other parameter, and `test_full` runs SYNTH-NC, a synthetic genome-wide cohort with five planted
drivers, on the full-size references. [Usage](docs/usage.md) describes every option and
[output](docs/output.md) every result file. [Validation](docs/validation.md) records the
regression runs of this release.

## Credits

khalofx/noncodingdrivers was originally written by Khalifa Juma.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](docs/CONTRIBUTING.md).

## Citations

If you use khalofx/noncodingdrivers for your analysis, please cite the repository
(<https://github.com/Khalofx/noncodingdrivers>) and the release you used, and the tools and data
listed in [`CITATIONS.md`](CITATIONS.md). The design of the element-level tests and post-filters
follows the PCAWG non-coding driver analysis (Rheinbay et al. 2020, _Nature_ 578:102–111,
doi:[10.1038/s41586-020-1965-x](https://doi.org/10.1038/s41586-020-1965-x)).

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
