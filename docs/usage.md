# khalofx/noncodingdrivers: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

The pipeline tests every element of a non-coding element catalogue for an excess of somatic
mutations in a cohort and explains the significant elements. Its inputs are pre-called somatic
variants on GRCh38 (VCF, nf-core/sarek output or MAF) and a reference set described by a datasheet.
Nothing in the pipeline filters a mutation or changes a p-value after the combination step: the
variant flags, post-filter flags, null calibration and annotation layers add files beside the
consensus table.

## Samplesheet input

The samplesheet is a comma-separated file with a header row. Each row is one input file of a
cohort. Rows with the same `cohort` are merged into one cohort mutation table before the callers
run.

```bash
--input '[path to samplesheet file]'
```

| Column          | Description                                                                                                                             |
| --------------- | --------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`        | Row identifier (`[A-Za-z0-9._-]`). For `vcf` and `sarek` rows it is also the sample name in the cohort table.                           |
| `entry`         | `vcf`, `sarek` or `maf` (below).                                                                                                        |
| `cohort`        | Cohort the row belongs to; empty makes the row its own cohort.                                                                          |
| `cancer_type`   | Optional label.                                                                                                                         |
| `assembly`      | `GRCh38` (default) or `hg38`; `GRCh37` or `hg19` for rows lifted to GRCh38 (with `--allow_hg19_liftover`).                              |
| `vcf`, `index`  | `vcf` and `sarek` rows: bgzipped VCF (`.vcf.gz`) and optional `.tbi`/`.csi`.                                                            |
| `maf`           | `maf` rows: MAF (`.maf` or `.maf.gz`), one or many tumours.                                                                             |
| `vcf_tumour_id` | Tumour column of a multi-sample VCF; required for `sarek` rows.                                                                         |
| `mutations`     | Internal: the canonical-table tier used by the test profile and the regression gates. It is not supported for analyses; leave it empty. |

### Entry tiers

- **`vcf`**: one sample's somatic calls. A VCF with several sample columns needs `vcf_tumour_id`.
- **`sarek`**: a tumour/normal VCF (for example nf-core/sarek output); `vcf_tumour_id` names the
  tumour column.
- **`maf`**: a MAF with any number of tumours (`Tumor_Sample_Barcode`); each tumour becomes one
  sample. Characters outside `[A-Za-z0-9._-]` in barcodes become `_`.

Every row goes through `INPUT_HARMONISATION`: contig names are mapped to UCSC names, maf rows are
converted with maf2vcf.pl, GRCh37/hg19 rows are lifted with Picard LiftoverVcf, then bcftools
normalises (multiallelic records split, indels left-aligned, records whose REF disagrees with the
FASTA dropped) and keeps PASS or missing-FILTER records on chr1-22, X and Y. The records are
written to the cohort table `sample chrom pos ref alt vtype`: deletions start at their first
deleted base with `alt` `-`, insertions sit on the base before them with `ref` `-`, and `vtype` is
`SNP`, `DNP`, `TNP`, `ONP`, `INS` or `DEL`. The per-sample PASS VCFs are published under
`harmonisation/vcf/` and are what MutEnricher (`--run_mutenricher`) reads.

### Assembly check

Every row is checked against the FASTA of its declared assembly before anything else reads it (EBI
`vcf_assembly_checker`). A row in which fewer than `--assembly_min_match` (default 0.99) of the REF
alleles match stops the run with `ASSEMBLY_MISMATCH`; for maf rows the check reads the MAF SNVs
themselves. A row with fewer than `--assembly_min_records` (default 50) checkable SNVs keeps its
declared assembly with an `ASSEMBLY_LOW_EVIDENCE` warning. Summaries are published under
`harmonisation/assembly_check/`.

GRCh37/hg19 rows stop the run with `ASSEMBLY_NOT_ENABLED` unless `--allow_hg19_liftover true` is
set; they then need the hg19 FASTA (`HG19_FASTA_MISSING` otherwise) and the chain (see
[references](references.md)). Rejected liftover records are published under
`harmonisation/liftover_rejects/`. Lifted rows are labelled as such in the evidence table
(`native_or_lifted`).

```csv title="samplesheet.csv"
sample,entry,cohort,cancer_type,assembly,mutations,vcf,index,maf,vcf_tumour_id
D01,vcf,LIVER,HCC,GRCh38,,D01.vcf.gz,D01.vcf.gz.tbi,,
D02,sarek,LIVER,HCC,GRCh38,,D02.mutect2.filtered.vcf.gz,D02.mutect2.filtered.vcf.gz.tbi,,D02_T
pcawg_liver,maf,LIVER,HCC,hg19,,,,liver.hg19.maf.gz,
```

## References

The GRCh38 reference set comes from the datasheet `assets/references/GRCh38.yml`: the analysis
FASTA (UCSC names), the CDS-excluded element catalogue, the GENCODE v50 GTF, the Ensembl VEP 116
cache, the ENCODE peak files, the OncodriveFML genome, the variant-flag masks and the gnomAD
archive, each with its checksum. Set `--reference_base` to the folder holding the files, or give
each file explicitly; an explicit parameter wins over its datasheet value. See
[references](references.md) for the checksums and the recipes that build each asset. CADD is never
in the datasheet: it is supplied by the user (`--oncodrivefml_score`).

## Callers and licences

- **ActiveDriverWGS** (`--run_activedriverwgs`, default `true`): SNVs only, against a +/- 50 kb
  background with trinucleotide context. Elements are run in shards of `--adwgs_shard_size` (5,000)
  and the FDR is recomputed over all shards.
- **OncodriveFML** (`--run_oncodrivefml`, default `false`): functional-impact bias of the element's
  mutations, scored with CADD v1.7 (`--oncodrivefml_score`, bgzipped with a `.tbi`). OncodriveFML
  is AGPL-3.0, and its authors' README adds that it may not be copied or distributed for commercial
  purposes or for redistribution without a licence from Universitat Pompeu Fabra; CADD is free for
  non-commercial use only. The pipeline therefore runs it only with
  `--accept_oncodrivefml_licence true`, and its image is not published: build it from
  `containers/oncodrivefml/` and point the process at it in a config file:

  ```bash
  docker build -t oncodrivefml-local:2.5.0 containers/oncodrivefml
  ```

  ```groovy title="oncodrivefml.config"
  process { withName: 'ONCODRIVEFML' { container = 'oncodrivefml-local:2.5.0' } }
  ```

- **MutEnricher** (`--run_mutenricher`, default `false`): QC only. Its native region table is
  published and never enters the combination; it reads the per-sample VCFs.

The per-caller p-values are combined per element with the empirical Brown's method
(`--combination_method ebm`, default) or ACAT (`acat`); both are always published, and the chosen
one sets `consensus_q` and `is_consensus` (Benjamini–Hochberg, `--fdr_threshold` 0.1). Robust rank
aggregation and the caller votes are sensitivity columns. `--skip_acat true` leaves out ACAT where
the package is unavailable.

### Software licences

| Software                                              | Licence                                                     |
| ----------------------------------------------------- | ----------------------------------------------------------- |
| ActiveDriverWGS                                       | GPL-3.0-only                                                |
| OncodriveFML                                          | AGPL-3.0-or-later, with a non-commercial notice (above)     |
| CADD v1.7 scores (user-supplied)                      | free for non-commercial use; commercial use needs a licence |
| MutEnricher                                           | MIT                                                         |
| EmpiricalBrownsMethod                                 | MIT                                                         |
| ACAT                                                  | GPL-3.0-only                                                |
| RobustRankAggreg                                      | GPL-2.0-only                                                |
| motifbreakR                                           | GPL-2.0-only                                                |
| universalmotif                                        | GPL-3.0-only                                                |
| JASPAR 2026 CORE vertebrates (bundled motifs)         | CC BY 4.0                                                   |
| bcftools, htslib                                      | MIT                                                         |
| bedtools                                              | MIT                                                         |
| csvtk                                                 | MIT                                                         |
| echtvar                                               | MIT                                                         |
| GNU Awk                                               | GPL-3.0-or-later                                            |
| Picard                                                | MIT                                                         |
| vcf2maf (maf2vcf.pl)                                  | Apache-2.0                                                  |
| EBI vcf-validator (vcf_assembly_checker)              | Apache-2.0                                                  |
| Ensembl VEP                                           | Apache-2.0                                                  |
| MultiQC                                               | GPL-3.0-or-later                                            |
| SigProfilerExtractor, SigProfilerMatrixGenerator (QC) | BSD-2-Clause                                                |
| gnomAD v4.1 data                                      | CC0 1.0 (attribution requested)                             |
| ENCODE Blacklist v2                                   | GPL-3.0 (repository licence)                                |
| GIAB v3.6 stratifications                             | U.S. Government work (NIST)                                 |
| Umap k100                                             | GPL                                                         |
| Pipeline code (`bin/`, local modules)                 | MIT                                                         |

## Flags, evidence and null calibration

None of these options filters a mutation or changes a p- or q-value; they add files.

- **Variant flags** (`--run_variant_flags`, default `true`) need the four mask BEDs and, for the
  germline flag, the gnomAD v4.1 genomes echtvar archive and the GRCh38 FASTA. With
  `--reference_base` they come from the datasheet; otherwise give `--blacklist_bed`,
  `--segdup_bed`, `--low_mappability_bed`, `--umap_bed`, `--gnomad_echtvar` and `--fasta`. Without
  a gnomAD archive, `--skip_germline_flags true` writes the germline columns as NA. A run that
  needs none of this sets `--run_variant_flags false`.
- **Post-filter flags and the evidence table** (`--run_driver_evidence`, default `true`).
  `--postfilter_mode pcawg` adds `q_postfilter` beside the unchanged q-values; `--evidence_scope all`
  writes every element instead of the nominal ones. Evidence tiers are assigned downstream from
  `driver_evidence.tsv`.
- **Null calibration** (`--run_null_calibration`, default `false`) re-runs the callers on
  `--null_replicates` neutral copies of each cohort (seed `--seed` + replicate), each costing about
  as much as the observed run, and reports the inflation factor of every caller and combination
  against `--null_lambda_band`.

## Annotation layers

They run on the significant elements only (`is_consensus`); a cohort without one publishes no
candidate or annotation file.

- Candidate mutations per patient on the fragments of each significant element, and the nearest
  gene from the element midpoint to the GENCODE gene TSSs with a strand-aware signed distance
  (always; `--gencode_gtf` from the datasheet, or a pre-built `--gene_tss_bed`).
- ENCODE regulatory context per tissue and layer (`--run_regulatory_context`;
  `--regulatory_peaks_dir`, layout `<dir>/<tissue>/<layer>.<accession>.bed.gz`).
- Ensembl VEP 116.3 on the candidate mutations (`--run_vep`; offline cache `--vep_cache`). The
  `--regulatory` option of VEP is never used.
- motifbreakR on the candidate SNVs (`--run_motifbreakr`; JASPAR 2026 CORE vertebrates bundled,
  `--jaspar_motifs` to replace them).
- SigProfiler cohort-level signature QC (`--run_sigprofiler`; `--sigprofiler_genome` for a
  pre-installed GRCh38 genome on offline nodes).

## Karakoram

`-profile karakoram` (`conf/karakoram.config`) sets SLURM on the `cpu` partition, Apptainer with
the project image cache, the validated reference folder as `--reference_base` and the CADD file as
`--oncodrivefml_score`. ActiveDriverWGS runs with 48 CPUs and 72 GB for up to 120 h. The
per-sample harmonisation tasks run with the local executor inside the head job, so the head job
needs about 10 CPUs and 72 GB. Compute nodes cannot rely on the Nextflow plugin registry: run with
`NXF_OFFLINE=true` and the nf-schema plugin cached, and export `APPTAINER_CACHEDIR` to a writable
folder (the site default is not writable). A digest-pinned image is cached under a new name, so the
private OncodriveFML image needs `apptainer registry login docker://ghcr.io` once per account.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run khalofx/noncodingdrivers --input ./samplesheet.csv --outdir ./results  -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/running/run-pipelines#configuring-pipelines), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run khalofx/noncodingdrivers -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull khalofx/noncodingdrivers
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [khalofx/noncodingdrivers releases page](https://github.com/khalofx/noncodingdrivers/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `test_full`
  - SYNTH-NC, a synthetic genome-wide cohort (60 donors, 231,472 SNVs, five planted drivers and
    one spectrum decoy) from the test-datasets repository
  - Needs the full-size references: `--reference_base` and `--oncodrivefml_score`, or a site
    profile that sets them (`-profile test_full,karakoram`)
- `karakoram`
  - The site profile of the Karakoram cluster (SLURM and Apptainer); see [Karakoram](#karakoram)
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow `24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). ActiveDriverWGS, OncodriveFML and MutEnricher have no conda package, so the callers cannot run under this profile, and ACAT needs `--skip_acat true` until an `r-acat` package is available. Use a container profile.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#set-max-resources) and [customise process resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#customize-process-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#update-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#modifying-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
