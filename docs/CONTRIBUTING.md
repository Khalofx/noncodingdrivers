---
title: Contributing
markdownPlugin: checklist
---

# `khalofx/noncodingdrivers`: Contributing guidelines

Hi there!
Thanks for taking an interest in improving khalofx/noncodingdrivers.

This page describes the recommended nf-core way to contribute to both khalofx/noncodingdrivers and nf-core pipelines in general, including:

- [General contribution guidelines](#general-contribution-guidelines): common procedures or guides across all nf-core pipelines.
- [Pipeline-specific contribution guidelines](#pipeline-specific-contribution-guidelines): procedures or guides specific to the development conventions of khalofx/noncodingdrivers.

## General contribution guidelines

### Contribution quick start

To contribute code to any nf-core pipeline:

- [ ] Ensure you have Nextflow, nf-core tools, and nf-test installed. See the [nf-core/tools repository](https://github.com/nf-core/tools) for instructions.
- [ ] Check whether a GitHub [issue](https://github.com/khalofx/noncodingdrivers/issues) about your idea already exists. If an issue does not exist, create one so that others are aware you are working on it.
- [ ] [Fork](https://help.github.com/en/github/getting-started-with-github/fork-a-repo) the [khalofx/noncodingdrivers repository](https://github.com/khalofx/noncodingdrivers) to your GitHub account.
- [ ] Create a branch on your forked repository and make your changes following [pipeline conventions](#pipeline-contribution-conventions) (if applicable).
- [ ] To fix major bugs, name your branch `patch` and follow the [patch release](#patch-release) process.
- [ ] Update relevant documentation within the `docs/` folder, use nf-core/tools to update `nextflow_schema.json`, and update `CITATIONS.md`.
- [ ] Run and/or update tests. See [Testing](#testing) for more information.
- [ ] [Lint](#lint-tests) your code with nf-core/tools.
- [ ] Submit a pull request (PR) against the `dev` branch and request a review.

If you are not used to this workflow with Git, see the [GitHub documentation](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests) or [Git resources](https://try.github.io/) for more information.

### GitHub Codespaces

You can contribute to khalofx/noncodingdrivers without installing a local development environment on your machine by using [GitHub Codespaces](https://github.com/codespaces).

[GitHub Codespaces](https://github.com/codespaces) is an online developer environment that runs in your browser, complete with VS Code and a terminal.
Most nf-core repositories include a devcontainer configuration, which creates a GitHub Codespaces environment specifically for Nextflow development.
The environment includes pre-installed nf-core tools, Nextflow, and a few other helpful utilities via a Docker container.

To get started, open the repository in [Codespaces](https://github.com/khalofx/noncodingdrivers/codespaces).

### Testing

Once you have made your changes, run the pipeline with nf-test to test them locally.
For additional information, use the `--verbose` flag to view the Nextflow console log output.

```bash
nf-test test --tag test --profile +docker --verbose
```

If you have added new functionality, ensure you update the test assertions in the `.nf.test` files in the `tests/` directory.
Update the snapshots with the following command:

```bash
nf-test test --tag test --profile +docker --verbose --update-snapshots
```

When you create a pull request with changes, GitHub Actions will run automatic tests.
Pull requests are typically reviewed when these tests are passing.

Two types of tests are typically run:

#### Lint tests

nf-core has a [set of guidelines](https://nf-co.re/docs/specifications/overview) which all pipelines must follow.
To enforce these, run linting with nf-core/tools:

```bash
nf-core pipelines lint <pipeline_directory>
```

If you encounter failures or warnings, follow the linked documentation printed to screen.
For more information about linting tests, see [nf-core/tools API documentation](https://nf-co.re/docs/nf-core-tools/api_reference/latest/pipeline_lint_tests/actions_awsfulltest).

#### Pipeline tests

Each nf-core pipeline should be set up with a minimal set of test data.
GitHub Actions runs the pipeline on this data to ensure it runs through and exits successfully.
If there are any failures then the automated tests fail.
These tests are run with the latest available version of Nextflow and the minimum required version specified in the pipeline code.

### Patch release

> [!WARNING]
> Only in the unlikely event of a release that contains a critical bug.

- [ ] Create a new branch `patch` on your fork based on `upstream/main` or `upstream/master`.
- [ ] Fix the bug and use nf-core/tools to bump the version to the next semantic version, for example, `1.2.3` → `1.2.4`.
- [ ] Open a Pull Request from `patch` directly to `main`/`master` with the changes.

### Pipeline contribution conventions

nf-core semi-standardises how you write code and other contributions to make the khalofx/noncodingdrivers code and processing logic more understandable for new contributors and to ensure quality.

#### Add a new pipeline step

To contribute a new step to the pipeline, follow the general nf-core coding procedure.
Please also refer to the [pipeline-specific contribution guidelines](#pipeline-specific-contribution-guidelines):

- [ ] Define the corresponding [input channel](#channel-naming-schemes) into your new process from the expected previous process channel.
- [ ] Install a module with nf-core/tools, or write a local module (see [default processes resource requirements](#default-processes-resource-requirements)), and add it to the target `<workflow>.nf`.
- [ ] Define the output channel if needed. Mix the version output channel into `ch_versions` and relevant files into `ch_multiqc`.
- [ ] Add new or updated parameters to `nextflow.config` with a [default value](#default-parameter-values).
- [ ] Add new or updated parameters and relevant help text to `nextflow_schema.json` with [nf-core/tools](#default-parameter-values).
- [ ] Add validation for relevant parameters to the pipeline utilisation section of `utils_nfcore_\_pipeline/main.nf` subworkflow.
- [ ] Perform local tests to validate that the new code works as expected.
  - [ ] If applicable, add a new test in the `tests` directory.
- [ ] Update `usage.md`, `output.md`, and `citation.md` as appropriate.
- [ ] [Lint](#lint-tests) the code with nf-core/tools.
- [ ] Update any diagrams or pipeline images as necessary.
- [ ] Update MultiQC config `assets/multiqc_config.yml` so relevant suffixes, file name cleanup, and module plots are in the appropriate order.
- [ ] If applicable, create a [MultiQC](https://seqera.io/multiqc/) module.
- [ ] Add a description of the output files and, if relevant, images from the MultiQC report to `docs/output.md`.

To update the minimum required Nextflow version, see the [Nextflow version bumping](#nextflow-version-bumping) section below. For more information about pipeline contributions, see [pipeline-specific contribution guidelines](#pipeline-specific-contribution-guidelines).

#### Channel naming schemes

Use the following naming schemes for channels to make the channel flow easier to understand:

- Initial process channel: `ch_output_from_<process>`
- Intermediate and terminal channels: `ch_<previousprocess>_for_<nextprocess>`

#### Default parameter values

Parameters should be initialised and defined with default values within the `params` scope in `nextflow.config`.
They should also be documented in the pipeline JSON schema.

To update `nextflow_schema.json`, run:

```bash
nf-core pipelines schema build
```

The schema builder interface that loads in your browser should automatically update the defaults in the parameter documentation.

#### Default processes resource requirements

If you write a local module, specify a default set of resource requirements for the process.

Sensible defaults for process resource requirements (CPUs, memory, time) should be defined in `conf/base.config`.
Specify these with generic `withLabel:` selectors, so they can be shared across multiple processes and steps of the pipeline.

nf-core provides a set of standard labels that you should follow where possible, as seen in the [nf-core pipeline template](https://github.com/nf-core/tools/blob/main/nf_core/pipeline-template/conf/base.config).
These labels define resource defaults for single-core processes, modules that require a GPU, and different levels of multi-core configurations with increasing memory requirements.

Values assigned within these labels can be dynamically passed to a tool using the `${task.cpus}` and `${task.memory}` Nextflow variables in the `script:` block of a module (see an example in the [modules repository](https://github.com/nf-core/modules/blob/bd1b6a40f55933d94b8c9ca94ec8c1ea0eaf4b82/modules/nf-core/samtools/bam2fq/main.nf#L30)).

#### Nextflow version bumping

If you use a new feature from core Nextflow, bump the minimum required Nextflow version in the pipeline with:

```bash
nf-core pipelines bump-version --nextflow . <min_nf_version>
```

#### Images and figures guidelines

If you update images or graphics, follow the nf-core [style guidelines](https://nf-co.re/docs/community/brand/workflow-schematics).

## Pipeline specific contribution guidelines

## Module standard (noncodingdrivers)

Every process follows one recipe, so later steps are mechanical. The full rules are R0–R10 in the build guide, Step 2.

1. **Choose the component.** Prefer an installed nf-core module (`nf-core modules install <tool/subtool>`), configured only through `ext.args` / `ext.prefix` in `conf/modules.config`. Otherwise write a local module around a packaged tool (Bioconda, conda-forge or CRAN), reusing the container of the closest nf-core module. Glue between tools is a format adapter (an inline `awk` or an aliased nf-core `gawk` with `ext.args2`) that never computes a statistic. Anything else is pipeline-owned code and needs an entry in the provenance exception table.
2. **Layout and naming.** Local modules are flat, `modules/local/<name>/{main.nf,meta.yml,environment.yml,tests/}`, `<name>` without punctuation, process `<NAME>`. Aliased nf-core modules get descriptive include names (`GAWK as COHORT_DEDUP`).
3. **Versions through topic channels.** One output per tool: `tuple val("${task.process}"), val('<tool>'), eval('<cmd>'), topic: versions, emit: versions_<tool>`, each listed under `topics: versions:` in `meta.yml`. No `versions.yml`. Version commands avoid `$` and double quotes; write a `<tool>.version.txt` and `cat` it when needed. MultiQC never emits on the topic.
4. **Inputs and options.** `tuple val(meta), path(x)`; references get `meta2`, `meta3`. Optional files are `[]`. Modules read no `params` and no non-standard meta keys; options travel in `task.ext.args`, output names in `task.ext.prefix`.
5. **Software.** `conda "${moduleDir}/environment.yml"` with channel- and version-pinned, build-free specs. No conda package means no conda directive; `environment.yml` then documents the runtime.
6. **Containers.** Pipeline-owned images are built by `.github/workflows/build-containers.yml` from `containers/<name>/Dockerfile` and published as `ghcr.io/khalofx/noncodingdrivers-<name>:<version>`. Docker Hub images are written with their registry (`docker.io/...`). Every image has `ps` (procps). Digests are pinned before release.
7. **Scripts.** `bin/` holds only the adapters listed in the guide; each has a shebang, is executable and is called by name from PATH. No `${projectDir}` inside a process.
8. **Stub, when, labels.** A `stub:` block that creates every output; `when: task.ext.when == null || task.ext.when`; exactly one standard label.
9. **Strict syntax.** Lower-case `channel` factories, `def` for workflow-local variables, explicit closure parameters (never `it`).
10. **Tests.** `tests/main.nf.test` tagged `modules`, `modules_local`, `<name>`, with a stub test that snapshots `process.out` and a real test on committed fixtures. Float-heavy outputs are asserted in Groovy (header, row count, planted elements); only versions and file names are snapshotted.

Checks before every commit:

```bash
nextflow lint .
nf-core modules lint --all --local --registry ghcr.io,quay.io,community.wave.seqera.io/library,docker.io
nf-test test modules/local/<name>/tests/main.nf.test --profile +singularity
pre-commit run --all-files
```
