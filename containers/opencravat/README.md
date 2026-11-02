# noncodingdrivers-opencravat image

Runtime of the `OPENCRAVAT` module: OpenCRAVAT 3.2.4 on Python 3.11 (conda-forge), gawk 5.3.1 for
the format adapters and procps for Nextflow's task metrics, with the OpenCRAVAT modules listed
below installed at build time. The pipeline runs the annotators `ccre_screen`, `repeat` and
`clinvar` on GRCh38 candidate mutations (`oc run -l hg38`).

Size (Apptainer build on Karakoram, 2 November 2026): 4,460 MB uncompressed, of which the
OpenCRAVAT modules are 3,873 MB (clinvar 1,691 MB, hg38wgs 797 MB, hg38 mapper 699 MB, repeat
350 MB); the compressed SIF is 1,586,651,136 bytes.

| File                      | Role                                                                                                                                                                              |
| ------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `environment.yml`         | the requested conda packages (Python 3.11, pip, gawk 5.3.1, procps-ng)                                                                                                            |
| `conda-lock.linux-64.txt` | explicit lock resolved from `environment.yml` (34 packages with URL and MD5)                                                                                                      |
| `pip-lock.txt`            | open-cravat 3.2.4 and its 63 PyPI dependencies, exact versions with the sha256 of each file (`pip install --no-deps --require-hashes`)                                            |
| `modules.freeze.json`     | the OpenCRAVAT modules and versions (output of `oc module freeze`), installed with `oc module install-freeze`                                                                     |
| `cravat-system.yml`       | the run-time OpenCRAVAT system configuration: every remote address is a closed local port                                                                                         |
| `setup-opencravat.sh`     | build step shared by both recipes: installs the frozen modules, removes the UCSC chain files, installs the run-time configuration and checks it                                   |
| `netaudit/`               | `sitecustomize.py`: records every IP connection attempt and host lookup of OpenCRAVAT when the module enables it                                                                  |
| `Dockerfile`              | Docker/OCI image, published by `.github/workflows/build-containers.yml` as `ghcr.io/khalofx/noncodingdrivers-opencravat:3.2.4` after `containers/smoke-test.sh opencravat` passes |
| `opencravat.def`          | Apptainer recipe from the same locks, for building the SIF directly on an HPC host                                                                                                |

## Telemetry

OpenCRAVAT 3.2.4 posts usage metrics after every `oc run` to `metrics_url` + `/job` (operating
system, memory, CPU count, file system, a machine id derived from the MAC address, the configured
email, genome, module names, variant and sample counts, run times). There is no switch that stops
it: `save_metrics` is not read on that path. When the POST fails, OpenCRAVAT keeps the payload in
its metrics folder and re-sends it after the next successful POST.

In this image every remote address of the configuration (`metrics_url`, `store_url`,
`publish_url`, `clingen_api_url`, `hgvs_api_url`) is `http://127.0.0.1:9`, a closed local port, and
the email prompt is off (`user_email_opt_out: true`), so no POST can succeed and no saved payload
is ever re-sent. The `OPENCRAVAT` module adds three safeguards per task: it routes every HTTP(S)
request to the same closed port, records every connection attempt (`netaudit`) and fails if one
is not local, and deletes the payload copy, which it writes to the task folder rather than the
image. The OpenCRAVAT store is reached only while the image is built.

## No liftover, no UCSC chain files

UCSC liftOver chain files are licensed for non-commercial use only. None is in the image: the two
that open-cravat bundles for `-l hg19`/`-l hg18` input (`hg19ToHg38`, `hg18ToHg38`) are deleted
(the module fixes `-l hg38`), and the `hg19` annotator, which ships `hg38ToHg19`, is not
installed. The `encode_tfbs` annotator, whose 2012 ENCODE TFBS data are hg19 coordinates looked up
through that annotator, is therefore not part of the image (decision of 2 November 2026).

## Module licences

The OpenCRAVAT code is MIT (KarchinLab); the module code is MIT
(KarchinLab/open-cravat-modules-karchinlab). Each module's data source and its terms:

| Module                       | Version    | Type            | Data                                                                                                  | Licence                                                                                                           |
| ---------------------------- | ---------- | --------------- | ----------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------- |
| `hg38`                       | 1.47.2     | mapper          | UCSC hg38 gene mapper: GENCODE 47 transcripts (module release note)                                   | GENCODE: released freely for research use (Ensembl/GENCODE); module code MIT                                      |
| `hg38wgs`                    | 1.0.0      | common          | UCSC hg38 whole-genome sequence (GRCh38)                                                              | UCSC data files: free for public and commercial use; module code MIT                                              |
| `ccre_screen`                | 1.0.2      | annotator       | SCREEN registry of candidate cis-regulatory elements, V2 (2019.10)                                    | ENCODE data: "free to download, analyse and publish without restrictions"; module code MIT                        |
| `clinvar`                    | 2026.10.01 | annotator       | NCBI ClinVar (2026.10.01)                                                                             | NCBI: no restrictions on use or distribution (some submitters may claim rights in their records); module code MIT |
| `go`                         | 2025.11.01 | annotator       | Gene Ontology (2025.11.01); installed as a dependency of the GO widgets and never run by the pipeline | CC BY 4.0 (Gene Ontology Consortium); module code MIT                                                             |
| `repeat`                     | 2020.10.16 | annotator       | UCSC RepeatMasker (rmsk) and Simple Repeats tracks, hg38 (2020.10.16)                                 | UCSC table data: free for public and commercial use; module code MIT                                              |
| `clingen-converter`          | 1.0.1      | converter       | none (code only)                                                                                      | MIT (open-cravat-modules-karchinlab)                                                                              |
| `cravat-converter`           | 1.1.2      | converter       | none (code only)                                                                                      | MIT (open-cravat-modules-karchinlab)                                                                              |
| `hgvs-converter`             | 1.1.0      | converter       | none (code only)                                                                                      | MIT (open-cravat-modules-karchinlab)                                                                              |
| `vcf-converter`              | 2.2.1      | converter       | none (code only)                                                                                      | MIT (open-cravat-modules-karchinlab)                                                                              |
| `casecontrol`                | 1.2.0      | postaggregator  | none (code only)                                                                                      | MIT (open-cravat-modules-karchinlab)                                                                              |
| `tagsampler`                 | 1.1.6      | postaggregator  | none (code only)                                                                                      | MIT (open-cravat-modules-karchinlab)                                                                              |
| `varmeta`                    | 1.0.0      | postaggregator  | none (code only)                                                                                      | MIT (open-cravat-modules-karchinlab)                                                                              |
| `vcfinfo`                    | 3.2.0      | postaggregator  | none (code only)                                                                                      | MIT (open-cravat-modules-karchinlab)                                                                              |
| `excelreporter`              | 2.1.1      | reporter        | none (code only)                                                                                      | MIT (open-cravat-modules-karchinlab)                                                                              |
| `textreporter`               | 2.1.0      | reporter        | none (code only)                                                                                      | MIT (open-cravat-modules-karchinlab)                                                                              |
| `wgbase`                     | 1.1.3      | webviewerwidget | none (code only)                                                                                      | MIT (open-cravat-modules-karchinlab)                                                                              |
| `wgcasecontrols`             | 1.0.1      | webviewerwidget | none (code only)                                                                                      | MIT (open-cravat-modules-karchinlab)                                                                              |
| `wgcasecontrolsummary`       | 1.0.1      | webviewerwidget | none (code only)                                                                                      | MIT (open-cravat-modules-karchinlab)                                                                              |
| `wgcircossummary`            | 2.2.0      | webviewerwidget | none (code only)                                                                                      | MIT (open-cravat-modules-karchinlab)                                                                              |
| `wgclinvar`                  | 1.1.1      | webviewerwidget | none (code only)                                                                                      | MIT (open-cravat-modules-karchinlab)                                                                              |
| `wgcodingvsnoncodingsummary` | 2.0.0      | webviewerwidget | none (code only)                                                                                      | MIT (open-cravat-modules-karchinlab)                                                                              |
| `wggo`                       | 1.2.0      | webviewerwidget | none (code only)                                                                                      | MIT (open-cravat-modules-karchinlab)                                                                              |
| `wggosummary`                | 2.4.0      | webviewerwidget | none (code only)                                                                                      | MIT (open-cravat-modules-karchinlab)                                                                              |
| `wglollipop`                 | 2.2.3      | webviewerwidget | none (code only)                                                                                      | MIT (open-cravat-modules-karchinlab)                                                                              |
| `wgncrna`                    | 1.1.0      | webviewerwidget | none (code only)                                                                                      | MIT (open-cravat-modules-karchinlab)                                                                              |
| `wgnote`                     | 3.0.0      | webviewerwidget | none (code only)                                                                                      | MIT (open-cravat-modules-karchinlab)                                                                              |
| `wgsosamplesummary`          | 2.2.0      | webviewerwidget | none (code only)                                                                                      | MIT (open-cravat-modules-karchinlab)                                                                              |
| `wgsosummary`                | 1.5.0      | webviewerwidget | none (code only)                                                                                      | MIT (open-cravat-modules-karchinlab)                                                                              |
| `wgvcfinfo`                  | 1.0.3      | webviewerwidget | none (code only)                                                                                      | MIT (open-cravat-modules-karchinlab)                                                                              |

The Python dependencies in `pip-lock.txt` are open-source packages from PyPI under their own
licences; the CI job records the installed list (`pip freeze`) for every build.

## Build the SIF on Karakoram (or any host with fakeroot and network access)

```bash
cd containers/opencravat
apptainer build --fakeroot noncodingdrivers-opencravat-3.2.4.sif opencravat.def
apptainer test noncodingdrivers-opencravat-3.2.4.sif
```

Nextflow uses a pre-built SIF without pulling when it sits in the Apptainer cache directory under
the name it derives from the module's container.

## Changing a pin

1. Conda: edit `environment.yml` and re-resolve `conda-lock.linux-64.txt` in a micromamba
   container (`micromamba create -p env -f environment.yml`, then
   `micromamba env export -p env --explicit --md5`).
2. PyPI: install `open-cravat==<version>` in that environment, write `pip freeze` without the
   conda-installed packages, download each file with `pip download --no-deps` and record its
   sha256 in `pip-lock.txt`.
3. Modules: `oc module install` the wanted set, then `oc module freeze > modules.freeze.json`; list
   every new module's data licence above, and check that it needs neither the `hg19` annotator nor
   a chain file.
4. Bump the image version in the module, the build-containers matrix and the labels. Published
   version tags are never overwritten.
