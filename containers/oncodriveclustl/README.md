# noncodingdrivers-oncodriveclustl image

Runtime of the `ONCODRIVECLUSTL` module: OncodriveCLUSTL 1.2.0 (PyPI) on python 3.12 (conda-forge),
plus procps for Nextflow's task metrics. OncodriveCLUSTL has no Bioconda or conda-forge package
(the bbglab channel stops at 1.1.3), so the pip part is locked with hashes.

| File                      | Role                                                                                                                                                                                        |
| ------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `environment.yml`         | conda part, requested packages: python 3.12, pip, setuptools, wheel, pycurl, procps-ng                                                                                                      |
| `conda-lock.linux-64.txt` | explicit lock resolved from `environment.yml` (every package with URL and MD5)                                                                                                              |
| `requirements.in`         | pip part, requested package: `oncodriveclustl==1.2.0` (its 13 dependencies are exact pins upstream)                                                                                         |
| `requirements.lock`       | every pip package with sha256 hashes, installed with `pip install --no-deps --no-build-isolation --require-hashes`                                                                          |
| `AGPL-3.0.txt`            | full text of the GNU AGPL v3 (the upstream LICENSE holds only the notice), installed in `/usr/share/licenses/oncodriveclustl/`                                                              |
| `NOTICE`                  | the one modification of the installed package and the corresponding source, installed next to the licence                                                                                   |
| `Dockerfile`              | Docker/OCI image, published by `.github/workflows/build-containers.yml` as `ghcr.io/khalofx/noncodingdrivers-oncodriveclustl:1.2.0` after `containers/smoke-test.sh oncodriveclustl` passes |
| `oncodriveclustl.def`     | Apptainer recipe from the same locks, for building the SIF directly on an HPC host                                                                                                          |

## Modification of the installed package

`oncodriveclustl/data/CGC_all_Oct15_10_29_09_2018.tsv`, a COSMIC Cancer Gene Census list bundled
in the PyPI wheel, is replaced by its header line, so the image redistributes no COSMIC content.
OncodriveCLUSTL reads the list only to fill its `CGC` column, which it writes as `Non Available`
for every genome other than hg19, so results on hg38 are unchanged. AGPL-3.0 section 5(a) notice:
`NOTICE`.

## Build the SIF on Karakoram (or any host with fakeroot)

```bash
cd containers/oncodriveclustl
apptainer build --fakeroot noncodingdrivers-oncodriveclustl-1.2.0.sif oncodriveclustl.def
apptainer test noncodingdrivers-oncodriveclustl-1.2.0.sif
```

## Changing a pin

1. Edit `environment.yml` or `requirements.in`.
2. Re-resolve both locks on linux-64 and commit them:
   ```bash
   apptainer exec -B $PWD docker://mambaorg/micromamba:2.0.5 \
       micromamba create -y -p $PWD/env -f environment.yml
   apptainer exec -B $PWD docker://mambaorg/micromamba:2.0.5 \
       micromamba env export -p $PWD/env --explicit --md5 > conda-lock.linux-64.txt
   apptainer exec -B $PWD docker://mambaorg/micromamba:2.0.5 bash -c \
       "env/bin/python -m venv tools && tools/bin/pip install pip-tools==7.4.1 && \
        tools/bin/pip-compile --generate-hashes --strip-extras --unsafe-package pycurl \
            --unsafe-package setuptools --unsafe-package wheel --unsafe-package pip \
            --output-file requirements.lock requirements.in"
   ```
3. Bump the image version in the module, the build-containers matrix and the labels.
   Published version tags are never overwritten.

## Licences

OncodriveCLUSTL AGPL-3.0-or-later (IRB Barcelona; no non-commercial clause, unlike OncodriveFML).
The pipeline runs it as a separate program, so the pipeline code stays MIT. The recipe files in this
folder are MIT, like the pipeline.
