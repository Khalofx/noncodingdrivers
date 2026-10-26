# noncodingdrivers-motifbreakr image

Runtime of the `MOTIFBREAKR` module: motifbreakR 2.24.0, universalmotif 1.28.0 and
BSgenome.Hsapiens.UCSC.hg38 1.4.5 on R 4.5.3 (Bioconda / conda-forge), plus procps for
Nextflow's task metrics.

| File                      | Role                                                                                                                                                                                 |
| ------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `environment.yml`         | the requested packages (pins)                                                                                                                                                        |
| `conda-lock.linux-64.txt` | explicit lock resolved from `environment.yml`: every package with URL and MD5 (285 packages). Both images install from this file only                                                |
| `Dockerfile`              | Docker/OCI image, published by `.github/workflows/build-containers.yml` as `ghcr.io/khalofx/noncodingdrivers-motifbreakr:2.24.0` after `containers/smoke-test.sh motifbreakr` passes |
| `motifbreakr.def`         | Apptainer recipe from the same lock, for building the SIF directly on an HPC host                                                                                                    |

## Build the SIF on Karakoram (or any host with fakeroot)

```bash
cd containers/motifbreakr
apptainer build --fakeroot noncodingdrivers-motifbreakr-2.24.0.sif motifbreakr.def
apptainer test noncodingdrivers-motifbreakr-2.24.0.sif
```

Nextflow finds a pre-built SIF without pulling when it sits in the Apptainer cache directory
under the name it derives from the module's container,
`ghcr.io-khalofx-noncodingdrivers-motifbreakr-2.24.0.img`.

## Changing a pin

1. Edit `environment.yml` (and its copy `modules/local/motifbreakr/environment.yml`).
2. Re-resolve the lock in a micromamba container and commit it:
   ```bash
   apptainer exec -B $PWD docker://mambaorg/micromamba:2.0.5 \
       micromamba create -y -p $PWD/env -f environment.yml
   apptainer exec -B $PWD docker://mambaorg/micromamba:2.0.5 \
       micromamba env export -p $PWD/env --explicit --md5 > conda-lock.linux-64.txt
   ```
3. Bump the image version in the module, the build-containers matrix and the labels.
   Published version tags are never overwritten.

## Licences

motifbreakR GPL-2, universalmotif GPL-3, the BSgenome packages Artistic-2.0; the UCSC hg38
sequence is freely redistributable. The recipe files in this folder are MIT, like the pipeline.
