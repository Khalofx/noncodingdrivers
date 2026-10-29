# khalofx/noncodingdrivers: References

## The reference datasheet

The pipeline analyses GRCh38. Its reference set is described by one datasheet,
[`assets/references/GRCh38.yml`](../assets/references/GRCh38.yml), validated against
[`assets/schema_references.json`](../assets/schema_references.json) when the run starts. The
datasheet records every file with its checksum and source; the same checksums are kept in
`tests/regression/PROVENANCE.md`, and `tests/test_references_datasheet.py` checks that the two agree.

- Paths that start with `${reference_base}` are completed with `--reference_base <dir>`. Without
  it those entries stay unset, and the matching parameters must be given explicitly.
- Paths that start with `assets/` ship with the pipeline (the JASPAR motifs and the contig
  rename map).
- An explicit parameter always wins over the datasheet value: `--elements_bed`, `--fasta`
  (with `--fasta_fai`, `--fasta_dict`), `--gencode_gtf`, `--vep_cache`, `--oncodrivefml_genome`,
  `--regulatory_peaks_dir`, `--jaspar_motifs`, `--rename_chrs`, `--liftover_chain`,
  `--liftover_source_fasta` (with `--liftover_source_fai`). When only `--fasta` is given, the
  index and dictionary are expected beside it (`<fasta>.fai`, `<fasta without extension>.dict`).
- `--references <file.yml>` replaces the datasheet; it must validate against the same schema.

The CADD score track for OncodriveFML (`--oncodrivefml_score`) is not in the datasheet: it is free
for non-commercial use only and is never bundled or downloaded by the pipeline.

| Datasheet key                                        | File under `reference_base`                                                   | Checksum (MD5)                                                         | Used by                                |
| ---------------------------------------------------- | ----------------------------------------------------------------------------- | ---------------------------------------------------------------------- | -------------------------------------- |
| `fasta`, `fasta_fai`                                 | `hg38.chr.fa`, `hg38.chr.fa.fai` (25 contigs, UCSC names)                     | `838cd9f56ca7d7dd345cb025d188c991`, `963c862ff26b0d914956faaf7192e54b` | vcf, sarek and maf rows                |
| `fasta_dict`                                         | `derived/hg38.chr.dict`                                                       | `eee72059c0f8b6407202baa38ba3be8a`                                     | Picard LiftoverVcf                     |
| `element_catalogue`                                  | `GRCh38_regulatory_build.ccre.noncoding.v50.bed` (595,725 rows, 544,737 ids)  | `39f838e6c81fdb0b152504f2d51939b5`                                     | all callers                            |
| `gtf`                                                | `gencode/gencode.v50.basic.annotation.gtf.gz`                                 | `11e77cf1d4b78b43f57d608e3f7f2202`                                     | nearest-gene table                     |
| `vep_cache`                                          | `vep_cache_116` (archive `homo_sapiens_vep_116_GRCh38.tar.gz`)                | archive `6e10959391069a9eebb8d8d165fd29be`                             | `--run_vep`                            |
| `oncodrivefml_genome`                                | `oncodrivefml_genome/datasets/genomereference/hg38-20161209`                  | official bgdata dataset                                                | OncodriveFML                           |
| `regulatory_peaks_dir`                               | `encode` (49 peak files, 12 tissues; per-file MD5 in `encode/PROVENANCE.tsv`) |                                                                        | `--run_regulatory_context`             |
| `liftover_chain_from_hg19`                           | `liftover/hg19ToHg38.over.chain.gz`                                           | `35887f73fe5e2231656504d1f6430900`                                     | GRCh37/hg19 rows                       |
| `liftover_source_fasta`, `liftover_source_fasta_fai` | `hg19/hg19.fa`, `hg19/hg19.fa.fai` (UCSC, 93 contigs)                         | `530d89d3ef07fdb2a9b3c701fb4ca486`, `e700bec606168273af0920bd1dd49eb9` | GRCh37/hg19 rows                       |
| `blacklist_bed`                                      | `masks/hg38-blacklist.v2.bed.gz` (ENCODE Blacklist v2, 636 regions)           | `83fe6bf8187a64dee8079b80f75ba289`                                     | variant flags, callable territory      |
| `segdup_bed`                                         | `masks/GRCh38_segdups.bed.gz` (GIAB v3.6, 7,432 regions)                      | `72ec3f7db18ef9318c64b1df271703fe`                                     | variant flags                          |
| `low_mappability_bed`                                | `masks/GRCh38_lowmappabilityall.bed.gz` (GIAB v3.6, 673,815 regions)          | `8a78766bb99abf1d25a5b307d49a0324`                                     | variant flags, mappability post-filter |
| `umap_bed`                                           | `masks/k100.umap.bed.gz` (Umap k100 single-read)                              | `1ef87107f449ccaf7e537de9740692fb`                                     | variant flags, callable territory      |
| `gnomad_echtvar`                                     | `gnomad/gnomad.v4.1.genomes.echtvar.zip` (gnomAD v4.1 genomes, echtvar 0.2.4) | `7bc89a78568a424d7694e2d351617bcd`                                     | germline flag                          |

Two derived tables are not in the datasheet because the pipeline builds them from the files
above in every run: the element coordinates for the anchor join (`ELEMENT_COORDS`, from the
catalogue) and the gene TSS BED for the nearest-gene table (`GENE_TSS`, from the GTF, unless
`--gene_tss_bed` is given).

## Asset recipes

Each asset is built once with packaged tools and pinned by checksum in the datasheet. Never edit
a reference in place: a changed file changes every downstream task hash. The commands are the ones
run on Karakoram (`reference_base` = `/storage/projects/khalifa.juma/noncodingdrivers/data/refs`);
the bedtools commands ran in the image of the nf-core `bedtools` modules (bedtools 2.31.1), the
samtools commands with samtools 1.23.1.

### Element catalogue

Ensembl Regulatory Build release 110 (build 20221007), primary chromosomes, element id from the
GFF `ID` attribute with `:` replaced by `_`, 0-based BED. Coding sequence (GENCODE v50 basic `CDS`
features, merged) is subtracted and remnants shorter than 50 bp are dropped. Rebuilt on
26 October 2026: 547,218 elements before the subtraction; the result is byte-identical to the
pinned catalogue (MD5 `39f838e6c81fdb0b152504f2d51939b5`) under both `LC_ALL=C` and `en_US.UTF-8`.

```bash
REF=/storage/projects/khalifa.juma/noncodingdrivers/data/refs
GFF=homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20221007.gff.gz   # MD5 ee3b05c9020f7da4244d563995a579bb
# ftp.ensembl.org over HTTPS was reset from Karakoram; the EBI copy is the same file
curl -L -o $GFF https://ftp.ebi.ac.uk/ensemblorg/pub/release-110/regulation/homo_sapiens/$GFF
md5sum $GFF

# 1. Regulatory Build GFF -> BED (547,218 elements)
zcat $GFF | awk -F'\t' -v OFS='\t' '!/^#/ && $1 ~ /^([0-9]+|X|Y|MT)$/ {
    id = ""; if (match($9, /ID=[^;]+/)) id = substr($9, RSTART + 3, RLENGTH - 3); gsub(":", "_", id)
    c = ($1 == "MT") ? "chrM" : "chr" $1; print c, $4 - 1, $5, id }' \
  | sort -k1,1 -k2,2n > regbuild110.bed

# 2. GENCODE v50 basic CDS, merged; subtract; keep fragments >= 50 bp
zcat $REF/gencode/gencode.v50.basic.annotation.gtf.gz \
  | awk -F'\t' -v OFS='\t' '$3 == "CDS" && $1 ~ /^chr/ { print $1, $4 - 1, $5 }' \
  | sort -k1,1 -k2,2n | bedtools merge -i - > gencode.v50.cds.merged.bed
bedtools subtract -a regbuild110.bed -b gencode.v50.cds.merged.bed \
  | awk -F'\t' -v OFS='\t' '($3 - $2) >= 50' | sort -k1,1 -k2,2n > catalogue.bed
md5sum catalogue.bed $REF/GRCh38_regulatory_build.ccre.noncoding.v50.bed   # both 39f838e6...
```

### Sequence dictionary (Picard)

```bash
samtools dict -o $REF/derived/hg38.chr.dict $REF/hg38.chr.fa   # MD5 eee72059c0f8b6407202baa38ba3be8a, 25 @SQ lines
```

### hg19 FASTA (GRCh37/hg19 rows only)

The build check and the MAF conversion of GRCh37/hg19 rows read the UCSC hg19 FASTA (chr-named).

```bash
mkdir -p $REF/hg19 && cd $REF/hg19
curl -L -O https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz    # MD5 806c02398f5ac5da8ffd6da2d1d5d1a9
curl -L -O https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/md5sum.txt && grep hg19.fa.gz md5sum.txt
gunzip -k hg19.fa.gz && samtools faidx hg19.fa   # hg19.fa MD5 530d89d3ef07fdb2a9b3c701fb4ca486; .fai e700bec606168273af0920bd1dd49eb9
```

### Liftover chain

UCSC `hg19ToHg38.over.chain.gz` (MD5 `35887f73fe5e2231656504d1f6430900`) from
`https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/`. The chain also maps to alt and
unplaced contigs that the 25-contig FASTA lacks; Picard LiftoverVcf runs with
`--WARN_ON_MISSING_CONTIG true`, so such records are rejected (`NoTarget`) and published with the
other rejects instead of stopping the run.

### OncodriveFML genome

The official bgdata dataset `datasets/genomereference/hg38-20161209`, staged once with its
`.downloaded` marker and copied unchanged into `oncodrivefml_genome/`. The pipeline reads it
offline.

### Ensembl VEP cache

`homo_sapiens_vep_116_GRCh38.tar.gz` (MD5 `6e10959391069a9eebb8d8d165fd29be`, BSD sum
`56036 26996736`) from `https://ftp.ebi.ac.uk/ensemblorg/pub/release-116/variation/indexed_vep_cache/`,
unpacked into `vep_cache_116/`.

### Variant-flag masks

Four BEDs under `masks/`, used unchanged (UCSC names, GRCh38/hg38). Retrieved 27 October 2026.

```bash
cd $REF/masks
GI=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/GRCh38@all
curl -L -O $GI/SegmentalDuplications/GRCh38_segdups.bed.gz          # MD5 72ec3f7db18ef9318c64b1df271703fe
curl -L -O $GI/Mappability/GRCh38_lowmappabilityall.bed.gz          # MD5 8a78766bb99abf1d25a5b307d49a0324
curl -L -O https://bismap.hoffmanlab.org/raw/hg38/k100.umap.bed.gz  # MD5 1ef87107f449ccaf7e537de9740692fb (track line first)
# ENCODE Blacklist v2: lists/hg38-blacklist.v2.bed.gz of Boyle-Lab/Blacklist at commit 61a04d2c5e49341d76735d485c61f0d1177d08a8
curl -L -O https://raw.githubusercontent.com/Boyle-Lab/Blacklist/61a04d2c5e49341d76735d485c61f0d1177d08a8/lists/hg38-blacklist.v2.bed.gz   # MD5 83fe6bf8187a64dee8079b80f75ba289
```

The callable territory of the element-level flags (Umap k100 unique sequence minus Blacklist v2)
is derived in every run (`CALLABLE_TERRITORY`, nf-core `bedtools/subtract`) and published under
`references/`.

### gnomAD v4.1 genomes echtvar archive

Built once from the 24 official per-chromosome sites VCFs of gnomAD v4.1 genomes (chr1-22, X, Y;
460.93 GiB in total), streamed from the gnomAD bucket and checked against the MD5 that Google
Cloud Storage publishes for each object (`x-goog-hash`), so the raw VCFs are never stored. INFO
is reduced to AC, AN, AF, nhomalt and `fafmax_faf95_max` (grpmax FAF95) with bcftools 1.23.1 (the
image of the nf-core `bcftools` modules); FILTER is kept. The archive is then encoded with
echtvar 0.2.4 (the image of nf-core `echtvar/encode`) and the field configuration below, which
stores the frequencies to 5e-7 (multiplier 2,000,000) and writes -1 for a missing value. Built
on 27 October 2026; per-chromosome source MD5s in `gnomad/SOURCE_MD5S.txt`; field configuration MD5 `06f2ebddcf38deebfbd8873ed5e40986`.

```bash
cd $REF/gnomad
cat > gnomad_v4.1_genomes.echtvar.json <<'EOF'
[
  {"field": "FILTER", "alias": "gnomad_filter"},
  {"field": "AC", "alias": "gnomad_ac", "missing_value": -1},
  {"field": "AN", "alias": "gnomad_an", "missing_value": -1},
  {"field": "nhomalt", "alias": "gnomad_nhomalt", "missing_value": -1},
  {"field": "AF", "alias": "gnomad_af", "multiplier": 2000000, "zigzag": true, "missing_value": -1},
  {"field": "fafmax_faf95_max", "alias": "gnomad_fafmax_faf95_max", "multiplier": 2000000, "zigzag": true, "missing_value": -1}
]
EOF
U=https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes
for c in chr{1..22} chrX chrY; do      # run 8 chromosomes at a time; each stream is CPU-bound in bcftools
  curl -sS --fail $U/gnomad.genomes.v4.1.sites.$c.vcf.bgz \
    | tee >(md5sum > stripped/$c.md5) \
    | bcftools annotate --threads 2 -x '^INFO/AC,INFO/AN,INFO/AF,INFO/nhomalt,INFO/fafmax_faf95_max' -Ob -o stripped/$c.bcf -
  # compare stripped/$c.md5 with the base64-decoded md5 of: curl -sI $U/gnomad.genomes.v4.1.sites.$c.vcf.bgz | grep x-goog-hash
done
echtvar encode gnomad.v4.1.genomes.echtvar.zip gnomad_v4.1_genomes.echtvar.json stripped/chr{1..22}.bcf stripped/chrX.bcf stripped/chrY.bcf
md5sum gnomad.v4.1.genomes.echtvar.zip      # 7bc89a78568a424d7694e2d351617bcd, 4,863,222,157 bytes, 759,336,320 variants (1,371 s)
```

echtvar annotates a variant gnomAD does not carry with `gnomad_filter=MISSING` and -1, and joins
several FILTER values with commas (`AC0,AS_VQSR`). gnomAD data are released under CC0 1.0; the
gnomAD terms ask for attribution and bind users not to attempt to reidentify participants.

### Contig rename map

`assets/references/GRCh38.rename_chrs.map` (479 rows, MD5 `91ac4e7b7e10ae19a8ffbe5270642711`) maps
the Ensembl/NCBI primary names (`1`..`22`, `X`, `Y`, `MT`) and the RefSeq accessions of the GRCh38
contigs (the map used for the frozen regression baseline) to the UCSC names of the analysis FASTA. It is
applied to vcf and sarek rows with `bcftools annotate --rename-chrs`; contigs not in the map keep
their name.
