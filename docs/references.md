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
COSMIC files (`--cosmic_ncv`, and a Census given as `--cancer_gene_list`) are not in the datasheet
either: they are licensed, user-supplied at run time and only the current release should be used
([usage](usage.md#cosmic)); the datasheet records only `cosmic_licence`.

The AlphaGenome AVI SNV scores (`alphagenome_avi`) are in the datasheet but are downloaded by the
user (Google sign-in; `docs/usage.md`) and never bundled. Until that download has been inspected,
the file name, the three `alphagenome_avi*_md5` keys and the download date are placeholders: the
MD5 keys hold `pending-user-download`, which the datasheet test accepts for these keys only.

| Datasheet key                                        | File under `reference_base`                                                                                                                    | Checksum (MD5)                                                         | Used by                                |
| ---------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------- | -------------------------------------- |
| `fasta`, `fasta_fai`                                 | `hg38.chr.fa`, `hg38.chr.fa.fai` (25 contigs, UCSC names)                                                                                      | `838cd9f56ca7d7dd345cb025d188c991`, `963c862ff26b0d914956faaf7192e54b` | vcf, sarek and maf rows                |
| `fasta_dict`                                         | `derived/hg38.chr.dict`                                                                                                                        | `eee72059c0f8b6407202baa38ba3be8a`                                     | Picard LiftoverVcf                     |
| `element_catalogue`                                  | `GRCh38_regulatory_build.ccre.noncoding.v50.bed` (595,725 rows, 544,737 ids)                                                                   | `39f838e6c81fdb0b152504f2d51939b5`                                     | all callers                            |
| `gtf`                                                | `gencode/gencode.v50.basic.annotation.gtf.gz`                                                                                                  | `11e77cf1d4b78b43f57d608e3f7f2202`                                     | nearest-gene table                     |
| `vep_cache`                                          | `vep_cache_116` (archive `homo_sapiens_vep_116_GRCh38.tar.gz`)                                                                                 | archive `6e10959391069a9eebb8d8d165fd29be`                             | `--run_vep`                            |
| `oncodrivefml_genome`                                | `oncodrivefml_genome/datasets/genomereference/hg38-20161209`                                                                                   | manifest `afe063fb3a742af2c16083e9a7e05d8a` (below)                    | OncodriveFML, OncodriveCLUSTL          |
| `regulatory_peaks_dir`                               | `encode` (49 peak files, 12 tissues; per-file MD5 in `encode/PROVENANCE.tsv`)                                                                  |                                                                        | `--run_regulatory_context`             |
| `liftover_chain_from_hg19`                           | `liftover/hg19ToHg38.over.chain.gz`                                                                                                            | `35887f73fe5e2231656504d1f6430900`                                     | GRCh37/hg19 rows                       |
| `liftover_source_fasta`, `liftover_source_fasta_fai` | `hg19/hg19.fa`, `hg19/hg19.fa.fai` (UCSC, 93 contigs)                                                                                          | `530d89d3ef07fdb2a9b3c701fb4ca486`, `e700bec606168273af0920bd1dd49eb9` | GRCh37/hg19 rows                       |
| `blacklist_bed`                                      | `masks/hg38-blacklist.v2.bed.gz` (ENCODE Blacklist v2, 636 regions)                                                                            | `83fe6bf8187a64dee8079b80f75ba289`                                     | variant flags, callable territory      |
| `segdup_bed`                                         | `masks/GRCh38_segdups.bed.gz` (GIAB v3.6, 7,432 regions)                                                                                       | `72ec3f7db18ef9318c64b1df271703fe`                                     | variant flags                          |
| `low_mappability_bed`                                | `masks/GRCh38_lowmappabilityall.bed.gz` (GIAB v3.6, 673,815 regions)                                                                           | `8a78766bb99abf1d25a5b307d49a0324`                                     | variant flags, mappability post-filter |
| `umap_bed`                                           | `masks/k100.umap.bed.gz` (Umap k100 single-read)                                                                                               | `1ef87107f449ccaf7e537de9740692fb`                                     | variant flags, callable territory      |
| `gnomad_echtvar`                                     | `gnomad/gnomad.v4.1.genomes.echtvar.zip` (gnomAD v4.1 genomes, echtvar 0.2.4)                                                                  | `7bc89a78568a424d7694e2d351617bcd`                                     | germline flag                          |
| `alphagenome_avi`                                    | `alphagenome/avi_scores_snvs.tsv.gz` (user-downloaded; name pending)                                                                           | pending the download (S7A.15)                                          | `--run_alphagenome_avi`                |
| `enhancer_gene_dir`                                  | `encode_re2g` (407 ENCODE-rE2G files, 12 tissues; per-file MD5 in `encode_re2g/PROVENANCE.tsv` and `assets/references/encode_re2g.GRCh38.tsv`) |                                                                        | `--run_enhancer_gene_links`            |
| `ccre_v4_bed`                                        | `screen_v4/GRCh38-cCREs.V4.ENCFF420VPZ.bed.gz` (SCREEN Registry V4 agnostic cCREs, 2,348,854 rows)                                             | `ea1820bd71581ba01c94fc600fbf794d`                                     | `--run_ccre_overlay`                   |
| `ccre_celltype_dir`                                  | `screen_v4/celltype` (V4 biosample cCREs; MD5 in `screen_v4/PROVENANCE.tsv`)                                                                   |                                                                        | `--ccre_celltype_tissues`              |

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
offline, for OncodriveFML and OncodriveCLUSTL. It is a directory of 27 files (4,427,776,705 bytes;
`oncodrivefml_genome_bytes`), so the datasheet records a manifest checksum
(`oncodrivefml_genome_manifest_md5`), the MD5 of the `md5sum` list of every file in C-locale order:

```bash
cd oncodrivefml_genome/datasets/genomereference/hg38-20161209
find . -type f | LC_ALL=C sort | xargs md5sum | md5sum   # afe063fb3a742af2c16083e9a7e05d8a
```

`package.tar.xz` (1,339,490,304 bytes) is the unextracted archive and is not read; it is part of
the manifest because it is part of the official dataset.

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

### AlphaGenome AVI SNV scores

"AVI SNV scores - AVI scores and Phred-scaled scores in Tabix format" (`avi_scores_snvs_tabix.zip`,
88.5 GB) from <https://deepmind.google.com/science/alphagenome/downloads>, listed there as a
permissive-use download (commercial and non-commercial use); a Google sign-in is needed. Unpacked
into `alphagenome/`; the score file must be bgzip TSV with the columns chrom, pos, ref, alt, raw,
phred (the order `AVI.pm` reads) and a tabix `.tbi` (the plugin does not accept `.csi`). No
versioned release identifier exists for the scores, so provenance is the download date and the
MD5s of the archive, the score file and the index (`alphagenome_avi_archive_md5`,
`alphagenome_avi_md5`, `alphagenome_avi_tbi_md5`). Pending (to-do S7A.15): none of these has been
recorded yet. The plugin itself (`AVI.pm`, MD5 `798aff5678f249fdb2f8b32c76e4733b`, Apache-2.0) is
part of the pinned ensembl-vep 116.3 image.

For tests a synthetic file in the same layout, with invented scores, is written by
`tests/data/avi/make_avi_fixture.sh <candidates.tsv> <out.tsv.gz>` (bgzip and tabix on PATH).

### Contig rename map

`assets/references/GRCh38.rename_chrs.map` (479 rows, MD5 `91ac4e7b7e10ae19a8ffbe5270642711`) maps
the Ensembl/NCBI primary names (`1`..`22`, `X`, `Y`, `MT`) and the RefSeq accessions of the GRCh38
contigs (the map used for the frozen regression baseline) to the UCSC names of the analysis FASTA. It is
applied to vcf and sarek rows with `bcftools annotate --rename-chrs`; contigs not in the map keep
their name.

### ENCODE-rE2G and SCREEN V4

Public ENCODE files, no login. The selection is `assets/references/encode_re2g.GRCh38.tsv`: every
released GRCh38 ENCODE-rE2G "thresholded element gene links" file (DNaseOnly model) of the 12
tissue folders, matched by the ENCODE organ slim, with `selected = yes` for the 11-file liver set
(decision D-E1) and for the normal tissue samples of the other tissues; and
`assets/references/screen_v4_celltype.GRCh38.tsv` for the V4 biosample cCRE files (liver: HepG2,
liver). The table was built on 2 November 2026 from the portal search
`type=File&output_type=thresholded element gene links&assembly=GRCh38&status=released`, keeping
the files whose alias is `encode_re2g_predictions_<biosample>_<ENCSR>_<model>_thresholded_predictions`
(1,458 of 3,078; the rest are ABC, EPIraction and GraphReg predictions). The recipe downloads the
selected files, checks each MD5 and size against the table, checks the rE2G layout (the first ten
columns `#chr start end name class TargetGene TargetGeneEnsemblID TargetGeneTSS isSelfPromoter
CellType` and a last column `Score`) and writes a `PROVENANCE.tsv` per folder with the columns of
`encode/PROVENANCE.tsv` (`n_peaks` holds the number of records). Run on Karakoram on 2 November
2026: 407 rE2G files (1.41 GB) and 3 V4 files, every MD5 matching.

```bash
REF=/storage/projects/khalifa.juma/noncodingdrivers/data/refs
PIPE=<checkout of this repository>
E=https://www.encodeproject.org/files
HEADER='tissue\tlayer\taccession\toutput_type\turl\tbytes\tmd5\tn_peaks\tformat_cols\tnon_chr_rows\n'
RE2G_COLS='#chr start end name class TargetGene TargetGeneEnsemblID TargetGeneTSS isSelfPromoter CellType'
set -euo pipefail

# fetch <url> <file> <md5> <bytes>: download unless present, then check MD5 and size
fetch() {
    [ -s "$2" ] || curl -sS --fail -L --retry 3 -o "$2" "$1"
    local got size
    got=$(md5sum < "$2" | cut -d' ' -f1); size=$(stat -c %s "$2")
    [ "$got" = "$3" ] && [ "$size" = "$4" ] || { echo "CHECKSUM MISMATCH $2" >&2; exit 1; }
}
# counts <file>: data rows, columns of the first data row, rows whose chromosome lacks 'chr'
counts() {
    zcat "$1" | awk -F'\t' '!/^#/ { n++; if (n == 1) nf = NF; if ($1 !~ /^chr/) nc++ } END { printf "%d\t%d\t%d", n, nf, nc }'
}

# 1. ENCODE-rE2G thresholded element-gene links (selected rows of the table)
mkdir -p $REF/encode_re2g && printf "$HEADER" > $REF/encode_re2g/PROVENANCE.tsv
while IFS=$'\t' read -r tissue acc md5 bytes name; do
    mkdir -p $REF/encode_re2g/$tissue; f=$REF/encode_re2g/$tissue/$name
    fetch "$E/$acc/@@download/$acc.bed.gz" "$f" "$md5" "$bytes"
    head=$( { zcat "$f" || true; } | head -1 )
    [ "$(echo "$head" | cut -f1-10 | tr '\t' ' ')" = "$RE2G_COLS" ] && [ "$(echo "$head" | awk -F'\t' '{ print $NF }')" = Score ] \
        || { echo "UNEXPECTED LAYOUT $f" >&2; exit 1; }
    printf '%s\tre2g\t%s\tthresholded element gene links\t%s\t%s\t%s\t%s\n' "$tissue" "$acc" "$E/$acc/@@download/$acc.bed.gz" "$bytes" "$md5" "$(counts "$f")" >> $REF/encode_re2g/PROVENANCE.tsv
done < <(awk -F'\t' -v OFS='\t' 'NR > 1 && $11 == "yes" { print $1, $4, $7, $8, $10 }' $PIPE/assets/references/encode_re2g.GRCh38.tsv)

# 2. SCREEN Registry V4: agnostic cCREs (ENCSR800VNX / ENCFF420VPZ) and the biosample cCRE files
mkdir -p $REF/screen_v4 && printf "$HEADER" > $REF/screen_v4/PROVENANCE.tsv
fetch $E/ENCFF420VPZ/@@download/ENCFF420VPZ.bed.gz $REF/screen_v4/GRCh38-cCREs.V4.ENCFF420VPZ.bed.gz ea1820bd71581ba01c94fc600fbf794d 32795453
printf 'all\tccre_v4\tENCFF420VPZ\tcandidate Cis-Regulatory Elements\t%s\t32795453\tea1820bd71581ba01c94fc600fbf794d\t%s\n' "$E/ENCFF420VPZ/@@download/ENCFF420VPZ.bed.gz" "$(counts $REF/screen_v4/GRCh38-cCREs.V4.ENCFF420VPZ.bed.gz)" >> $REF/screen_v4/PROVENANCE.tsv
while IFS=$'\t' read -r tissue acc md5 bytes name; do
    mkdir -p $REF/screen_v4/celltype/$tissue; f=$REF/screen_v4/celltype/$tissue/$name
    fetch "$E/$acc/@@download/$acc.bed.gz" "$f" "$md5" "$bytes"
    printf '%s\tccre_v4_celltype\t%s\tcandidate Cis-Regulatory Elements\t%s\t%s\t%s\t%s\n' "$tissue" "$acc" "$E/$acc/@@download/$acc.bed.gz" "$bytes" "$md5" "$(counts "$f")" >> $REF/screen_v4/PROVENANCE.tsv
done < <(awk -F'\t' -v OFS='\t' 'NR > 1 { print $1, $3, $5, $6, $7 }' $PIPE/assets/references/screen_v4_celltype.GRCh38.tsv)
```

To add a biosample (for example HuH-7 for liver), set `selected` to `yes` in the table and rerun
the recipe; a file staged outside the table is used with a warning. The pipeline checks the V4
agnostic file against `ccre_v4_md5` when it comes from the datasheet. ENCODE data are free to
download, analyse and publish (ENCODE data use policy); cite the ENCODE Consortium (2020) and the
producing labs.

### COSMIC (user-supplied)

COSMIC files are never staged by the pipeline or pinned in the datasheet (licence; current release
only). Download them with your own account as in [usage](usage.md#cosmic) and keep them outside the
repository, for example `${reference_base}/cosmic/<release>/`. On Karakoram the user's v104
Non-Coding Variants file (tar MD5 `08a3a2970b72a1190c0130c072831fd7`, as published; extracted
`.tsv.gz` MD5 `bde811cfcf6428d2671710c48c7c23dd`) is linked under `data/refs/cosmic/v104/`.
