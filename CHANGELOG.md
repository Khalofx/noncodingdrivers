# khalofx/noncodingdrivers: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0 - [unreleased]

Initial release of khalofx/noncodingdrivers, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Entry tiers `vcf`, `sarek` and `maf` (`assets/schema_input.json`) through the
  `INPUT_HARMONISATION` subworkflow: contig renaming (nf-core `bcftools/annotate`), MAF split per
  tumour and conversion with maf2vcf.pl (local `MAF2VCF`, vcf2maf 1.6.22), sorting, normalisation,
  PASS and primary-chromosome filtering and the canonical cohort table (nf-core `bcftools/sort`,
  `bcftools/norm`, `bcftools/view`, `bcftools/query`, gawk); per-sample PASS VCFs published. The
  canonical `mutations` tier is kept for the test profile and the regression gates only.
- Assembly check of every row against the FASTA of its declared build (local `ASSEMBLYCHECK`, EBI
  vcf_assembly_checker 0.10.2): below `--assembly_min_match` (0.99) the run stops with
  `ASSEMBLY_MISMATCH` before any caller.
- Opt-in GRCh37/hg19 liftover (`--allow_hg19_liftover`; nf-core `picard/liftovervcf`) with the
  rejected records published; hg19 rows without the opt-in stop with `ASSEMBLY_NOT_ENABLED`.
- Reference datasheet `assets/references/GRCh38.yml` (schema `assets/schema_references.json`): the
  validated GRCh38 reference set with checksums, paths relative to `--reference_base`; explicit
  parameters take precedence. Asset recipes in `docs/references.md`.
- Driver discovery (`DRIVER_DISCOVERY`) on a CDS-excluded Ensembl Regulatory Build catalogue:
  ActiveDriverWGS 1.2.1 (sharded, FDR over all shards) and OncodriveFML 2.5.0 with CADD v1.7
  (licence-gated: `--accept_oncodrivefml_licence`); MutEnricher 1.3.3 as a QC-only caller
  (`--run_mutenricher`). Element coordinates as fragment spans with `n_fragments`.
- P-value combination (`COMBINEPVALUES`): the empirical Brown's method (default) or ACAT sets
  `consensus_q` (Benjamini–Hochberg, `--fdr_threshold`); ACAT, robust rank aggregation and caller
  votes are published as sensitivity columns; `--skip_acat` where the ACAT package is unavailable.
- OncodriveCLUSTL 1.2.0 (`--run_oncodriveclustl`, off by default; local `ONCODRIVECLUSTL`, public
  image `ghcr.io/khalofx/noncodingdrivers-oncodriveclustl:1.2.0` from `containers/oncodriveclustl`
  with the bundled COSMIC CGC list reduced to its header), reported alongside the consensus:
  `p_oncodriveclustl`/`q_oncodriveclustl` (`P_`/`Q_ANALYTICAL`) carried in the consensus table,
  cluster counts in the evidence table. The regions are prefiltered to the ids with at least
  `--oncodriveclustl_element_mutations` SNVs inside the task (`bin/oncodriveclustl_prefilter.py`; exact) and
  run in shards of `--oncodriveclustl_shard_size` ids with one cohort signature and q-values over
  all shards (`bin/merge_oncodriveclustl_shards.py`). It reads the OncodriveFML bgdata genome.
- `--consensus_callers` (default `activedriverwgs,oncodrivefml`): the callers whose p-values are
  combined (`combine_pvalues.R --callers`). Every other enabled caller is reported alongside: its
  `p_`/`q_` columns are carried in the consensus and evidence tables, not combined. Validation
  stops unknown names, MutEnricher, duplicates, a named non-default caller that is not run and an
  empty combined set; combining a non-default caller prints the opt-in warning. Caller roles in
  `driver_evidence.tsv` (`caller_role_<caller>`), the MultiQC run summary and the methods text,
  written only when a caller is reported alongside, so the default outputs are unchanged.
- Variant flags (`VARIANT_FLAGS`, on by default; flag, never filter): ENCODE Blacklist v2, GIAB v3.6
  segmental duplications and low mappability, Umap k100 uniqueness (nf-core `bedtools/intersect`,
  `bedtools/coverage`, `bedtools/subtract`) and gnomAD v4.1 genomes germline flags (nf-core
  `echtvar/anno`; archive built with `echtvar/encode`); element-level fractions and callable
  fraction.
- PCAWG-style post-filter flags (`POSTFILTERFLAGS`: minimum support, mappability, maximum sample
  fraction; `--postfilter_mode pcawg` adds `q_postfilter`) and the long-form evidence table
  (`WRITEDRIVEREVIDENCE`, `evidence/<cohort>.driver_evidence.tsv`, layout
  `assets/schema_driver_evidence.json`). Evidence tiers are assigned downstream from this table.
- Null calibration (`NULL_CALIBRATION`, `--run_null_calibration`, off by default): neutral copies
  of each cohort (each SNV moved within +/- 50 kb to the same trinucleotide, same patient) re-run
  through the callers; inflation factor, null discoveries and QQ points per caller and
  combination, with a one-sided `--null_lambda_band` gate. Reported only.
- Annotation layers of the significant elements (`EXPLAIN`; none filters or re-ranks): candidate
  mutations per patient on the element fragments (nf-core `bedtools/intersect`), nearest GENCODE
  gene TSS with a strand-aware signed distance (nf-core `bedtools/closest`), ENCODE regulatory
  context per tissue and layer (`--run_regulatory_context`; nf-core `bedtools/coverage`,
  `bedtools/map`), Ensembl VEP 116.3 on candidate mutations (`--run_vep`; nf-core
  `ensemblvep/vep`, offline cache) and motifbreakR 2.24.0 on candidate SNVs (`--run_motifbreakr`;
  JASPAR 2026 CORE vertebrates bundled).
- AlphaGenome Variant Impact scores of the candidate SNVs (`--run_alphagenome_avi`, off by
  default, needs `--run_vep`): the Ensembl VEP AVI plugin of the ensembl-vep 116.3 image adds
  `AVI_PHRED` and `AVI_RAW` to the VEP table; gawk adapters write `candidates/avi/<cohort>.avi.tsv`
  (one row per candidate variant with `avi_status`) and `<cohort>.avi.element_summary.tsv`. The
  user-downloaded score file (`--alphagenome_avi`, datasheet key `alphagenome_avi`, checksums
  pending the download) is checked with its `.tbi` before any task (`AVI_NEEDS_VEP`,
  `AVI_FILE_MISSING`, `AVI_INDEX_MISSING`), and `AVI_TABLE` fails when the plugin did not load.
  Annotation only; with the layer off the VEP command is unchanged.
- OpenCRAVAT 3.2.4 on the candidate mutations (`--run_opencravat`, off by default; local
  `OPENCRAVAT` with gawk adapters): `oc run -l hg38` with the annotators `ccre_screen`, `repeat`
  and `clinvar`, joined back to the candidates by chrom, pos, ref and alt
  (`opencravat/<cohort>.opencravat.tsv`). Telemetry is neutralised: every remote address of the
  OpenCRAVAT configuration is a closed local port, every connection attempt is audited
  (`<cohort>.opencravat.network.tsv`) and a non-local one fails the task, and the locally saved
  metrics payload is deleted (`<cohort>.opencravat.telemetry.tsv`). No liftover: the image holds
  no `hg19` annotator and no UCSC chain file, so `encode_tfbs` is not included.
- Optional cohort-level signature QC (`--run_sigprofiler`; nf-core `sigprofiler`).
- MultiQC 1.35 report with a per-cohort run summary (`RUNSUMMARYMQC`).
- Pipeline-owned images `ghcr.io/khalofx/noncodingdrivers-{utils,activedriverwgs,mutenricher,combinepvalues,motifbreakr,opencravat}`,
  built from `containers/` by `.github/workflows/build-containers.yml` and pinned by digest; the
  OncodriveFML image is built locally from `containers/oncodrivefml/`.
- Profiles `test` (synthetic TERT-window cohort), `test_full` (SYNTH-NC, synthetic genome-wide
  cohort with planted drivers) on `Khalofx/noncodingdrivers-test-datasets`, and the `karakoram`
  site profile.
- Regression baseline records (`tests/regression/`) and gate tools (`benchmarks/`); pytest and ruff
  CI.
- Context layers (`CONTEXT_LAYERS`, after `EXPLAIN`; each off by default, annotation only, own
  tables): COSMIC known-cancer context (`--run_cosmic_context`, local `COSMICCONTEXT`,
  `bin/cosmic_context.py` ported from `benchmarks/cosmic_compare.py`): exact GRCh38 matches of the
  candidate mutations in the user's own COSMIC Non-Coding Variants file, given at run time
  (`--cosmic_ncv`, `--cosmic_release`, `--cosmic_ncv_md5`; never in the datasheet), with file-name
  release and assembly checks, MD5 and provenance JSON and exact sample-name self-overlap flags;
  counts in `cosmic/`, COSMIC records only in `cosmic/restricted/` with a NOTICE.
- Cancer gene-list labels (`--cancer_gene_list`, local `CANCERGENECONTEXT`): nearest gene, nearest
  protein-coding gene (nf-core `bedtools/closest` on protein_coding TSSs) and ENCODE-rE2G targets
  against the COSMIC Cancer Gene Census or a plain list (`cancer_genes/`).
- ENCODE-rE2G enhancer-gene links (`--run_enhancer_gene_links`, `--enhancer_gene_tissues`; nf-core
  `bedtools/intersect` and gawk): links per element and file with the self-promoter split and the
  agreement with the nearest gene (`enhancer_gene/`). Coverage table
  `assets/references/encode_re2g.GRCh38.tsv` (12 tissues supported; 407 DNaseOnly files selected:
  the 11-file liver set and the normal tissue samples of the other tissues); an unsupported or
  unstaged tissue stops the run, never a substitute.
- SCREEN Registry V4 cCRE overlay (`--run_ccre_overlay`, optional `--ccre_celltype_tissues`): V4
  classes, covered fraction and dominant class per element (`ccre/`). Datasheet keys
  `enhancer_gene_dir`, `ccre_v4_bed` (ENCFF420VPZ, MD5 checked at start), `ccre_v4_md5`,
  `ccre_v4_rows`, `ccre_v4_source`, `ccre_celltype_dir`, `encode_licence`, `cosmic_licence`;
  staging recipe in `docs/references.md`.

### `Dependencies`

| Tool                  | Version |
| --------------------- | ------- |
| ActiveDriverWGS       | 1.2.1   |
| OncodriveFML          | 2.5.0   |
| MutEnricher           | 1.3.3   |
| EmpiricalBrownsMethod | 1.34.0  |
| ACAT                  | 0.91    |
| RobustRankAggreg      | 1.2.1   |
| bcftools              | 1.25    |
| bedtools              | 2.31.1  |
| csvtk                 | 0.37.0  |
| echtvar               | 0.2.4   |
| gawk                  | 5.3.1   |
| Picard                | 3.5.0   |
| vcf2maf               | 1.6.22  |
| vcf-validator         | 0.10.2  |
| Ensembl VEP           | 116.3   |
| motifbreakR           | 2.24.0  |
| universalmotif        | 1.28.0  |
| VEP AVI plugin        | 116     |
| OpenCRAVAT            | 3.2.4   |
| SigProfilerExtractor  | 1.2.2   |
| MultiQC               | 1.35    |
