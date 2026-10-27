//
// INPUT_HARMONISATION: vcf, sarek and maf samplesheet rows -> per-sample PASS VCFs (GRCh38,
// normalised, primary chromosomes) and canonical mutation tables (sample chrom pos ref alt
// vtype). Every step is a packaged tool (bcftools, maf2vcf.pl, vcf_assembly_checker, Picard
// LiftoverVcf); the gawk programs in assets/awk only convert formats.
//
//   vcf, sarek : BCFTOOLS_ANNOTATE (UCSC contig names) -> ASSEMBLYCHECK
//   maf        : MAF_CHECK_VCF (SNV rows as sites VCF) -> ASSEMBLYCHECK -> MAF_SPLIT (one MAF per
//                tumour) -> MAF2VCF -> BCFTOOLS_SORT
//   GRCh37/hg19 rows (--allow_hg19_liftover): PICARD_LIFTOVERVCF to GRCh38; rejects are published
//   all       : BCFTOOLS_NORM (split multiallelics, left-align, drop REF mismatches)
//               -> BCFTOOLS_VIEW (PASS or missing FILTER, chr1-22/X/Y, tumour column)
//               -> BCFTOOLS_QUERY -> QUERY_TO_CANONICAL (VCF alleles -> MAF conventions)
//
// The assembly check runs once per samplesheet row, against the FASTA of the declared build,
// before anything else reads the row; a row below the match threshold stops the run with
// ASSEMBLY_MISMATCH (decision D28).
//
include { BCFTOOLS_ANNOTATE as RENAME_CHROMS } from '../../../modules/nf-core/bcftools/annotate/main'
include { GAWK as MAF_CHECK_VCF              } from '../../../modules/nf-core/gawk/main'
include { ASSEMBLYCHECK                      } from '../../../modules/local/assemblycheck'
include { GAWK as MAF_SPLIT                  } from '../../../modules/nf-core/gawk/main'
include { MAF2VCF                            } from '../../../modules/local/maf2vcf'
include { BCFTOOLS_SORT as MAF_VCF_SORT      } from '../../../modules/nf-core/bcftools/sort/main'
include { PICARD_LIFTOVERVCF                 } from '../../../modules/nf-core/picard/liftovervcf/main'
include { BCFTOOLS_NORM                      } from '../../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_VIEW as PASS_PRIMARY      } from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_QUERY                     } from '../../../modules/nf-core/bcftools/query/main'
include { GAWK as QUERY_TO_CANONICAL         } from '../../../modules/nf-core/gawk/main'

workflow INPUT_HARMONISATION {

    take:
    ch_rows        // channel: [ val(meta), path(vcf), path(index), path(maf) ]  vcf, sarek and maf rows
    ch_fasta       // channel: [ val(meta2), path(fasta), path(fai) ]            value; GRCh38 analysis FASTA
    ch_dict        // channel: [ val(meta3), path(dict) ]                        value; its sequence dictionary (liftover)
    ch_source      // channel: [ val(meta4), path(fasta), path(fai) ]            value; hg19 FASTA, or [ [:], [], [] ]
    ch_chain       // channel: [ val(meta5), path(chain) ]                       value; hg19 -> GRCh38 chain, or [ [:], [] ]
    rename_map     // path: bcftools --rename-chrs map (RefSeq / Ensembl names -> UCSC)
    min_match      // val: minimum fraction of matching REF alleles (assembly check)
    min_records    // val: fewer checked SNVs -> low_evidence

    main:

    def maf_split_awk = file("${projectDir}/assets/awk/maf_split.awk", checkIfExists: true)
    def canonical_awk = file("${projectDir}/assets/awk/query_to_canonical.awk", checkIfExists: true)

    def ch_vcf_rows = ch_rows
        .filter { meta, _vcf, _index, _maf -> meta.entry in [ 'vcf', 'sarek' ] }
        .map { meta, vcf, index, _maf -> [ meta + [ sample: meta.id ], vcf, index ?: [] ] }
    def ch_maf_rows = ch_rows
        .filter { meta, _vcf, _index, _maf -> meta.entry == 'maf' }
        .map { meta, _vcf, _index, maf -> [ meta, maf ] }

    //
    // vcf and sarek rows: UCSC contig names, then the build check on the renamed VCF
    //
    RENAME_CHROMS(ch_vcf_rows.map { meta, vcf, index -> [ meta, vcf, index, [], [], [], [], rename_map ] })
    def ch_renamed = RENAME_CHROMS.out.vcf.join(RENAME_CHROMS.out.index, failOnMismatch: true)

    //
    // maf rows: the build check reads the MAF SNV rows themselves (maf2vcf.pl would drop rows
    // that disagree with its FASTA and so hide a wrong build)
    //
    MAF_CHECK_VCF(ch_maf_rows.map { meta, maf -> [ meta, maf, 'check.vcf' ] }, maf_split_awk, false)

    // Each item is paired with the FASTA of its declared build (hg19 rows: the liftover source FASTA)
    def ch_check = ch_renamed.map { meta, vcf, _tbi -> [ meta, vcf ] }
        .mix(MAF_CHECK_VCF.out.output)
        .combine(ch_fasta)
        .combine(ch_source)
        .multiMap { meta, vcf, m2, fasta, fai, m4, sfasta, sfai ->
            data: [ meta, vcf ]
            fasta: isLifted(meta) ? [ m4, sfasta, sfai ] : [ m2, fasta, fai ]
        }
    ASSEMBLYCHECK(ch_check.data, ch_check.fasta, min_match, min_records)
    def ch_checked = ASSEMBLYCHECK.out.summary

    //
    // Only rows that passed the check go on. maf rows: one MAF per tumour, then VCF
    //
    def ch_maf_ok = ch_maf_rows.join(ch_checked, failOnMismatch: false).map { meta, maf, _summary -> [ meta, maf ] }
    MAF_SPLIT(ch_maf_ok.map { meta, maf -> [ meta, maf, 'sample.maf' ] }, maf_split_awk, true)
    def ch_tumour_mafs = MAF_SPLIT.out.output
        .flatMap { meta, mafs ->
            (mafs instanceof List ? mafs : [ mafs ]).collect { maf ->
                def sample = maf.name - ~/\.sample\.maf$/
                [ meta + [ id: sample, sample: sample, row: meta.id, vcf_tumour_id: null ], maf ]
            }
        }
    def ch_m2v = ch_tumour_mafs
        .combine(ch_fasta)
        .combine(ch_source)
        .multiMap { meta, maf, m2, fasta, fai, m4, sfasta, sfai ->
            data: [ meta, maf ]
            fasta: isLifted(meta) ? [ m4, sfasta, sfai ] : [ m2, fasta, fai ]
        }
    MAF2VCF(ch_m2v.data, ch_m2v.fasta)
    MAF_VCF_SORT(MAF2VCF.out.vcf)

    def ch_vcf_ok = ch_renamed.join(ch_checked, failOnMismatch: false).map { meta, vcf, _tbi, _summary -> [ meta, vcf ] }
    def ch_declared = ch_vcf_ok.mix(MAF_VCF_SORT.out.vcf)
        .branch { meta, _vcf ->
            liftover: isLifted(meta)
            grch38: true
        }

    //
    // GRCh37/hg19 rows -> GRCh38 (opt-in, --allow_hg19_liftover); rejected records are published
    //
    PICARD_LIFTOVERVCF(
        ch_declared.liftover,
        ch_dict,
        ch_fasta.map { meta2, fasta, _fai -> [ meta2, fasta ] },
        ch_chain
    )

    //
    // Common path on GRCh38
    //
    def ch_grch38 = ch_declared.grch38
        .mix(PICARD_LIFTOVERVCF.out.vcf_lifted)
        .map { meta, vcf -> [ meta + [ assembly: 'GRCh38' ], vcf, [] ] }
    BCFTOOLS_NORM(ch_grch38, ch_fasta.map { meta2, fasta, _fai -> [ meta2, fasta ] })
    PASS_PRIMARY(BCFTOOLS_NORM.out.vcf.map { meta, vcf -> [ meta, vcf, [] ] }, [], [], [])
    def ch_pass = PASS_PRIMARY.out.vcf.join(PASS_PRIMARY.out.index, failOnMismatch: true)
    BCFTOOLS_QUERY(ch_pass, [], [], [])
    QUERY_TO_CANONICAL(BCFTOOLS_QUERY.out.output.map { meta, txt -> [ meta, txt, 'canonical.tsv' ] }, canonical_awk, false)

    emit:
    mutations      = QUERY_TO_CANONICAL.out.output        // channel: [ val(meta), path(tsv) ]       one canonical table per sample
    vcf            = ch_pass                              // channel: [ val(meta), path(vcf), path(tbi) ] per-sample PASS VCFs (GRCh38)
    assembly_check = ch_checked                           // channel: [ val(meta), path(tsv) ]       one summary per samplesheet row
    unlifted       = PICARD_LIFTOVERVCF.out.vcf_unlifted  // channel: [ val(meta), path(vcf) ]       liftover rejects
    skipped        = MAF2VCF.out.skipped                  // channel: [ val(meta), path(tsv) ]       MAF rows maf2vcf.pl dropped
}

//
// GRCh37/hg19 rows are lifted to GRCh38; their build check and MAF conversion use the hg19 FASTA
//
def isLifted(meta) {
    return meta.assembly_declared in [ 'GRCh37', 'hg19' ]
}
