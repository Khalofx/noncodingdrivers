//
// Subworkflow with functionality specific to the khalofx/noncodingdrivers pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //

    def before_text = ""
    def after_text = ""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command,
        null
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Licence and parameter gates for the callers (before any task runs)
    //
    def refs = resolveReferences()
    validateDriverParameters(refs)
    validateContextParameters(refs)

    //
    // Create channel from input file provided through params.input:
    // one [ meta, mutations, vcf, index, maf ] tuple per samplesheet row
    //
    def rows = samplesheetToList(input, "${projectDir}/assets/schema_input.json")
    validateInputRows(rows, refs)
    if (params.run_mutenricher && rows.every { row -> row[0].entry == 'mutations' }) {
        log.warn("--run_mutenricher: MutEnricher reads per-sample VCFs (vcf, sarek or maf rows); " +
                 "this samplesheet has only mutations rows, so MutEnricher is not run.")
    }
    if (params.run_mutenricher) {
        rows.groupBy { row -> row[0].cohort ?: row[0].id }.each { cohort, cohort_rows ->
            def entries = cohort_rows.collect { row -> row[0].entry }.unique()
            if ('mutations' in entries && entries.size() > 1) {
                log.warn("--run_mutenricher: cohort '${cohort}' mixes mutations rows with vcf, sarek or maf rows; " +
                         "MutEnricher (QC-only) sees only the donors of the VCF and MAF rows.")
            }
        }
    }
    channel
        .fromList(rows)
        .map { meta, mutations, vcf, index, maf ->
            // An empty cohort makes the row its own cohort. meta.assembly is the analysis assembly
            // (GRCh38); meta.assembly_declared keeps the samplesheet value (GRCh37/hg19 rows are lifted)
            def declared = meta.assembly_declared ?: 'GRCh38'
            [ meta + [ cohort: meta.cohort ?: meta.id, assembly: 'GRCh38', assembly_declared: declared ], mutations, vcf, index, maf ]
        }
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.isBound() ? multiqc_reports.getVal() : [],
            )
        }

        completionSummary(monochrome_logs)

    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Validate samplesheet rows before any task runs: entry tiers, assembly gates and the
// references each tier needs. GRCh38/hg38 rows are analysed as declared (the assembly check
// of INPUT_HARMONISATION confirms it); GRCh37/hg19 rows are lifted to GRCh38 only with
// --allow_hg19_liftover. A mutations table is the canonical GRCh38 tier and is not lifted.
//
def validateInputRows(rows, refs) {
    def hg19 = [ 'GRCh37', 'hg19' ]
    rows.each { row ->
        def meta = row[0]
        def assembly = meta.assembly_declared ?: 'GRCh38'
        if (assembly in hg19) {
            if (meta.entry == 'mutations') {
                error("ASSEMBLY_NOT_SUPPORTED: samplesheet row '${meta.id}': a mutations table is the canonical GRCh38 tier and is not lifted. " +
                      "Supply the calls as a vcf or maf row with assembly ${assembly} (and --allow_hg19_liftover true), or lift the table first.")
            }
            if (!params.allow_hg19_liftover) {
                error("ASSEMBLY_NOT_ENABLED: samplesheet row '${meta.id}' declares assembly ${assembly}. The pipeline analyses GRCh38; " +
                      "set --allow_hg19_liftover true to lift the row with Picard LiftoverVcf (rejected records are published), or supply GRCh38 calls.")
            }
        }
        if (meta.entry == 'sarek' && !meta.vcf_tumour_id) {
            error("VCF_TUMOUR_ID_MISSING: samplesheet row '${meta.id}' (entry sarek): set vcf_tumour_id to the tumour column of the tumour/normal VCF.")
        }
        if (meta.entry == 'vcf' && !meta.vcf_tumour_id) {
            def samples = vcfSampleNames(row[2])
            if (samples.size() > 1) {
                error("VCF_TUMOUR_ID_MISSING: samplesheet row '${meta.id}': the VCF has ${samples.size()} sample columns (${samples.join(', ')}); " +
                      "set vcf_tumour_id to the tumour column.")
            }
        }
    }
    def variant_rows = rows.findAll { row -> row[0].entry != 'mutations' }
    if (variant_rows && !refs.fasta) {
        error("FASTA_MISSING: vcf, sarek and maf rows are harmonised against the GRCh38 FASTA. Set --fasta (UCSC names, .fai beside it) " +
              "or --reference_base for the datasheet assets/references/${params.genome}.yml.")
    }
    def lifted_rows = variant_rows.findAll { row -> (row[0].assembly_declared ?: 'GRCh38') in hg19 }
    if (lifted_rows) {
        def ids = lifted_rows.collect { row -> row[0].id }.join(', ')
        if (!refs.liftover_source_fasta) {
            error("HG19_FASTA_MISSING: rows ${ids} declare GRCh37/hg19. Their build check and MAF conversion read the hg19 FASTA: " +
                  "set --liftover_source_fasta (UCSC hg19, chr-named, .fai beside it) or --reference_base.")
        }
        if (!refs.liftover_chain) {
            error("LIFTOVER_CHAIN_MISSING: rows ${ids} declare GRCh37/hg19: set --liftover_chain (UCSC hg19ToHg38.over.chain.gz) or --reference_base.")
        }
        if (!refs.fasta_dict) {
            error("FASTA_DICT_MISSING: Picard LiftoverVcf needs the sequence dictionary of the GRCh38 FASTA: set --fasta_dict or --reference_base.")
        }
    }
}

//
// Sample column names of a VCF (plain or bgzipped), read from its #CHROM header line
//
def vcfSampleNames(vcf) {
    def path = vcf instanceof Path ? vcf : file(vcf.toString())
    def stream = path.name.endsWith('.gz') ? new java.util.zip.GZIPInputStream(path.newInputStream()) : path.newInputStream()
    def reader = new BufferedReader(new InputStreamReader(stream))
    def header = reader.lines().filter { line -> !line.startsWith('##') }.findFirst().orElse('')
    reader.close()
    return header.startsWith('#CHROM') ? header.tokenize('\t').drop(9) : []
}

//
// Reference set of the run. The datasheet assets/references/<genome>.yml (or --references) is
// validated by samplesheetToList against assets/schema_references.json. `${reference_base}` is
// substituted here (nf-schema 2.5 has no parameter substitution); entries that need it stay unset
// when --reference_base is not given. Relative entries (assets/...) resolve against the pipeline
// directory. An explicit parameter always wins over the datasheet value (pickRef); a companion
// file (.fai, .dict) follows its primary file when only the primary is given explicitly.
//
def resolveReferences() {
    def sheet = [:]
    if (params.genome || params.references) {
        def sheet_path = params.references ?: "${projectDir}/assets/references/${params.genome}.yml"
        def entries = samplesheetToList(sheet_path, "${projectDir}/assets/schema_references.json")
        def first = entries[0]
        sheet = first instanceof Map ? first : first[0]
    }
    def refs = [
        genome            : sheet.genome ?: params.genome,
        datasheet_version : sheet.datasheet_version,
        elements_bed      : pickRef(sheet, 'elements_bed', 'element_catalogue'),
        gencode_gtf       : pickRef(sheet, 'gencode_gtf', 'gtf'),
        jaspar_motifs     : pickRef(sheet, 'jaspar_motifs', 'jaspar_motifs') ?: "${projectDir}/assets/motifs/jaspar_core_vertebrates_2026.pfm",
        vep_cache         : pickRef(sheet, 'vep_cache', 'vep_cache'),
        oncodrivefml_genome  : pickRef(sheet, 'oncodrivefml_genome', 'oncodrivefml_genome'),
        regulatory_peaks_dir : pickRef(sheet, 'regulatory_peaks_dir', 'regulatory_peaks_dir'),
        rename_chrs       : pickRef(sheet, 'rename_chrs', 'rename_map'),
        liftover_chain    : pickRef(sheet, 'liftover_chain', 'liftover_chain_from_hg19'),
        blacklist_bed       : pickRef(sheet, 'blacklist_bed', 'blacklist_bed'),
        segdup_bed          : pickRef(sheet, 'segdup_bed', 'segdup_bed'),
        low_mappability_bed : pickRef(sheet, 'low_mappability_bed', 'low_mappability_bed'),
        umap_bed            : pickRef(sheet, 'umap_bed', 'umap_bed'),
        gnomad_echtvar      : pickRef(sheet, 'gnomad_echtvar', 'gnomad_echtvar'),
        alphagenome_avi     : pickRef(sheet, 'alphagenome_avi', 'alphagenome_avi'),
        enhancer_gene_dir   : pickRef(sheet, 'enhancer_gene_dir', 'enhancer_gene_dir'),
        ccre_v4_bed         : pickRef(sheet, 'ccre_v4_bed', 'ccre_v4_bed'),
        ccre_celltype_dir   : pickRef(sheet, 'ccre_celltype_dir', 'ccre_celltype_dir'),
    ]
    // checksum of the datasheet's V4 file; not applied to an explicit --ccre_v4_bed
    refs.ccre_v4_md5 = params.ccre_v4_bed ? null : sheet.ccre_v4_md5
    // FASTA and its companions
    refs.fasta = pickRef(sheet, 'fasta', 'fasta')
    if (params.fasta) {
        refs.fasta_fai  = params.fasta_fai ?: "${params.fasta}.fai"
        refs.fasta_dict = params.fasta_dict ?: params.fasta.toString().replaceAll(/\.(fa|fasta|fna)$/, '') + '.dict'
    } else {
        refs.fasta_fai  = pickRef(sheet, 'fasta_fai', 'fasta_fai')
        refs.fasta_dict = pickRef(sheet, 'fasta_dict', 'fasta_dict')
    }
    refs.liftover_source_fasta = pickRef(sheet, 'liftover_source_fasta', 'liftover_source_fasta')
    refs.liftover_source_fai = params.liftover_source_fasta
        ? (params.liftover_source_fai ?: "${params.liftover_source_fasta}.fai")
        : pickRef(sheet, 'liftover_source_fai', 'liftover_source_fasta_fai')
    return refs
}

//
// An explicit parameter wins; otherwise the datasheet value with ${reference_base} substituted
//
def pickRef(Map sheet, String param_name, String key) {
    def explicit = params[param_name]
    if (explicit != null && explicit.toString() != '') {
        return explicit.toString()
    }
    def value = sheet[key]
    if (value == null || value.toString() == '') {
        return null
    }
    def text = value.toString()
    def placeholder = '${reference_base}'
    if (text.contains(placeholder)) {
        if (!params.reference_base) {
            return null
        }
        return text.replace(placeholder, params.reference_base.toString().replaceAll(/\/+$/, ''))
    }
    if (text.startsWith('/') || text ==~ /^[A-Za-z][A-Za-z0-9+.-]*:\/\/.*/) {
        return text
    }
    return "${projectDir}/${text}"
}

//
// Licence and parameter gates for the driver callers
//
def validateDriverParameters(refs) {
    if (!refs.elements_bed) {
        error("No element catalogue: set --elements_bed (headerless BED: chrom start end element_id) or --reference_base for the datasheet assets/references/${params.genome}.yml.")
    }
    if (!params.run_activedriverwgs && !params.run_oncodrivefml) {
        error("No driver caller selected: set --run_activedriverwgs true and/or --run_oncodrivefml true.")
    }
    if (params.run_oncodrivefml && !params.accept_oncodrivefml_licence) {
        error("OncodriveFML is licence-gated: it is AGPL-3.0 software whose rights holder, Universitat Pompeu Fabra, " +
              "permits academic, non-commercial and personal use only without a separate licence. " +
              "Set --accept_oncodrivefml_licence true to confirm that your use complies, or --run_oncodrivefml false.")
    }
    if (params.run_oncodrivefml && !params.oncodrivefml_score) {
        error("--run_oncodrivefml needs --oncodrivefml_score: a bgzip-compressed, tabix-indexed functional-impact score track " +
              "(CADD GRCh38 v1.7 whole-genome SNVs, free for non-commercial use) with its .tbi beside it.")
    }
    if (params.run_mutenricher && !params.mutenricher_mappable_regions) {
        error("--run_mutenricher needs --mutenricher_mappable_regions: a GRCh38 callable/mappable regions BED, " +
              "BGZF-compressed and tabix-indexed (.tbi beside it).")
    }
    if (params.combination_method == 'acat' && params.skip_acat) {
        error("--combination_method acat cannot be used with --skip_acat true: ACAT would not be computed.")
    }
    if (params.run_regulatory_context && !refs.regulatory_peaks_dir) {
        error("--run_regulatory_context needs --regulatory_peaks_dir (or --reference_base for the datasheet): one folder per tissue holding <layer>.<accession>.bed.gz peak files.")
    }
    if (params.run_variant_flags) {
        def masks = [ blacklist_bed: refs.blacklist_bed, segdup_bed: refs.segdup_bed, low_mappability_bed: refs.low_mappability_bed, umap_bed: refs.umap_bed ]
        def missing = masks.findAll { _name, path -> !path }.keySet()
        if (missing) {
            def names = missing.collect { name -> "--${name}" }.join(', ')
            error("VARIANT_FLAG_MASKS_MISSING: --run_variant_flags needs the mask BEDs ${names} " +
                  "(or --reference_base for the datasheet assets/references/${params.genome}.yml), or set --run_variant_flags false.")
        }
        if (!params.skip_germline_flags && !refs.gnomad_echtvar) {
            error("GNOMAD_ECHTVAR_MISSING: the germline flag reads gnomAD v4.1 genomes through echtvar: set --gnomad_echtvar " +
                  "(or --reference_base for the datasheet), or --skip_germline_flags true.")
        }
        if (!params.skip_germline_flags && !refs.fasta) {
            error("GERMLINE_FASTA_MISSING: the germline flag anchors indels on the GRCh38 FASTA for the gnomAD lookup: set --fasta " +
                  "(UCSC names, .fai beside it) or --reference_base, or --skip_germline_flags true.")
        }
    }
    if (params.run_null_calibration && !refs.fasta) {
        error("NULL_FASTA_MISSING: --run_null_calibration keeps each SNV's trinucleotide, read from the GRCh38 FASTA: set --fasta " +
              "(UCSC names, .fai beside it) or --reference_base.")
    }
    if (params.run_vep && !refs.vep_cache) {
        error("--run_vep needs --vep_cache (or --reference_base for the datasheet): an offline Ensembl VEP cache directory (homo_sapiens/${params.vep_cache_version}_GRCh38 inside it).")
    }
    // VEP exits 0 when the AVI plugin cannot load (no plugin folder, no .tbi): check here, before any task
    if (params.run_alphagenome_avi) {
        if (!params.run_vep) {
            error("AVI_NEEDS_VEP: --run_alphagenome_avi reads the AlphaGenome AVI scores through the Ensembl VEP AVI plugin: set --run_vep true as well.")
        }
        if (!refs.alphagenome_avi) {
            error("AVI_FILE_MISSING: --run_alphagenome_avi needs --alphagenome_avi (or --reference_base for the datasheet): the AlphaGenome AVI SNV score file \"AVI SNV scores\" (bgzip TSV, tabix-indexed), downloaded by the user (docs/usage.md).")
        }
        if (!file(refs.alphagenome_avi).exists()) {
            error("AVI_FILE_MISSING: the AlphaGenome AVI score file ${refs.alphagenome_avi} does not exist.")
        }
        if (!file("${refs.alphagenome_avi}.tbi").exists()) {
            error("AVI_INDEX_MISSING: ${refs.alphagenome_avi}.tbi not found. The VEP AVI plugin reads only a tabix .tbi index (not .csi) and VEP does not fail without it: index the file with tabix -s 1 -b 2 -e 2.")
        }
    }
    validateConsensusCallers()
}

//
// Callers that are harmonised and joined on element_id when enabled (MutEnricher is QC-only and
// never joined), with their --run_ flag, and the default combined set of --consensus_callers
//
def harmonisedCallerFlags() {
    return [
        activedriverwgs: params.run_activedriverwgs,
        oncodrivefml   : params.run_oncodrivefml,
        oncodriveclustl: params.run_oncodriveclustl
    ]
}

def defaultConsensusCallers() {
    return [ 'activedriverwgs', 'oncodrivefml' ]
}

def consensusCallerNames() {
    return (params.consensus_callers ?: '').toString().tokenize(',').collect { name -> name.trim() }.findAll { name -> name }
}

//
// Role of every enabled harmonised caller: 'combined' (named in --consensus_callers) or
// 'reported_alongside' (p_/q_ columns carried in the consensus table, not combined)
//
def callerRoles() {
    def named = consensusCallerNames()
    return harmonisedCallerFlags()
        .findAll { _caller, enabled -> enabled }
        .collectEntries { caller, _enabled -> [ (caller): caller in named ? 'combined' : 'reported_alongside' ] }
}

//
// --consensus_callers: names checked here (the schema pattern also restricts them). A default
// caller that is not run is dropped silently (OncodriveFML is off by default); a non-default
// caller must be run when it is named, and combining it prints the opt-in warning.
//
def validateConsensusCallers() {
    def flags   = harmonisedCallerFlags()
    def named   = consensusCallerNames()
    def default_callers = defaultConsensusCallers()
    if (!named) {
        error("--consensus_callers is empty: name at least one caller whose p-values are combined (default: ${default_callers.join(',')}).")
    }
    if ('mutenricher' in named) {
        error("--consensus_callers names mutenricher: MutEnricher is QC-only (its native table is published, never harmonised or combined).")
    }
    def unknown = named.findAll { name -> !flags.containsKey(name) }
    if (unknown) {
        error("--consensus_callers names unknown caller(s) ${unknown.join(', ')}: choose from ${flags.keySet().join(', ')}.")
    }
    def duplicates = named.findAll { name -> named.count(name) > 1 }.unique()
    if (duplicates) {
        error("--consensus_callers names ${duplicates.join(', ')} more than once.")
    }
    def not_run = named.findAll { name -> !(name in default_callers) && !flags[name] }
    if (not_run) {
        error("--consensus_callers names ${not_run.join(', ')} but --run_${not_run[0]} is false: enable the caller or remove it from --consensus_callers.")
    }
    def combined = named.findAll { name -> flags[name] }
    if (!combined) {
        error("None of the callers named in --consensus_callers (${named.join(', ')}) is run: enable one of them with its --run_ flag.")
    }
    combined.findAll { name -> !(name in default_callers) }.each { name ->
        log.warn("--consensus_callers combines ${name}, which the default reports alongside. ${name} has not passed calibration on real cohorts: " +
                 "the opt-in gate (null calibration and spike-in recovery on SYNTH-NC and Liver-HCC-public, docs/validation.md) is not met, " +
                 "and combining it changes the EBM covariance set and every consensus q-value.")
    }
}

//
// Meta map of one ENCODE peak file, <dir>/<tissue>/<layer>.<accession>.bed[.gz]. The signal
// column follows the file layout, read from its first data line: narrowPeak (10 or 20 columns)
// keeps signalValue (fold enrichment) in column 7, DNase hotspot files (8 columns) keep
// max_density in column 5. Any other layout stops the run rather than reading an arbitrary column.
//
def peakFileMeta(peak) {
    def tissue    = peak.parent.name
    def tokens    = peak.name.tokenize('.')
    def layer     = tokens[0]
    def accession = tokens.size() > 2 ? tokens[1] : 'NA'
    def stream    = peak.name.endsWith('.gz') ? new java.util.zip.GZIPInputStream(peak.newInputStream()) : peak.newInputStream()
    def first     = stream.withReader { reader ->
        new BufferedReader(reader).lines()
            .filter { line -> line && !(line.startsWith('#') || line.startsWith('track') || line.startsWith('browser')) }
            .findFirst()
            .orElse(null)
    }
    def ncol       = first ? first.tokenize('\t').size() : 0
    def signal_col = ncol in [ 10, 20 ] ? 7 : (ncol == 8 ? 5 : 0)
    if (!signal_col) {
        error("Peak file ${peak}: ${ncol} columns is not a supported layout. Supported: narrowPeak (10 or 20 columns, " +
              "signalValue in column 7) and DNase hotspot (8 columns, max_density in column 5).")
    }
    return [ id: "${tissue}.${layer}.${accession}", tissue: tissue, layer: layer, accession: accession, signal_col: signal_col ]
}

//
// Context layers (CONTEXT_LAYERS): parameter and file checks before any task runs
//
def validateContextParameters(refs) {
    if (params.run_cosmic_context) {
        if (!params.cosmic_ncv) {
            error("COSMIC_NCV_MISSING: --run_cosmic_context needs --cosmic_ncv: your own copy of the COSMIC Non-Coding Variants file, " +
                  "Cosmic_NonCodingVariants_v<release>_GRCh38.tsv.gz (extracted from the download tar, name unchanged). COSMIC data are " +
                  "licensed (registration for academic use, a QIAGEN licence for commercial use) and are never bundled or downloaded by the pipeline.")
        }
        if (!(params.cosmic_release ==~ /^v\d+$/)) {
            error("COSMIC_RELEASE_MISSING: --run_cosmic_context needs --cosmic_release, the COSMIC release of --cosmic_ncv (e.g. v104).")
        }
        def name  = file(params.cosmic_ncv).name
        def match = name =~ /^Cosmic_NonCodingVariants_(v\d+)_(GRCh3[78])\.tsv(\.gz)?$/
        if (!match.matches()) {
            error("COSMIC_FILE_NAME: --cosmic_ncv '${name}' is not COSMIC's file name Cosmic_NonCodingVariants_v<release>_GRCh38.tsv.gz. " +
                  "Give the .tsv.gz extracted from the download tar (not the .tar) under its own name, so that release and assembly can be checked.")
        }
        if (match.group(2) != 'GRCh38') {
            error("COSMIC_ASSEMBLY: --cosmic_ncv '${name}' is the ${match.group(2)} file set. The pipeline analyses GRCh38 and never lifts COSMIC: " +
                  "download the GRCh38 Non-Coding Variants file.")
        }
        if (match.group(1) != params.cosmic_release) {
            error("COSMIC_RELEASE_MISMATCH: --cosmic_ncv '${name}' is release ${match.group(1)}, but --cosmic_release is ${params.cosmic_release}.")
        }
    } else if (params.cosmic_ncv || params.cosmic_release || params.cosmic_ncv_md5) {
        log.warn("--cosmic_ncv/--cosmic_release/--cosmic_ncv_md5 are set but --run_cosmic_context is false: the COSMIC layer is not run.")
    }
    if (params.cancer_gene_list && !params.gene_tss_bed && !refs.gencode_gtf) {
        log.warn("--cancer_gene_list needs the nearest-gene table (--gencode_gtf or --gene_tss_bed): no cancer gene labels are written.")
    }
    if (params.run_enhancer_gene_links) {
        enhancerGeneFiles(refs).warnings.each { message -> log.warn(message) }
    } else if (params.enhancer_gene_tissues) {
        log.warn("--enhancer_gene_tissues is set but --run_enhancer_gene_links is false: the ENCODE-rE2G layer is not run.")
    }
    if (params.run_ccre_overlay) {
        if (!refs.ccre_v4_bed) {
            error("CCRE_V4_MISSING: --run_ccre_overlay needs --ccre_v4_bed (or --reference_base for the datasheet): the SCREEN Registry V4 agnostic cCREs, " +
                  "ENCODE ENCFF420VPZ (GRCh38-cCREs.V4.ENCFF420VPZ.bed.gz; recipe in docs/references.md).")
        }
        def v4 = file(refs.ccre_v4_bed, checkIfExists: true)
        if (refs.ccre_v4_md5) {
            def md5 = fileMd5(v4)
            if (md5 != refs.ccre_v4_md5) {
                error("CCRE_V4_CHECKSUM: ${v4} has MD5 ${md5}, the datasheet expects ${refs.ccre_v4_md5} (ENCFF420VPZ). Re-stage it with the recipe in docs/references.md.")
            }
        }
        if (params.ccre_celltype_tissues) {
            ccreCelltypeFiles(refs)
        }
    } else if (params.ccre_celltype_tissues) {
        log.warn("--ccre_celltype_tissues is set but --run_ccre_overlay is false: the SCREEN V4 layer is not run.")
    }
}

def fileMd5(path) {
    def digest = java.security.MessageDigest.getInstance('MD5')
    path.withInputStream { stream -> stream.eachByte(1 << 20) { buffer, n -> digest.update(buffer, 0, n) } }
    return digest.digest().encodeHex().toString()
}

//
// ENCODE-rE2G files of the requested tissues: <enhancer_gene_dir>/<tissue>/<biosample>.<accession>.bed.gz.
// A tissue is supported when assets/references/encode_re2g.GRCh38.tsv lists ENCODE-rE2G
// predictions for it; an unsupported or unstaged tissue stops the run (never a substitute).
// Returns [ tissues, files: [ [ meta, path ] ], n_files: 'tissue:N,...', warnings ].
//
def enhancerGeneFiles(refs) {
    if (!refs.enhancer_gene_dir) {
        error("ENHANCER_GENE_DIR_MISSING: --run_enhancer_gene_links needs --enhancer_gene_dir (or --reference_base for the datasheet): one folder per tissue " +
              "holding <biosample>.<accession>.bed.gz ENCODE-rE2G thresholded link files (recipe in docs/references.md).")
    }
    def rows      = file("${projectDir}/assets/references/encode_re2g.GRCh38.tsv", checkIfExists: true).splitCsv(header: true, sep: '\t')
    def supported = rows.findAll { row -> row.accession ==~ /^ENCFF\w+$/ }.collect { row -> row.tissue }.unique().sort()
    def dir       = file(refs.enhancer_gene_dir, checkIfExists: true)
    def staged    = dir.listFiles().findAll { entry -> entry.isDirectory() }.collect { entry -> entry.name }.sort()
    def warnings  = []
    def requested = params.enhancer_gene_tissues
        ? params.enhancer_gene_tissues.toString().tokenize(',').collect { tissue -> tissue.trim() }.unique()
        : staged.findAll { tissue -> tissue in supported }
    def unsupported = requested.findAll { tissue -> !(tissue in supported) }
    if (unsupported) {
        error("ENHANCER_GENE_TISSUE_UNSUPPORTED: ENCODE-rE2G has no thresholded enhancer-gene predictions for ${unsupported.join(', ')} " +
              "(assets/references/encode_re2g.GRCh38.tsv). Supported tissues: ${supported.join(', ')}. No other tissue is used as a substitute: " +
              "remove it from --enhancer_gene_tissues or set --run_enhancer_gene_links false.")
    }
    if (!requested) {
        error("ENHANCER_GENE_TISSUE_NOT_STAGED: no supported tissue folder is staged under ${dir} (supported: ${supported.join(', ')}); " +
              "stage the files with the recipe in docs/references.md or set --enhancer_gene_tissues.")
    }
    def files = []
    requested.each { tissue ->
        def folder = dir.resolve(tissue)
        def beds   = folder.exists() ? folder.listFiles().findAll { entry -> entry.name.endsWith('.bed.gz') }.sort { entry -> entry.name } : []
        if (!beds) {
            error("ENHANCER_GENE_TISSUE_NOT_STAGED: tissue '${tissue}' has ENCODE-rE2G predictions but no <biosample>.<accession>.bed.gz file is staged in ${folder}. " +
                  "Stage them with the recipe in docs/references.md.")
        }
        def listed   = rows.findAll { row -> row.tissue == tissue }.collect { row -> row.file_name }
        def selected = rows.findAll { row -> row.tissue == tissue && row.selected == 'yes' }.collect { row -> row.file_name }
        def names    = beds.collect { bed -> bed.name }
        def unknown  = names.findAll { name -> !(name in listed) }
        def missing  = selected.findAll { name -> !(name in names) }
        if (unknown) {
            warnings << "ENCODE-rE2G tissue '${tissue}': ${unknown.size()} staged file(s) not in assets/references/encode_re2g.GRCh38.tsv are used as staged: ${unknown.join(', ')}".toString()
        }
        if (missing) {
            warnings << "ENCODE-rE2G tissue '${tissue}': ${missing.size()} of the ${selected.size()} selected file(s) are not staged: ${missing.join(', ')}".toString()
        }
        beds.each { bed -> files << [ contextFileMeta(bed, tissue), bed ] }
    }
    def n_files = requested.collect { tissue -> "${tissue}:${files.count { item -> item[0].tissue == tissue }}" }.join(',')
    return [ tissues: requested, files: files, n_files: n_files, warnings: warnings ]
}

//
// SCREEN V4 biosample-specific cCRE files of --ccre_celltype_tissues:
// <ccre_celltype_dir>/<tissue>/<biosample>.<accession>.bed.gz
//
def ccreCelltypeFiles(refs) {
    if (!refs.ccre_celltype_dir) {
        error("CCRE_CELLTYPE_DIR_MISSING: --ccre_celltype_tissues needs --ccre_celltype_dir (or --reference_base for the datasheet): one folder per tissue " +
              "holding <biosample>.<accession>.bed.gz SCREEN V4 biosample cCRE files.")
    }
    def dir   = file(refs.ccre_celltype_dir, checkIfExists: true)
    def files = []
    params.ccre_celltype_tissues.toString().tokenize(',').collect { tissue -> tissue.trim() }.unique().each { tissue ->
        def folder = dir.resolve(tissue)
        def beds   = folder.exists() ? folder.listFiles().findAll { entry -> entry.name.endsWith('.bed.gz') }.sort { entry -> entry.name } : []
        if (!beds) {
            error("CCRE_CELLTYPE_TISSUE_NOT_STAGED: no SCREEN V4 biosample cCRE file is staged for '${tissue}' in ${folder} (recipe in docs/references.md). " +
                  "No other tissue is used as a substitute.")
        }
        beds.each { bed -> files << [ contextFileMeta(bed, tissue), bed ] }
    }
    return [ files: files ]
}

//
// Meta map of a staged context file <tissue>/<biosample>.<accession>.bed.gz
//
def contextFileMeta(bed, tissue) {
    def tokens = bed.name.tokenize('.')
    if (tokens.size() != 4 || !(tokens[1] ==~ /^ENCFF\w+$/)) {
        error("Context file ${bed}: expected <biosample>.<ENCODE accession>.bed.gz (biosample without dots), e.g. HepG2.ENCFF960GTR.bed.gz.")
    }
    return [ id: "${tissue}.${tokens[0]}.${tokens[1]}", tissue: tissue, biosample: tokens[0], accession: tokens[1] ]
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // In-text citations of the tools of this run (MultiQC methods description)
    def citation_text = [
            "Tools used in the workflow included:",
            "bcftools (Danecek et al. 2021),",
            "BEDTools (Quinlan and Hall 2010),",
            params.run_activedriverwgs ? "ActiveDriverWGS (Zhu et al. 2020)," : "",
            params.run_oncodrivefml ? "OncodriveFML (Mularoni et al. 2016) with CADD v1.7 scores (Schubach et al. 2024)," : "",
            params.run_oncodriveclustl ? "OncodriveCLUSTL (Arnedo-Pac et al. 2019)," : "",
            params.run_mutenricher ? "MutEnricher (Soltis et al. 2020)," : "",
            "the empirical Brown's method (Poole et al. 2016),",
            params.skip_acat ? "" : "ACAT (Liu et al. 2019),",
            "RobustRankAggreg (Kolde et al. 2012),",
            params.run_variant_flags ? "echtvar (Pedersen et al. 2023) with gnomAD v4.1 (Chen et al. 2024)," : "",
            params.run_vep ? "Ensembl VEP (McLaren et al. 2016)," : "",
            params.run_motifbreakr ? "motifbreakR (Coetzee et al. 2015) with JASPAR (Rauluseviciute et al. 2024)," : "",
            params.run_alphagenome_avi ? "AlphaGenome Variant Impact scores (Avsec et al. 2026; Cheng et al. 2026) through the VEP AVI plugin," : "",
            params.run_opencravat ? "OpenCRAVAT (Pagel et al. 2020)," : "",
            params.run_cosmic_context ? "COSMIC (Sondka et al. 2024)," : "",
            params.cancer_gene_list ? "a cancer gene list (e.g. the COSMIC Cancer Gene Census, Sondka et al. 2018)," : "",
            params.run_enhancer_gene_links ? "ENCODE-rE2G enhancer-gene predictions (Gschwind et al. 2026; ENCODE Project Consortium 2020)," : "",
            params.run_ccre_overlay ? "the SCREEN Registry V4 of candidate cis-regulatory elements (Moore et al. 2026)," : "",
            "MultiQC (Ewels et al. 2016)",
            ".",
            callerRolesText()
        ].join(' ').replaceAll(/ +/, ' ').trim()

    return citation_text
}

def callerRolesText() {
    // Which callers were combined and which were reported alongside. Empty when every enabled
    // caller is a combined default caller, so the default methods text is unchanged.
    def names = [ activedriverwgs: 'ActiveDriverWGS', oncodrivefml: 'OncodriveFML', oncodriveclustl: 'OncodriveCLUSTL' ]
    def roles = callerRoles()
    def combined  = roles.findAll { _caller, role -> role == 'combined' }.keySet() as List
    def alongside = roles.findAll { _caller, role -> role == 'reported_alongside' }.keySet() as List
    if (!alongside && combined.every { caller -> caller in defaultConsensusCallers() }) {
        return ''
    }
    def text = "Element p-values of ${combined.collect { caller -> names[caller] }.join(' and ')} were combined"
    if (alongside) {
        text += "; ${alongside.collect { caller -> names[caller] }.join(' and ')} ${alongside.size() > 1 ? 'were' : 'was'} reported alongside and not combined"
    }
    return text + "."
}

def toolBibliographyText() {
    // Bibliography of the tools of this run (MultiQC methods description)
    def reference_text = ([
            "<li>Danecek, P., et al. (2021). Twelve years of SAMtools and BCFtools. GigaScience, 10(2). doi: <a href=\"https://doi.org/10.1093/gigascience/giab008\">10.1093/gigascience/giab008</a></li>",
            "<li>Quinlan, A. R., & Hall, I. M. (2010). BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics, 26(6). doi: <a href=\"https://doi.org/10.1093/bioinformatics/btq033\">10.1093/bioinformatics/btq033</a></li>",
            params.run_activedriverwgs ? "<li>Zhu, H., et al. (2020). Candidate cancer driver mutations in distal regulatory elements and long-range chromatin interaction networks. Molecular Cell, 77(6). doi: <a href=\"https://doi.org/10.1016/j.molcel.2019.12.027\">10.1016/j.molcel.2019.12.027</a></li>" : "",
            params.run_oncodrivefml ? "<li>Mularoni, L., et al. (2016). OncodriveFML: a general framework to identify coding and non-coding regions with cancer driver mutations. Genome Biology, 17. doi: <a href=\"https://doi.org/10.1186/s13059-016-0994-0\">10.1186/s13059-016-0994-0</a></li>" : "",
            params.run_oncodrivefml ? "<li>Schubach, M., et al. (2024). CADD v1.7: using protein language models, regulatory CNNs and other nucleotide-level scores to improve genome-wide variant predictions. Nucleic Acids Research, 52(D1). doi: <a href=\"https://doi.org/10.1093/nar/gkad989\">10.1093/nar/gkad989</a></li>" : "",
            params.run_oncodriveclustl ? "<li>Arnedo-Pac, C., et al. (2019). OncodriveCLUSTL: a sequence-based clustering method to identify cancer drivers. Bioinformatics, 35(22), 4788–4790. doi: <a href=\"https://doi.org/10.1093/bioinformatics/btz501\">10.1093/bioinformatics/btz501</a> (erratum: <a href=\"https://doi.org/10.1093/bioinformatics/btz588\">10.1093/bioinformatics/btz588</a>)</li>" : "",
            params.run_mutenricher ? "<li>Soltis, A. R., et al. (2020). MutEnricher: a flexible toolset for somatic mutation enrichment analysis of tumor whole genomes. BMC Bioinformatics, 21. doi: <a href=\"https://doi.org/10.1186/s12859-020-03695-z\">10.1186/s12859-020-03695-z</a></li>" : "",
            "<li>Poole, W., et al. (2016). Combining dependent P-values with an empirical adaptation of Brown's method. Bioinformatics, 32(17). doi: <a href=\"https://doi.org/10.1093/bioinformatics/btw438\">10.1093/bioinformatics/btw438</a></li>",
            params.skip_acat ? "" : "<li>Liu, Y., et al. (2019). ACAT: a fast and powerful p value combination method for rare-variant analysis in sequencing studies. American Journal of Human Genetics, 104(3). doi: <a href=\"https://doi.org/10.1016/j.ajhg.2019.01.002\">10.1016/j.ajhg.2019.01.002</a></li>",
            "<li>Kolde, R., et al. (2012). Robust rank aggregation for gene list integration and meta-analysis. Bioinformatics, 28(4). doi: <a href=\"https://doi.org/10.1093/bioinformatics/btr709\">10.1093/bioinformatics/btr709</a></li>",
            params.run_variant_flags ? "<li>Pedersen, B. S., et al. (2023). Echtvar: compressed variant representation for rapid annotation and filtering of SNPs and indels. Nucleic Acids Research, 51(1). doi: <a href=\"https://doi.org/10.1093/nar/gkac931\">10.1093/nar/gkac931</a></li>" : "",
            params.run_variant_flags ? "<li>Chen, S., et al. (2024). A genomic mutational constraint map using variation in 76,156 human genomes. Nature, 625. doi: <a href=\"https://doi.org/10.1038/s41586-023-06045-0\">10.1038/s41586-023-06045-0</a></li>" : "",
            params.run_vep ? "<li>McLaren, W., et al. (2016). The Ensembl Variant Effect Predictor. Genome Biology, 17. doi: <a href=\"https://doi.org/10.1186/s13059-016-0974-4\">10.1186/s13059-016-0974-4</a></li>" : "",
            params.run_motifbreakr ? "<li>Coetzee, S. G., Coetzee, G. A., & Hazelett, D. J. (2015). motifbreakR: an R/Bioconductor package for predicting variant effects at transcription factor binding sites. Bioinformatics, 31(23). doi: <a href=\"https://doi.org/10.1093/bioinformatics/btv470\">10.1093/bioinformatics/btv470</a></li>" : "",
            params.run_motifbreakr ? "<li>Rauluseviciute, I., et al. (2024). JASPAR 2024: 20th anniversary of the open-access database of transcription factor binding profiles. Nucleic Acids Research, 52(D1). doi: <a href=\"https://doi.org/10.1093/nar/gkad1059\">10.1093/nar/gkad1059</a></li>" : "",
        ] +
        // Step 7A layers: listed only when on, so the text of a run without them is unchanged
        (params.run_alphagenome_avi ? [
            "<li>Avsec, Ž., et al. (2026). Advancing regulatory variant effect prediction with AlphaGenome. Nature, 649(8099), 1206–1218. doi: <a href=\"https://doi.org/10.1038/s41586-025-10014-0\">10.1038/s41586-025-10014-0</a></li>",
            "<li>Cheng, J., et al. (2026). AlphaGenome Atlas: in silico mutagenesis of the entire human genome improves prioritization and interpretation of non-coding variants. medRxiv (preprint). doi: <a href=\"https://doi.org/10.64898/2026.09.16.26363192\">10.64898/2026.09.16.26363192</a></li>"
        ] : []) +
        (params.run_opencravat ? [
            "<li>Pagel, K. A., et al. (2020). Integrated informatics analysis of cancer-related variants. JCO Clinical Cancer Informatics, (4), 310–317. doi: <a href=\"https://doi.org/10.1200/CCI.19.00132\">10.1200/CCI.19.00132</a></li>"
        ] : []) +
        (params.run_cosmic_context ? [
            "<li>Sondka, Z., et al. (2024). COSMIC: a curated database of somatic variants and clinical data for cancer. Nucleic Acids Research, 52(D1). doi: <a href=\"https://doi.org/10.1093/nar/gkad986\">10.1093/nar/gkad986</a></li>",
        ] : []) +
        (params.cancer_gene_list ? [
            "<li>Sondka, Z., et al. (2018). The COSMIC Cancer Gene Census: describing genetic dysfunction across all human cancers. Nature Reviews Cancer, 18(11). doi: <a href=\"https://doi.org/10.1038/s41568-018-0060-1\">10.1038/s41568-018-0060-1</a></li>",
        ] : []) +
        (params.run_enhancer_gene_links ? [
            "<li>Gschwind, A. R., et al. (2026). An encyclopedia of human enhancer-gene regulatory interactions. Nature, 657. doi: <a href=\"https://doi.org/10.1038/s41586-026-10781-4\">10.1038/s41586-026-10781-4</a></li>",
        ] : []) +
        ((params.run_enhancer_gene_links || params.run_ccre_overlay) ? [
            "<li>ENCODE Project Consortium, et al. (2020). Expanded encyclopaedias of DNA elements in the human and mouse genomes. Nature, 583. doi: <a href=\"https://doi.org/10.1038/s41586-020-2493-4\">10.1038/s41586-020-2493-4</a></li>",
        ] : []) +
        (params.run_ccre_overlay ? [
            "<li>Moore, J. E., et al. (2026). An expanded registry of candidate cis-regulatory elements. Nature. doi: <a href=\"https://doi.org/10.1038/s41586-025-09909-9\">10.1038/s41586-025-09909-9</a></li>",
        ] : []) +
        [
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19), 3047–3048. doi: <a href=\"https://doi.org/10.1093/bioinformatics/btw354\">10.1093/bioinformatics/btw354</a></li>"
        ]).join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
