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
        false
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
                multiqc_reports.getVal(),
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
    ]
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
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

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

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
