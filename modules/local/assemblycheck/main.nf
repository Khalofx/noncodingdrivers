process ASSEMBLYCHECK {
    tag "${meta.id}"
    label 'process_single'

    // Build check of one samplesheet row before harmonisation (decision D28): EBI
    // vcf_assembly_checker (vcf-validator 0.10.2, Bioconda) compares every REF allele with the
    // FASTA of the declared build. Records on contigs absent from that FASTA are not checked
    // and are counted separately. The task stops with ASSEMBLY_MISMATCH when fewer than
    // min_match of the checked records agree. With fewer than min_records checked records
    // the build cannot be judged and the row passes with status low_evidence (logged).
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vcf-validator:0.10.2--h7a61d87_0' :
        'biocontainers/vcf-validator:0.10.2--h7a61d87_0' }"

    input:
    tuple val(meta), path(vcf)                // VCF or bgzipped VCF; for MAF rows the SNV projection of the MAF
    tuple val(meta2), path(fasta), path(fai)  // FASTA of the declared build and its samtools index
    val(min_match)                            // minimum fraction of matching REF alleles, e.g. 0.99
    val(min_records)                          // fewer checked records: low_evidence, e.g. 50

    output:
    tuple val(meta), path("${prefix}.assembly_check.tsv"),      emit: summary
    tuple val(meta), path("${prefix}.assembly_mismatches.txt"), emit: mismatches
    tuple val("${task.process}"), val('vcf-validator'), eval("vcf_assembly_checker --version 2>&1 | sed -n 's/.*version //p'"), topic: versions, emit: versions_vcfvalidator

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def reader = vcf.name.endsWith('.gz') ? "gzip -cd ${vcf}" : "cat ${vcf}"
    """
    # Records on contigs of the FASTA index only; the others are counted, not checked
    ${reader} \\
        | awk -F'\\t' 'NR == FNR { contig[\$1] = 1; next } /^#/ || (\$1 in contig) { print; next } { off++ } END { print off + 0 > "off_reference.count" }' ${fai} - \\
        > checked_records.vcf

    mkdir -p report
    vcf_assembly_checker \\
        -i checked_records.vcf \\
        -f ${fasta} \\
        -r summary,text \\
        -o report \\
        ${args} \\
        > check.log 2>&1 || true
    cat check.log >&2

    # The full text report lists every mismatch (hundreds of MB for a wrong build): keep the first 1,000
    cat report/*text_assembly_report* 2>/dev/null | head -n 1000 > ${prefix}.assembly_mismatches.txt || true
    rm -rf report checked_records.vcf

    awk -v id='${meta.id}' -v entry='${meta.entry ?: ''}' -v declared='${meta.assembly_declared ?: ''}' \\
        -v reference='${meta2.id ?: fasta.name}' -v min_match='${min_match}' -v min_records='${min_records}' \\
        -v off="\$(cat off_reference.count)" '
        /Number of matches:/ { split(\$NF, a, "/"); matches = a[1] + 0; checked = a[2] + 0; found = 1 }
        END {
            if (!found) {
                print "ASSEMBLY_CHECK_ERROR: vcf_assembly_checker gave no match count for row " id " (see the log above)" > "/dev/stderr"
                exit 2
            }
            frac = checked > 0 ? matches / checked : 0
            status = checked < min_records ? "low_evidence" : (frac < min_match ? "ASSEMBLY_MISMATCH" : "pass")
            OFS = "\\t"
            print "sample", "entry", "assembly_declared", "reference", "records_checked", "records_matching", "match_fraction", "records_off_reference", "min_match", "status"
            print id, entry, declared, reference, checked, matches, sprintf("%.6f", frac), off, min_match, status
            if (status == "low_evidence")
                printf "ASSEMBLY_LOW_EVIDENCE: row %s: only %d SNV records could be checked against %s (fewer than %d); the declared assembly %s is used as given\\n", id, checked, reference, min_records, declared > "/dev/stderr"
            if (status == "ASSEMBLY_MISMATCH") {
                printf "ASSEMBLY_MISMATCH: row %s (entry %s): %d of %d checked REF alleles (%.2f%%) match %s, the FASTA of the declared assembly %s; at least %.0f%% are required. The calls are probably on another build: correct the assembly column of the samplesheet.\\n", id, entry, matches, checked, 100 * frac, reference, declared, 100 * min_match > "/dev/stderr"
                exit 3
            }
        }' check.log > ${prefix}.assembly_check.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'sample\\tentry\\tassembly_declared\\treference\\trecords_checked\\trecords_matching\\tmatch_fraction\\trecords_off_reference\\tmin_match\\tstatus\\n' > ${prefix}.assembly_check.tsv
    printf '${meta.id}\\t${meta.entry ?: ''}\\t${meta.assembly_declared ?: ''}\\tstub\\t0\\t0\\t0\\t0\\t${min_match}\\tlow_evidence\\n' >> ${prefix}.assembly_check.tsv
    touch ${prefix}.assembly_mismatches.txt
    """
}
