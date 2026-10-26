# Format adapter (no statistics): split one MAF into one MAF per Tumor_Sample_Barcode for
# maf2vcf.pl. Run through the nf-core gawk module with output redirection disabled; every
# sample is written to <barcode>.sample.maf (characters outside [A-Za-z0-9._-] become '_',
# in the file name and in the barcode columns).
#
#  - strips CR line endings and '#' comment lines
#  - renames header fields case-insensitively to the GDC MAF spelling (PCAWG writes
#    Start_position / End_position; maf2vcf.pl expects Start_Position / End_Position)
#  - writes UCSC contig names (chr1..chr22, chrX, chrY): the declared-build FASTA used by
#    maf2vcf.pl (hg19 UCSC or the GRCh38 analysis set) is chr-named. Other contigs
#    (chrM, unplaced, decoys) are dropped and counted on stderr.
#
# With -v mode=check the program writes no files. It prints the SNV rows (single-base
# Reference_Allele and Tumor_Seq_Allele2, A/C/G/T) as a sites-only VCF on stdout, with the
# MAF reference allele as REF and no FASTA lookup. That VCF is what the assembly check
# (vcf_assembly_checker) compares with the declared build: maf2vcf.pl drops records whose
# reference allele disagrees with its FASTA, so its output cannot reveal a build mismatch.
BEGIN {
    FS = OFS = "\t"
    n = split("Hugo_Symbol Chromosome Start_Position End_Position Variant_Type Reference_Allele " \
              "Tumor_Seq_Allele1 Tumor_Seq_Allele2 Tumor_Sample_Barcode Matched_Norm_Sample_Barcode " \
              "t_ref_count t_alt_count t_depth n_ref_count n_alt_count n_depth", canon, " ")
    for (i = 1; i <= n; i++) spell[tolower(canon[i])] = canon[i]
}
{ sub(/\r$/, "") }
/^#/ || /^[[:space:]]*$/ { next }
!have_header {
    for (i = 1; i <= NF; i++) {
        key = tolower($i)
        if (key in spell) $i = spell[key]
        col[$i] = i
    }
    split("Chromosome Start_Position End_Position Reference_Allele Tumor_Seq_Allele2 Tumor_Sample_Barcode", req, " ")
    for (r in req) if (!(req[r] in col)) {
        printf "MAF_COLUMN_MISSING: %s not in the MAF header of %s\n", req[r], FILENAME > "/dev/stderr"
        failed = 1; exit 1
    }
    header = $0; have_header = 1
    if (mode == "check") {
        print "##fileformat=VCFv4.2"
        print "##source=noncodingdrivers maf_split.awk (mode=check): MAF SNV rows, REF = Reference_Allele"
        print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    }
    next
}
{
    chrom = $col["Chromosome"]; sub(/^chr/, "", chrom)
    if (chrom !~ /^([1-9]|1[0-9]|2[0-2]|X|Y)$/) { dropped[chrom]++; n_dropped++; next }
    $col["Chromosome"] = "chr" chrom
    if (mode == "check") {
        ref = toupper($col["Reference_Allele"]); alt = toupper($col["Tumor_Seq_Allele2"])
        if (ref ~ /^[ACGT]$/ && alt ~ /^[ACGT]$/ && ref != alt) { print "chr" chrom, $col["Start_Position"], ".", ref, alt, ".", ".", "."; n_kept++ }
        else n_other++
        next
    }
    barcode = $col["Tumor_Sample_Barcode"]
    safe = barcode; gsub(/[^A-Za-z0-9._-]/, "_", safe)
    if (safe == "") { printf "MAF_SAMPLE_MISSING: empty Tumor_Sample_Barcode at line %d of %s\n", FNR, FILENAME > "/dev/stderr"; failed = 1; exit 1 }
    if ((safe in owner) && owner[safe] != barcode) {
        printf "MAF_SAMPLE_COLLISION: '%s' and '%s' map to the same file name %s\n", owner[safe], barcode, safe > "/dev/stderr"
        failed = 1; exit 1
    }
    # maf2vcf.pl refuses barcodes with '/' or blanks: the rows carry the file-name-safe form,
    # which is also the sample name of the canonical table
    $col["Tumor_Sample_Barcode"] = safe
    if ("Matched_Norm_Sample_Barcode" in col) gsub(/[^A-Za-z0-9._-]/, "_", $col["Matched_Norm_Sample_Barcode"])
    out = safe ".sample.maf"
    if (!(safe in owner)) { owner[safe] = barcode; print header > out; n_samples++ }
    print > out
    n_kept++
}
END {
    if (failed) exit 1
    if (!have_header) { printf "MAF_EMPTY: no header line in %s\n", FILENAME > "/dev/stderr"; exit 1 }
    if (mode == "check") {
        printf "maf_split (check): %d SNV rows written, %d other rows (indels, MNVs) not checked, %d rows on other contigs dropped\n", n_kept, n_other, n_dropped > "/dev/stderr"
        exit 0
    }
    printf "maf_split: %d rows kept in %d samples; %d rows on other contigs dropped\n", n_kept, n_samples, n_dropped > "/dev/stderr"
    for (k in dropped) printf "maf_split: dropped contig %s: %d rows\n", k, dropped[k] > "/dev/stderr"
}
