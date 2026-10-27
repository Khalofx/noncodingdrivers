# Format adapter (no statistics): bcftools query output -> canonical mutation table.
# Input lines: CHROM POS REF ALT (one ALT per record; bcftools norm -m -any upstream).
# Pass the sample name with -v sample=<id> (conf/modules.config).
#
# VCF alleles are anchored; the canonical table uses the MAF conventions of the mutations tier:
#   SNV / MNV : pos, REF, ALT unchanged; vtype SNP, DNP, TNP or ONP by length
#   deletion  : REF=AT ALT=A at p   -> pos p+1, ref T, alt '-', vtype DEL
#   insertion : REF=A ALT=AT at p   -> pos p,   ref '-', alt T, vtype INS
#   anything else (complex length change without a shared first base): pos, REF, ALT
#   unchanged, vtype COMPLEX (counted on stderr; callers drop non-SNV rows they cannot use).
# Records with a missing or symbolic ALT ('.', '*', '<...>') are skipped and counted.
BEGIN {
    FS = OFS = "\t"
    if (sample == "") { print "QUERY_SAMPLE_MISSING: run with -v sample=<id>" > "/dev/stderr"; exit 1 }
    print "sample", "chrom", "pos", "ref", "alt", "vtype"
}
{
    chrom = $1; pos = $2; ref = toupper($3); alt = toupper($4)
    if (alt == "." || alt == "*" || alt ~ /^</ || alt ~ /[^ACGTN]/ || ref ~ /[^ACGTN]/) { skipped++; next }
    lr = length(ref); la = length(alt)
    if (lr == la) {
        vtype = (lr == 1) ? "SNP" : (lr == 2) ? "DNP" : (lr == 3) ? "TNP" : "ONP"
    } else if (substr(ref, 1, 1) == substr(alt, 1, 1) && (lr == 1 || la == 1)) {
        if (lr > la) { pos = pos + 1; ref = substr(ref, 2); alt = "-"; vtype = "DEL" }
        else         { ref = "-"; alt = substr(alt, 2); vtype = "INS" }
    } else {
        vtype = "COMPLEX"; complex++
    }
    print sample, chrom, pos, ref, alt, vtype
    n++
}
END {
    printf "query_to_canonical: %s: %d rows written, %d complex, %d skipped\n", sample, n, complex, skipped > "/dev/stderr"
}
