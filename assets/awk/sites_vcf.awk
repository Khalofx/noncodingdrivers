# SITES_VCF: distinct variants -> sites-only VCF for echtvar anno (gnomAD lookup).
# Format adapter, no statistics. Inputs, in this order:
#   1. FASTA index (.fai): one ##contig line per sequence
#   2. indel anchor bases (bedtools getfasta -nameOnly -tab): key, base
#   3. sites BED4 from variant_sites.awk (mode=sites): chrom, start, end, key
# ID = the canonical key, so the annotation joins back to the table rows.
#   SNV, MNV, COMPLEX : POS pos,   REF ref,          ALT alt
#   deletion          : POS pos-1, REF anchor + ref, ALT anchor
#   insertion         : POS pos,   REF anchor,       ALT anchor + alt
# A variant on a contig absent from the FASTA index (e.g. chrM or an unplaced contig with a
# windowed FASTA) and an indel without an anchor base are not written (bcftools sort needs every
# contig in the header); the counts go to stderr and the germline columns of those rows stay NA
# (variant_flags.awk).
function header() {
    print "##fileformat=VCFv4.2"
    for (i = 1; i <= nc; i++) print contig[i]
    print "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"
    hdr = 1
}
ARGIND == 1 { contig[++nc] = "##contig=<ID=" $1 ",length=" $2 ">"; known[$1] = 1; next }
ARGIND == 2 { anchor[$1] = toupper($2); next }
{
    if (!hdr) header()
    key = $4; n = split(key, k, ":")
    if (n != 4) { print "sites_vcf.awk: malformed key " key > "/dev/stderr"; bad = 1; exit 1 }
    chrom = k[1]; pos = k[2] + 0; ref = k[3]; alt = k[4]
    if (!(chrom in known)) { offcontig++; next }
    if (alt == "-" || ref == "-") {
        if (!(key in anchor)) { skipped++; next }
        a = anchor[key]
        if (alt == "-") print chrom, pos - 1, key, a ref, a, ".", ".", "."
        else            print chrom, pos, key, a, a alt, ".", ".", "."
    } else {
        print chrom, pos, key, ref, alt, ".", ".", "."
    }
}
END {
    if (bad) exit 1
    if (!hdr) header()
    if (offcontig) print "sites_vcf.awk: " offcontig " variant(s) on contigs absent from the FASTA index were not written" > "/dev/stderr"
    if (skipped) print "sites_vcf.awk: " skipped " indel(s) without an anchor base were not written" > "/dev/stderr"
}
