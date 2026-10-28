# VARIANT_SITES / INDEL_ANCHORS: canonical cohort mutation table -> BED4 of distinct variants.
# Format adapter, no statistics. Run with -v mode=sites or -v mode=anchors (conf/modules.config).
#   key = chrom:pos:ref:alt (the canonical, MAF-style alleles); each key is written once.
#   mode=sites   : the variant's span on GRCh38, for the mask intersections
#                  SNV, MNV, COMPLEX, deletion: [pos-1, pos-1+length(ref))
#                  insertion (ref '-'): the base before the insertion, [pos-1, pos)
#   mode=anchors : the base the VCF form of an indel is anchored on (bedtools getfasta)
#                  deletion: the base before the first deleted base, [pos-2, pos-1)
#                  insertion: the base at pos, [pos-1, pos)
BEGIN {
    if (mode != "sites" && mode != "anchors") { print "variant_sites.awk: run with -v mode=sites|anchors" > "/dev/stderr"; bad = 1; exit 1 }
}
NR == 1 {
    for (i = 1; i <= NF; i++) col[$i] = i
    split("chrom pos ref alt", need, " ")
    for (k in need) if (!(need[k] in col)) { print "variant_sites.awk: table lacks " need[k] > "/dev/stderr"; bad = 1; exit 1 }
    next
}
{
    chrom = $col["chrom"]; pos = $col["pos"] + 0; ref = $col["ref"]; alt = $col["alt"]
    key = chrom ":" pos ":" ref ":" alt
    if (seen[key]++) next
    if (mode == "sites") {
        if (ref == "-") print chrom, pos - 1, pos, key
        else            print chrom, pos - 1, pos - 1 + length(ref), key
    } else if (alt == "-" && pos > 1) {
        print chrom, pos - 2, pos - 1, key
    } else if (ref == "-") {
        print chrom, pos - 1, pos, key
    }
}
