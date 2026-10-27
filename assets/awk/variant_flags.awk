# VARIANT_FLAGS_TABLE: per-row variant flags of the cohort mutation table (flag, never filter).
# Format adapter: joins the tool outputs on the variant key and applies the germline rule.
# Inputs, in this order (keys: chrom:pos:ref:alt, from variant_sites.awk):
#   1. bedtools intersect -c, ENCODE Blacklist v2 : chrom start end key count
#   2. bedtools intersect -c, GIAB segdups        : chrom start end key count
#   3. bedtools intersect -c, GIAB low mappability: chrom start end key count
#   4. bedtools coverage, Umap k100               : chrom start end key count covered length fraction
#   5. (optional) bcftools query of the echtvar-annotated sites: key gnomad_filter gnomad_af gnomad_fafmax_faf95_max
#   last. the cohort mutation table (header: sample chrom pos ref alt vtype ...)
# germline_pop_common = true when gnomAD FILTER is PASS and fafmax_faf95_max > faf_threshold
# (-v faf_threshold=..., --germline_faf_threshold). Without input 5 the germline columns are NA;
# a variant absent from the sites VCF (indel without anchor base) is NA; a variant gnomAD does not
# carry is false. echtvar writes gnomad_filter=MISSING and -1 for the numbers of a variant gnomAD
# does not carry, and joins several FILTER values with commas (AC0,AS_VQSR); MISSING and -1 are
# written here as NA.
BEGIN {
    if (faf_threshold == "") { print "variant_flags.awk: run with -v faf_threshold=<value>" > "/dev/stderr"; bad = 1; exit 1 }
    if (ARGC - 1 != 5 && ARGC - 1 != 6) { print "variant_flags.awk: expected 5 or 6 input files, got " ARGC - 1 > "/dev/stderr"; bad = 1; exit 1 }
    germline = (ARGC - 1 == 6)
}
ARGIND < ARGC - 1 {
    if (ARGIND <= 3) { flag[ARGIND, $4] = ($5 > 0) ? "true" : "false"; next }
    if (ARGIND == 4) { umap[$4] = $8; next }
    if (ARGIND == 5) { gfilter[$1] = $2; gaf[$1] = $3; gfaf[$1] = $4; next }
}
function na(v) { return (v == "." || v == "" || v + 0 < 0) ? "NA" : v }
FNR == 1 && ARGIND == ARGC - 1 {
    for (i = 1; i <= NF; i++) col[$i] = i
    split("sample chrom pos ref alt vtype", need, " ")
    for (k in need) if (!(need[k] in col)) { print "variant_flags.awk: table lacks " need[k] > "/dev/stderr"; bad = 1; exit 1 }
    print "sample", "chrom", "pos", "ref", "alt", "vtype", "in_blacklist_v2", "in_segdup", "in_low_mappability",
          "umap_unique_fraction", "gnomad_filter", "gnomad_af", "gnomad_fafmax_faf95_max", "germline_pop_common"
    next
}
{
    key = $col["chrom"] ":" $col["pos"] ":" $col["ref"] ":" $col["alt"]
    if (!((1, key) in flag)) { print "variant_flags.awk: no mask result for " key > "/dev/stderr"; bad = 1; exit 1 }
    if (!germline || !(key in gfilter)) { f = "NA"; af = "NA"; faf = "NA"; common = "NA" }
    else {
        f = (gfilter[key] == "MISSING" || gfilter[key] == "." || gfilter[key] == "") ? "NA" : gfilter[key]
        af = na(gaf[key]); faf = na(gfaf[key])
        common = (f == "PASS" && faf != "NA" && faf + 0 > faf_threshold + 0) ? "true" : "false"
    }
    print $col["sample"], $col["chrom"], $col["pos"], $col["ref"], $col["alt"], $col["vtype"],
          flag[1, key], flag[2, key], flag[3, key], umap[key], f, af, faf, common
}
END { if (bad) exit 1 }
