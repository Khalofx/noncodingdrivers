# OPENCRAVAT_INPUT: candidate table -> OpenCRAVAT native ("cravat") input, one line per distinct
# sample x variant: chrom, pos, strand (+), ref, alt, sample. Alleles are passed as in the table
# ('-' for MAF-style indels); OpenCRAVAT normalises them itself and keeps the input coordinates in
# its "Original Input" columns, which OPENCRAVAT_TABLE joins back on (chrom, pos, ref, alt).
NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
{
    k = $c["sample"] SUBSEP $c["chrom"] SUBSEP $c["pos"] SUBSEP $c["ref"] SUBSEP $c["alt"]
    if (k in seen) next
    seen[k] = 1
    print $c["chrom"], $c["pos"], "+", $c["ref"], $c["alt"], $c["sample"]
}
