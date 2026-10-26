# MUTATION_BED: canonical cohort mutation table -> BED6+3
#   chrom, pos-1, pos, sample, 0, ., ref, alt, vtype
# One base per mutation at its table position (indels too: the point rule of v0.1, which
# assigned a mutation to an element when its 1-based position lies in the BED interval).
NR == 1 {
    for (i = 1; i <= NF; i++) col[$i] = i
    split("sample chrom pos ref alt vtype", need, " ")
    for (k in need) if (!(need[k] in col)) { print "mutation_bed.awk: table lacks " need[k] > "/dev/stderr"; exit 1 }
    next
}
{ print $col["chrom"], $col["pos"] - 1, $col["pos"], $col["sample"], 0, ".", $col["ref"], $col["alt"], $col["vtype"] }
