# FLAGS_BED: per-row variant flags -> BED for the element assignment (bedtools intersect -wa -wb).
# Format adapter, no statistics. One base per row at its table position: the point rule that
# assigns a mutation to an element (mutation_bed.awk).
#   chrom, pos-1, pos, sample, key, in_blacklist_v2, in_segdup, in_low_mappability,
#   umap_unique_fraction, germline_pop_common
NR == 1 {
    for (i = 1; i <= NF; i++) col[$i] = i
    split("sample chrom pos ref alt in_blacklist_v2 in_segdup in_low_mappability umap_unique_fraction germline_pop_common", need, " ")
    for (k in need) if (!(need[k] in col)) { print "flags_bed.awk: table lacks " need[k] > "/dev/stderr"; exit 1 }
    next
}
{
    print $col["chrom"], $col["pos"] - 1, $col["pos"], $col["sample"],
          $col["chrom"] ":" $col["pos"] ":" $col["ref"] ":" $col["alt"],
          $col["in_blacklist_v2"], $col["in_segdup"], $col["in_low_mappability"],
          $col["umap_unique_fraction"], $col["germline_pop_common"]
}
