# SIGPROFILER_INPUT: canonical cohort mutation table -> the columns the nf-core sigprofiler
# module reads (Indiv chr from to ref alt NV), single-base substitutions only.
NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; print "Indiv", "chr", "from", "to", "ref", "alt", "NV"; next }
length($col["ref"]) == 1 && length($col["alt"]) == 1 && $col["ref"] != "-" && $col["alt"] != "-" {
    print $col["sample"], $col["chrom"], $col["pos"], $col["pos"], $col["ref"], $col["alt"], 1
}
