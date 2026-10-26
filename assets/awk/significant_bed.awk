# SIGNIFICANT_BED: consensus table + element catalogue BED -> BED4 of the FRAGMENTS of every
# element with is_consensus == True. Fragments, not the element span: the span of a fragmented
# element also covers the excluded coding bases between its fragments, which no caller tested.
# Inputs, in this order: <cohort>.consensus.tsv  <elements>.bed (headerless: chrom start end id)
FNR == 1 { file++ }
file == 1 && FNR == 1 {
    for (i = 1; i <= NF; i++) col[$i] = i
    if (!("element_id" in col) || !("is_consensus" in col)) {
        print "significant_bed.awk: consensus table lacks element_id or is_consensus" > "/dev/stderr"
        exit 1
    }
    next
}
file == 1 { if ($col["is_consensus"] == "True") sig[$col["element_id"]] = 1; next }
($4 in sig) { print $1, $2, $3, $4 }
