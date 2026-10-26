# MIDPOINT_BED: consensus table -> one base at the midpoint of each significant element span,
#   chrom, mid-1, mid, element_id, 0, ., span_start, span_end
# with mid = floor((start + end) / 2) on the BED span, as v0.1 (evidence_cards._nearest_gene).
NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
$col["is_consensus"] == "True" {
    s = $col["start"] + 0; e = $col["end"] + 0; m = int((s + e) / 2)
    print $col["chrom"], m - 1, m, $col["element_id"], 0, ".", s, e
}
