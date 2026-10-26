# NEAREST_TABLE: bedtools closest -d -t first output (element midpoint BED6+2, TSS BED6+2,
# distance) -> one row per significant element. PROXIMITY ONLY: the nearest TSS is not evidence
# that the element regulates that gene. distance_to_tss is signed in the direction of
# transcription (negative = upstream of the TSS), as v0.1: + strand mid - tss, - strand tss - mid.
BEGIN { print "element_id", "chrom", "start", "end", "nearest_gene", "nearest_gene_ensembl_id", "gene_type", "nearest_tss", "gene_strand", "distance_to_tss", "abs_distance_to_tss", "relative_position", "nearest_gene_status" }
{
    if ($12 == "." || $17 == -1) { print $4, $1, $7, $8, "", "", "", "", "", "", "", "", "no_same_chromosome_gene"; next }
    m = $3 + 0; t = $11 + 0; s = $7 + 0; e = $8 + 0
    delta = m - t; signed = ($14 == "+") ? delta : -delta
    rel = (s <= t && t < e) ? "overlapping_tss" : (signed < 0 ? "upstream" : "downstream")
    print $4, $1, s, e, $15, $12, $16, t, $14, signed, (delta < 0 ? -delta : delta), rel, "assigned"
}
