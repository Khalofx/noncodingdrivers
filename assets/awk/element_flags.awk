# ELEMENT_FLAGS_TABLE: element-level support and mask fractions, one row per element id.
# Format adapter: counts and fractions only. Inputs, in this order:
#   1. bedtools coverage of the catalogue fragments by the callable territory (Umap k100 minus
#      ENCODE Blacklist v2): chrom start end element_id count covered length fraction
#   2. bedtools intersect -wa -wb of FLAGS_BED with the catalogue: 10 flag columns + chrom start end element_id
# A row (sample, variant) that touches two fragments of one element counts once.
#   n_mutations, n_patients, max_sample_fraction (largest share of one patient)
#   frac_in_blacklist_v2, frac_in_segdup, frac_in_low_mappability: share of the element's mutations
#   mean_umap_unique_fraction; frac_germline_pop_common (over rows with a germline result)
#   element_bp, callable_bp, callable_fraction: element territory in callable sequence
# Elements without a mutation keep their territory columns; their mutation columns are 0 or NA.
function frac(a, b) { return b > 0 ? sprintf("%.6g", a / b) : "NA" }
ARGIND == 1 {
    id = $4
    if (!(id in bp)) order[++n_el] = id
    bp[id] += $7; cbp[id] += $6
    next
}
{
    id = $14; row = $4 SUBSEP $5
    if ((id, row) in seen) next
    seen[id, row] = 1
    nm[id]++
    if (!((id, $4) in pat)) { pat[id, $4] = 0; np[id]++ }
    pat[id, $4]++
    if (pat[id, $4] > maxp[id]) maxp[id] = pat[id, $4]
    if ($6 == "true") nbl[id]++
    if ($7 == "true") nsd[id]++
    if ($8 == "true") nlm[id]++
    um[id] += $9
    if ($10 == "true" || $10 == "false") { ng[id]++; if ($10 == "true") ngc[id]++ }
}
END {
    print "element_id", "n_mutations", "n_patients", "max_sample_fraction", "frac_in_blacklist_v2", "frac_in_segdup",
          "frac_in_low_mappability", "mean_umap_unique_fraction", "frac_germline_pop_common",
          "element_bp", "callable_bp", "callable_fraction"
    for (i = 1; i <= n_el; i++) {
        id = order[i]; m = nm[id] + 0
        print id, m, np[id] + 0, frac(maxp[id], m), frac(nbl[id], m), frac(nsd[id], m), frac(nlm[id], m),
              frac(um[id], m), frac(ngc[id], ng[id]), bp[id], cbp[id], frac(cbp[id], bp[id])
    }
}
