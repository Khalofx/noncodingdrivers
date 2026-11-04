# AVI_ELEMENT_SUMMARY: <cohort>.avi.tsv (AVI_TABLE) -> one row per element with candidate variants:
# n_candidate_variants, n_snv (variants the plugin can score), n_scored, max_avi_phred with its
# top_variant (first in table order on ties), and n_phred_ge_20 (AVI_PHRED 20 = top 1 % of all
# possible substitutions by the plugin's definition; a description, not a threshold). Elements
# without a scored variant have max_avi_phred and top_variant NA.
FNR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
{
    m = split($c["element_ids"], el, ",")
    for (j = 1; j <= m; j++) {
        e = el[j]; nv[e]++
        if ($c["avi_status"] != "not_scored_non_snv") nsnv[e]++
        if ($c["avi_status"] == "scored") {
            nsc[e]++; p = $c["avi_phred"] + 0
            if (!(e in mx) || p > mx[e]) { mx[e] = p; mxs[e] = $c["avi_phred"]; top[e] = $c["variant_id"] }
            if (p >= 20) n20[e]++
        }
    }
}
END {
    print "element_id", "n_candidate_variants", "n_snv", "n_scored", "max_avi_phred", "top_variant", "n_phred_ge_20"
    PROCINFO["sorted_in"] = "@ind_str_asc"
    for (e in nv) print e, nv[e], nsnv[e] + 0, nsc[e] + 0, (e in mxs ? mxs[e] : "NA"), (e in top ? top[e] : "NA"), n20[e] + 0
}
