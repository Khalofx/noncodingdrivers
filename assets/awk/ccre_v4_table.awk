# CCRE_TABLE / CCRE_SUMMARY: bedtools intersect -wao of the element fragments (BED4) with the
# SCREEN Registry V4 agnostic cCREs (BED9+1: accession in column 4, class in column 10), so the
# cCRE accession is $8, its class $14 and the overlap $NF; fragments without a cCRE carry '.'.
#   -v mode=table    one row per element x cCRE (overlap summed over the element's fragments)
#   -v mode=summary  one row per element: element and covered bases (union of the cCREs within
#                    the fragments), covered fraction, number of cCREs, class with most bases
#                    (ties by name) and bases per class
{
    e = $4; frag = $1 SUBSEP $2 SUBSEP $3
    if (!((e, frag) in seenfrag)) { seenfrag[e, frag] = 1; ebp[e] += $3 - $2; element[e] = 1 }
    if ($8 == "." || $NF + 0 == 0) next
    if ($14 == "" || $14 == ".") { printf "ccre_v4_table.awk: line %d has no cCRE class in column 10 of the V4 file\n", NR > "/dev/stderr"; exit 1 }
    ov[e, $8] += $NF; cls[e, $8] = $14; clsbp[e, $14] += $NF
    s = ($6 > $2 ? $6 : $2); t = ($7 < $3 ? $7 : $3)
    ni[e]++; is[e, ni[e]] = s; ie[e, ni[e]] = t; ichr[e, ni[e]] = $1
}
function covered(e,    i, j, k, n, tmp, cs, ce, cc, tot, ord) {
    n = ni[e]; tot = 0
    for (i = 1; i <= n; i++) { ord[i] = i }
    for (i = 2; i <= n; i++) {
        tmp = ord[i]; j = i - 1
        while (j > 0 && (ichr[e, ord[j]] > ichr[e, tmp] || (ichr[e, ord[j]] == ichr[e, tmp] && is[e, ord[j]] > is[e, tmp]))) { ord[j + 1] = ord[j]; j-- }
        ord[j + 1] = tmp
    }
    cc = ""
    for (i = 1; i <= n; i++) {
        k = ord[i]
        if (ichr[e, k] != cc || is[e, k] > ce) { if (cc != "") tot += ce - cs; cc = ichr[e, k]; cs = is[e, k]; ce = ie[e, k] }
        else if (ie[e, k] > ce) ce = ie[e, k]
    }
    if (cc != "") tot += ce - cs
    return tot
}
END {
    PROCINFO["sorted_in"] = "@ind_str_asc"
    if (mode == "table") {
        print "element_id", "ccre_accession", "ccre_class", "overlap_bp", "element_bp"
        for (k in ov) { split(k, f, SUBSEP); print f[1], f[2], cls[k], ov[k], ebp[f[1]] }
    } else if (mode == "summary") {
        print "element_id", "element_bp", "n_ccre", "covered_bp", "covered_frac", "dominant_class", "ccre_v4_classes_bp"
        for (e in element) {
            n = 0; for (k in ov) { split(k, f, SUBSEP); if (f[1] == e) n++ }
            delete C; nc = 0
            for (k in clsbp) { split(k, f, SUBSEP); if (f[1] == e) { nc++; C[nc] = f[2] } }
            for (i = 2; i <= nc; i++) { tmp = C[i]; j = i - 1; while (j > 0 && (clsbp[e, C[j]] < clsbp[e, tmp] || (clsbp[e, C[j]] == clsbp[e, tmp] && C[j] > tmp))) { C[j + 1] = C[j]; j-- } C[j + 1] = tmp }
            lst = ""; for (i = 1; i <= nc; i++) lst = lst (i > 1 ? "," : "") C[i] ":" clsbp[e, C[i]]
            cb = covered(e)
            print e, ebp[e], n, cb, sprintf("%.4f", ebp[e] ? cb / ebp[e] : 0), (nc ? C[1] : "none"), lst
        }
    } else {
        print "ccre_v4_table.awk: set -v mode=table or -v mode=summary" > "/dev/stderr"; exit 1
    }
}
