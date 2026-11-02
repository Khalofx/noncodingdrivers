# RE2G_SUMMARY: one row per significant element x requested tissue, from
#   <cohort>.significant.sorted.bed     element fragments (the element universe)
#   <cohort>.enhancer_gene_links.tsv    RE2G_LINKS
#   <cohort>.nearest_gene.tsv           NEAREST_TABLE (optional)
# -v tissues=a,b  -v n_files=a:11,b:52  (staged files per tissue)
# Self-promoter links (isSelfPromoter TRUE: the element is the gene's own promoter) are kept
# apart from distal links. A target's score is its maximum over the tissue's files; ties are
# broken by gene name; top scores are printed as in the file, the target lists to 3 decimals. nearest_gene_status compares the nearest gene (symbol or Ensembl id
# without version) with the distal targets, nearest_gene_status_any with all targets:
# top_is_nearest, nearest_in_targets, nearest_not_in_targets, no_link, or no_nearest_gene.
BEGIN {
    nt = split(tissues, tlist, ",")
    m = split(n_files, nf, ",")
    for (i = 1; i <= m; i++) { split(nf[i], kv, ":"); files[kv[1]] = kv[2] }
}
function strip(id) { sub(/\.[0-9]+$/, "", id); return id }
FNR == 1 {
    kind = FILENAME ~ /enhancer_gene_links\.tsv$/ ? "links" : (FILENAME ~ /nearest_gene\.tsv$/ ? "nearest" : "fragments")
    if (kind != "fragments") { for (i = 1; i <= NF; i++) col[kind, $i] = i; next }
}
kind == "fragments" { element[$4] = 1; next }
kind == "nearest" { has_nearest = 1; ng[$col["nearest", "element_id"]] = $col["nearest", "nearest_gene"]; ne[$col["nearest", "element_id"]] = $col["nearest", "nearest_gene_ensembl_id"]; next }
kind == "links" {
    e = $col["links", "element_id"]; t = $col["links", "tissue"]; g = $col["links", "target_gene"]
    ss = $col["links", "re2g_score"]; s = ss + 0; sp = $col["links", "is_self_promoter"] == "TRUE"
    nrow[e, t]++; withlink[e, t, $col["links", "accession"]] = 1
    gid[g] = strip($col["links", "target_gene_ensembl_id"])
    if (!((e, t, g) in best) || s > best[e, t, g]) { best[e, t, g] = s; bests[e, t, g] = ss }
    if (sp) selfp[e, t, g] = 1
    else if (!((e, t, g) in dist) || s > dist[e, t, g]) { dist[e, t, g] = s; dists[e, t, g] = ss }
}
# targets of (e, t) from array A, sorted by score desc then gene name, into out[]; returns count
function ranked(A, e, t, out,    k, g, n, i, j, tmp, p) {
    n = 0
    for (k in A) { split(k, p, SUBSEP); if (p[1] == e && p[2] == t) { n++; out[n] = p[3] } }
    for (i = 2; i <= n; i++) {
        tmp = out[i]; j = i - 1
        while (j > 0 && (A[e, t, out[j]] < A[e, t, tmp] || (A[e, t, out[j]] == A[e, t, tmp] && out[j] > tmp))) { out[j + 1] = out[j]; j-- }
        out[j + 1] = tmp
    }
    return n
}
function is_nearest(g, e) { return (ng[e] != "" && g == ng[e]) || (ne[e] != "" && gid[g] == strip(ne[e])) }
function status(n, out, e,    i) {
    if (!has_nearest || ng[e] == "" && ne[e] == "") return "no_nearest_gene"
    if (n == 0) return "no_link"
    if (is_nearest(out[1], e)) return "top_is_nearest"
    for (i = 2; i <= n; i++) if (is_nearest(out[i], e)) return "nearest_in_targets"
    return "nearest_not_in_targets"
}
END {
    print "element_id", "tissue", "nearest_gene", "nearest_gene_ensembl_id", "n_files", "n_files_with_link", "n_link_rows", "self_promoter_genes", "n_distal_targets", "distal_targets_by_max_score", "top_distal_target", "top_distal_score", "n_targets", "top_target_any", "top_any_score", "nearest_gene_status", "nearest_gene_status_any", "threshold", "threshold_source"
    PROCINFO["sorted_in"] = "@ind_str_asc"
    for (e in element) {
        for (ti = 1; ti <= nt; ti++) {
            t = tlist[ti]
            delete D; delete A
            nd = ranked(dist, e, t, D)
            na = ranked(best, e, t, A)
            nw = 0
            for (k in withlink) { split(k, p, SUBSEP); if (p[1] == e && p[2] == t) nw++ }
            sps = ""
            for (k in selfp) { split(k, p, SUBSEP); if (p[1] == e && p[2] == t) sps = sps (sps == "" ? "" : ",") p[3] }
            dl = ""
            for (i = 1; i <= nd; i++) dl = dl (i > 1 ? "," : "") sprintf("%s:%.3f", D[i], dist[e, t, D[i]])
            print e, t, (has_nearest ? ng[e] : "NA"), (has_nearest ? ne[e] : "NA"), (t in files ? files[t] : "NA"), nw, nrow[e, t] + 0, sps, nd, dl,
                (nd ? D[1] : ""), (nd ? dists[e, t, D[1]] : ""), na, (na ? A[1] : ""), (na ? bests[e, t, A[1]] : ""),
                status(nd, D, e), status(na, A, e), "NA", "ENCODE-rE2G thresholded predictions (threshold at 70% recall of CRISPR-validated links; not in the file)"
        }
    }
}
