# GENE_TSS: GENCODE GTF -> sorted TSS BED6+2 of every `gene` record (all biotypes)
#   chrom, tss-1, tss, gene_id, 0, strand, gene_name, gene_type
# A gene's TSS is its start on + and its end on - (the v0.1 strand fix, d994e30). Rows are
# sorted by chrom (byte order, as `sort -k1,1`), TSS and gene_id, so `bedtools closest -t first`
# resolves genes sharing one TSS by the lowest gene_id.
BEGIN { FS = "\t"; OFS = "\t" }
/^#/ { next }
$3 == "gene" {
    chrom = ($1 ~ /^chr/) ? $1 : "chr" $1
    tss = ($7 == "+") ? $4 : $5
    id = ""; name = ""; type = ""
    if (match($9, /gene_id "[^"]+"/))   id   = substr($9, RSTART + 9,  RLENGTH - 10)
    if (match($9, /gene_name "[^"]+"/)) name = substr($9, RSTART + 11, RLENGTH - 12)
    if (match($9, /gene_type "[^"]+"/)) type = substr($9, RSTART + 11, RLENGTH - 12)
    if (name == "") name = (id == "" ? "NA" : id)
    if (id == "") id = "NA"
    n++; c[n] = chrom; t[n] = tss; row[n] = chrom OFS (tss - 1) OFS tss OFS id OFS 0 OFS $7 OFS name OFS type; g[n] = id
}
function cmp(i1, v1, i2, v2) {
    if (c[v1] != c[v2]) return c[v1] < c[v2] ? -1 : 1
    if (t[v1] != t[v2]) return t[v1] < t[v2] ? -1 : 1
    if (g[v1] != g[v2]) return g[v1] < g[v2] ? -1 : 1
    return 0
}
END {
    for (i = 1; i <= n; i++) idx[i] = i
    m = asort(idx, ord, "cmp")
    for (i = 1; i <= m; i++) print row[ord[i]]
}
