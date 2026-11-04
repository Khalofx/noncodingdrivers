# RE2G_LINKS: per-file bedtools intersect -wo outputs (element fragment BED4, ENCODE-rE2G link
# record, overlap bp) -> one row per element x link file x rE2G element x target gene.
# Input files: <cohort>.<tissue>.<biosample>.<accession>.re2g.bed (prefix set in modules.config).
# rE2G thresholded files (DNaseOnly model): columns 1-10 are chr start end name class TargetGene
# TargetGeneEnsemblID TargetGeneTSS isSelfPromoter CellType and the last one is Score (checked
# by the staging recipe); here they start at column 5, Score is $(NF-1) and the overlap $NF.
# Overlaps of several fragments of one element with the same link are summed.
FNR == 1 {
    n = split(FILENAME, a, "/"); base = a[n]
    n = split(base, p, "."); acc = p[n - 2]; bios = p[n - 3]; tissue = p[n - 4]
}
{
    cls = $9; sp = $13; score = $(NF - 1)
    if (!(cls == "promoter" || cls == "genic" || cls == "intergenic") || !(sp == "TRUE" || sp == "FALSE") || score !~ /^[0-9.eE+-]+$/) {
        printf "re2g_links.awk: %s line %d is not an ENCODE-rE2G thresholded record (class %s, isSelfPromoter %s, Score %s)\n", FILENAME, FNR, cls, sp, score > "/dev/stderr"
        exit 1
    }
    k = $4 SUBSEP tissue SUBSEP bios SUBSEP acc SUBSEP $5 ":" $6 "-" $7 SUBSEP $10
    if (!(k in bp)) { klass[k] = cls; ens[k] = $11; tss[k] = $12; self[k] = sp; sc[k] = score }
    bp[k] += $NF
}
function by_link(i1, v1, i2, v2,    a, b, c) {
    split(i1, a, SUBSEP); split(i2, b, SUBSEP)
    if (a[1] != b[1]) return a[1] < b[1] ? -1 : 1                       # element
    if (a[2] != b[2]) return a[2] < b[2] ? -1 : 1                       # tissue
    if (sc[i1] + 0 != sc[i2] + 0) return sc[i1] + 0 > sc[i2] + 0 ? -1 : 1   # score, descending
    for (c = 6; c >= 3; c--) if (a[c] != b[c]) return a[c] < b[c] ? -1 : 1   # gene, link, accession, biosample
    return 0
}
END {
    print "element_id", "tissue", "biosample", "accession", "re2g_element", "re2g_class", "target_gene", "target_gene_ensembl_id", "target_gene_tss", "is_self_promoter", "re2g_score", "overlap_bp"
    PROCINFO["sorted_in"] = "by_link"
    for (k in bp) {
        split(k, f, SUBSEP)
        print f[1], f[2], f[3], f[4], f[5], klass[k], f[6], ens[k], tss[k], self[k], sc[k], bp[k]
    }
}
