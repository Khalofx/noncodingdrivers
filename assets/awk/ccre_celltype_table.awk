# CCRE_CELLTYPE_TABLE: per-file bedtools intersect -wo outputs of the element fragments (BED4)
# with SCREEN V4 biosample-specific cCRE files (BED9+2: accession in column 4, the biosample's
# class in column 10) -> one row per element x file x cCRE, overlap summed over fragments.
# Input files: <cohort>.<tissue>.<biosample>.<accession>.ccre_celltype.bed
FNR == 1 {
    n = split(FILENAME, a, "/"); base = a[n]
    n = split(base, p, "."); acc = p[n - 2]; bios = p[n - 3]; tissue = p[n - 4]
}
{ k = $4 SUBSEP tissue SUBSEP bios SUBSEP acc SUBSEP $8; ov[k] += $NF; cls[k] = $14 }
END {
    print "element_id", "tissue", "biosample", "accession", "ccre_accession", "celltype_class", "overlap_bp"
    PROCINFO["sorted_in"] = "@ind_str_asc"
    for (k in ov) { split(k, f, SUBSEP); print f[1], f[2], f[3], f[4], f[5], cls[k], ov[k] }
}
