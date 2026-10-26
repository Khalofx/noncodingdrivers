# REGCTX_TABLE: per-peak-file bedtools outputs -> long table, one row per significant element x
# staged tissue x layer. A layer that was not staged for a tissue has no row (absent, never 0).
# Inputs: <cohort>.<tissue>.<layer>.<accession>.coverage.bed  (bedtools coverage: A + count, bp, len, frac)
#         <cohort>.<tissue>.<layer>.<accession>.signal-c<N>.bed (bedtools map -o max on column N)
# Rows of a fragmented element are summed over its fragments (covered bp, length, peak count).
# max_signal is signalValue (narrowPeak column 7, fold enrichment) or max_density (DNase hotspot
# column 5): comparable WITHIN a layer only.
FNR == 1 {
    n = split(FILENAME, a, "/"); base = a[n]
    n = split(base, p, "."); kind = p[n - 1]; acc = p[n - 2]; layer = p[n - 3]; tissue = p[n - 4]
}
{
    k = $4 SUBSEP tissue SUBSEP layer SUBSEP acc
    if (kind == "coverage") { seen[k] = 1; npk[k] += $5; cov[k] += $6; len[k] += $7; nfr[k]++ }
    else if (kind ~ /^signal-c/) {
        src[k] = (kind == "signal-c5") ? "max_density(col5)" : "signalValue(col" substr(kind, 9) ")"
        if ($5 != "." && (!(k in mx) || $5 + 0 > mx[k] + 0)) mx[k] = $5
    }
}
END {
    print "element_id", "tissue", "layer", "accession", "n_fragments", "n_peaks", "covered_bp", "element_bp", "covered_frac", "max_signal", "signal_source", "status"
    PROCINFO["sorted_in"] = "@ind_str_asc"
    for (k in seen) {
        split(k, f, SUBSEP)
        print f[1], f[2], f[3], f[4], nfr[k], npk[k], cov[k], len[k], sprintf("%.4f", len[k] ? cov[k] / len[k] : 0), (k in mx ? mx[k] : "NA"), (k in src ? src[k] : "NA"), (npk[k] > 0 ? "overlap" : "no_overlap")
    }
}
