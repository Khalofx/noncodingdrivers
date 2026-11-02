# OPENCRAVAT_TABLE: candidate table + OpenCRAVAT text report (variant level) -> one row per distinct
# candidate variant with the OpenCRAVAT columns. Inputs (in this order): <cohort>.candidates.tsv,
# <cohort>.opencravat.report.tsv. The report's variant section has a group row (module titles) and
# a column row; columns are named <module>__<column> (lower case, e.g. ccre_screen__classification).
# The join key is the report's "Original Input" chrom, pos, ref and alt, i.e. the candidate's own
# alleles (assets/awk/opencravat_input.awk). Dropped: the "Original Input" columns (the join key),
# the internal UID and, should a report hold them, back-lifted "hg19 coordinates" columns (the
# pipeline image has no hg19 annotator, so its reports have none). Empty cells (no annotation) are written as NA. A candidate
# variant without a report row is not_in_report. Annotation only: nothing is filtered or ranked.
function slug(s) { s = tolower(s); gsub(/[^a-z0-9]+/, "_", s); gsub(/^_+|_+$/, "", s); return s }
function module_of(t) { return (t in modname) ? modname[t] : slug(t) }
BEGIN {
    modname["Variant Annotation"] = "base"; modname["hg19 coordinates"] = "hg19"; modname["Original Input"] = "original_input"
    modname["Candidate cis-Regulatory Elements by ENCODE (SCREEN)"] = "ccre_screen"; modname["ClinVar"] = "clinvar"
    modname["ENCODE TFBS"] = "encode_tfbs"; modname["Repeat Sequences"] = "repeat"
}
ARGIND == 1 && FNR == 1 { for (i = 1; i <= NF; i++) cc[$i] = i; next }
ARGIND == 1 {
    key = $cc["chrom"] ":" $cc["pos"] ":" $cc["ref"] ":" $cc["alt"]
    if (!(key in seen)) {
        seen[key] = 1; order[++n] = key
        chrom[key] = $cc["chrom"]; pos[key] = $cc["pos"]; ref[key] = $cc["ref"]; alt[key] = $cc["alt"]; vtype[key] = $cc["vtype"]
    }
    if (!((key, $cc["element_id"]) in has_el)) { has_el[key, $cc["element_id"]] = 1; els[key] = (key in els) ? els[key] "," $cc["element_id"] : $cc["element_id"] }
    if (!((key, $cc["sample"]) in has_s)) { has_s[key, $cc["sample"]] = 1; ns[key]++ }
    next
}
ARGIND == 2 && /^#Report level:/ { section = $0; sub(/^#Report level: */, "", section); row = 0; next }
ARGIND == 2 && (/^#/ || NF == 0) { next }
ARGIND == 2 && section == "variant" {
    row++
    if (row == 1) {                                   # group row: a title starts each module's block
        for (i = 1; i <= NF; i++) { if ($i != "") g = module_of($i); grp[i] = g }
        nf = NF; next
    }
    if (row == 2) {                                   # column row
        if (NF > nf) nf = NF
        for (i = 1; i <= nf; i++) {
            m = (i in grp) ? grp[i] : g; name = m "__" slug($i)
            if (name in used) name = name "_" (++dup[name] + 1)
            used[name] = 1; cname[i] = name; cmod[i] = m
            if (m == "original_input") {
                t = slug($i)
                if (t == "chrom") ic = i; else if (t == "pos") ip = i; else if (t == "reference_allele") ir = i; else if (t == "alternate_allele") ia = i
            }
            if (m != "hg19" && m != "original_input" && name != "base__uid") keep[++nk] = i
        }
        if (!(ic && ip && ir && ia)) { print "OPENCRAVAT_TABLE: no Original Input chrom/pos/ref/alt columns in the report: cannot join it to the candidates" > "/dev/stderr"; failed = 1; exit 1 }
        header_seen = 1; next
    }
    k = $ic ":" $ip ":" $ir ":" $ia
    if (k in annot) next
    s = ""; for (j = 1; j <= nk; j++) { v = $keep[j]; s = s OFS (v == "" ? "NA" : v) }
    annot[k] = s; nrep++
}
END {
    if (failed) exit 1
    if (!header_seen) { print "OPENCRAVAT_TABLE: no variant-level section in the OpenCRAVAT text report" > "/dev/stderr"; exit 1 }
    h = "variant_id" OFS "chrom" OFS "pos" OFS "ref" OFS "alt" OFS "vtype" OFS "element_ids" OFS "n_samples" OFS "oc_status"
    for (j = 1; j <= nk; j++) h = h OFS cname[keep[j]]
    print h
    na = ""; for (j = 1; j <= nk; j++) na = na OFS "NA"
    for (i = 1; i <= n; i++) {
        k = order[i]
        print k, chrom[k], pos[k], ref[k], alt[k], vtype[k], els[k], ns[k], (k in annot ? "annotated" : "not_in_report") (k in annot ? annot[k] : na)
    }
}
