# AVI_TABLE: candidate table + Ensembl VEP tab output with the AVI plugin columns -> one row per
# distinct candidate variant with its AlphaGenome Variant Impact (AVI) scores.
# Inputs (in this order): <cohort>.candidates.tsv, <cohort>.candidates.vep.tab (gunzipped by GAWK).
# The VEP Uploaded_variation is the candidate key chrom:pos:ref:alt (assets/awk/vep_input.awk).
# The plugin scores SNVs only (one base A/C/G/T on both sides); every other variant is
# not_scored_non_snv. An SNV without a score in the AVI file is not_in_avi_file. Annotation and
# hypothesis only: no threshold is applied and nothing is filtered.
# The run fails when the VEP header has no AVI_PHRED column: the plugin did not load (VEP itself
# only warns and exits 0 without its plugin folder or the .tbi).
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
ARGIND == 2 && /^##/ { next }
ARGIND == 2 && /^#/ {
    sub(/^#/, ""); for (i = 1; i <= NF; i++) vc[$i] = i
    header_seen = 1
    if (!("AVI_PHRED" in vc)) {
        print "AVI_TABLE: the VEP table has no AVI_PHRED column: the AVI plugin did not load (check --alphagenome_avi, its .tbi and --dir_plugins)" > "/dev/stderr"
        failed = 1; exit 1
    }
    next
}
ARGIND == 2 {
    v = $vc["Uploaded_variation"]; in_vep[v] = 1
    p = $vc["AVI_PHRED"]; r = ("AVI_RAW" in vc) ? $vc["AVI_RAW"] : "-"
    if (p != "-" && p != "" && !(v in phred)) { phred[v] = p; raw[v] = (r == "-" || r == "") ? "NA" : r }
}
END {
    if (failed) exit 1
    if (!header_seen) { print "AVI_TABLE: no VEP header line (#Uploaded_variation ...) in the VEP table" > "/dev/stderr"; exit 1 }
    print "variant_id", "chrom", "pos", "ref", "alt", "vtype", "element_ids", "n_samples", "avi_phred", "avi_raw", "avi_status"
    for (i = 1; i <= n; i++) {
        k = order[i]
        snv = (toupper(ref[k]) ~ /^[ACGT]$/ && toupper(alt[k]) ~ /^[ACGT]$/)
        if (k in phred) { status = "scored"; ph = phred[k]; rw = raw[k] }
        else { ph = "NA"; rw = "NA"; status = !snv ? "not_scored_non_snv" : ((k in in_vep) ? "not_in_avi_file" : "not_in_vep_output") }
        print k, chrom[k], pos[k], ref[k], alt[k], vtype[k], els[k], ns[k], ph, rw, status
    }
}
