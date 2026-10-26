# VEP_INPUT: candidate table -> Ensembl VEP default ("ensembl") input, one line per distinct
# variant: chrom, start, end, REF/ALT, +, id (id = chrom:pos:ref:alt as in the table).
# Port of v0.1 bin/vep_annotate.py to_vep_line: '-' alleles (MAF style) are native; an anchored
# VCF-style indel drops its shared first base. An insertion after `pos` has start = pos + 1, end = pos.
NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
{
    chrom = $col["chrom"]; p0 = $col["pos"]; r0 = $col["ref"]; a0 = $col["alt"]
    key = chrom ":" p0 ":" r0 ":" a0
    if (key in seen) next
    seen[key] = 1
    pos = p0 + 0; ref = toupper(r0); alt = toupper(a0)
    if (ref != "-" && alt != "-" && length(ref) != length(alt) && substr(ref, 1, 1) == substr(alt, 1, 1)) {
        ref = substr(ref, 2); alt = substr(alt, 2)
        if (ref == "") ref = "-"
        if (alt == "") alt = "-"
        if (ref != "-") pos++
    }
    if (ref == "-") { start = pos + 1; end = pos } else { start = pos; end = pos + length(ref) - 1 }
    print chrom, start, end, ref "/" alt, "+", key
}
