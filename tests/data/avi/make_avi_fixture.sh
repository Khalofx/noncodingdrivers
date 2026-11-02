#!/usr/bin/env bash
# Synthetic AVI score file over the SNVs of a candidate table, for the Karakoram integration run
# of --run_alphagenome_avi (docs/validation.md). Same layout as the AlphaGenome "AVI SNV scores"
# download (bgzip TSV: chrom, pos, ref, alt, raw, phred; tabix .tbi) but INVENTED, deterministic
# scores (phred = (pos mod 400) / 10, raw = phred / 10): not AlphaGenome output, licence-clean.
# usage: make_avi_fixture.sh <candidates.tsv> <out.tsv.gz>      (bgzip and tabix on PATH)
set -euo pipefail
cand="${1:?candidates.tsv}"
out="${2:?out.tsv.gz}"
{
    printf '#chrom\tpos\tref\talt\tavi_raw\tavi_phred\n'
    awk -F '\t' -v OFS='\t' '
        NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
        toupper($c["ref"]) ~ /^[ACGT]$/ && toupper($c["alt"]) ~ /^[ACGT]$/ {
            k = $c["chrom"] SUBSEP $c["pos"] SUBSEP $c["ref"] SUBSEP $c["alt"]
            if (k in seen) next
            seen[k] = 1
            ph = ($c["pos"] % 400) / 10
            print $c["chrom"], $c["pos"], toupper($c["ref"]), toupper($c["alt"]), sprintf("%.3f", ph / 10), sprintf("%.1f", ph)
        }' "$cand" | LC_ALL=C sort -t "$(printf '\t')" -k1,1 -k2,2n -k3,3 -k4,4
} | bgzip -c > "$out"
tabix -f -s 1 -b 2 -e 2 "$out"
