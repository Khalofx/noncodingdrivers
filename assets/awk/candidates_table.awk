# CANDIDATE_TABLE: bedtools intersect -wa -wb output (mutation BED6+3, then fragment BED4)
# -> patient-resolved candidate table, one row per sample x mutation x significant element.
BEGIN { print "sample", "chrom", "pos", "ref", "alt", "vtype", "element_id", "fragment_start", "fragment_end" }
{ print $4, $1, $3, $7, $8, $9, $13, $11, $12 }
