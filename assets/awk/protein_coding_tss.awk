# PC_TSS: gene TSS BED6+2 (assets/awk/gene_tss.awk) -> the protein_coding rows only, in the same
# (sorted) order, for the nearest protein-coding gene of CANCERGENECONTEXT.
$8 == "protein_coding"
