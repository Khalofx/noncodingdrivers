# Motifs for the MOTIFBREAKR module

`jaspar_core_vertebrates_2026.pfm` is the complete JASPAR CORE vertebrates non-redundant
collection of JASPAR 2026 (release 11): 1,019 position frequency matrices in JASPAR format
(`>MATRIX_ID NAME` followed by four `A/C/G/T [ counts ]` rows).

|          |                                                                                                            |
| -------- | ---------------------------------------------------------------------------------------------------------- |
| Source   | https://jaspar.elixir.no/download/data/2026/CORE/JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt |
| MD5      | `c63218acf88eb5fe10255b7c8ef23d48`                                                                         |
| Profiles | 1,019                                                                                                      |
| Licence  | CC BY 4.0 (JASPAR); cite the JASPAR 2026 paper (Nucleic Acids Res. 54(D1):D184, PMID 41325984)             |

The module reads the file with universalmotif (`read_jaspar`, then `convert_type(..., "PPM")`
with a pseudocount of 0.8 per base) and scores both alleles with motifbreakR. The collection
also holds non-human vertebrate profiles (lower-case names such as `Pax7`), so read TF names with
care. Use `--jaspar_motifs` to supply another JASPAR-format file.
