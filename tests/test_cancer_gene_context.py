"""bin/cancer_gene_context.py: gene-list labels of the nearest gene, the nearest protein-coding
gene and the ENCODE-rE2G targets. Gene lists here are synthetic (invented tiers and roles)."""
import csv
import gzip
import subprocess
import sys
from pathlib import Path

import pytest

BIN = Path(__file__).resolve().parents[1] / "bin" / "cancer_gene_context.py"

NEAREST = [
    ["element_id", "chrom", "start", "end", "nearest_gene", "nearest_gene_ensembl_id", "gene_type", "nearest_tss",
     "gene_strand", "distance_to_tss", "abs_distance_to_tss", "relative_position", "nearest_gene_status"],
    ["TERT_promoter_fixture", "chr5", "1294800", "1295300", "ENSG00000300381", "ENSG00000300381.1", "lncRNA",
     "1295139", "+", "-89", "89", "overlapping_tss", "assigned"],
    ["frag_F", "chr5", "1400000", "1400600", "GENEX", "ENSG00000900001.3", "protein_coding", "1414230", "-", "-13930",
     "13930", "upstream", "assigned"],
    ["lonely", "chr21", "100", "200", "", "", "", "", "", "", "", "", "no_same_chromosome_gene"],
]
# bedtools closest -d -t first: midpoint BED6+2, protein-coding TSS BED6+2, distance
NEAREST_PC = [
    ["chr5", "1295049", "1295050", "TERT_promoter_fixture", "0", ".", "1294800", "1295300",
     "chr5", "1295067", "1295068", "ENSG00000164362.22", "0", "-", "TERT", "protein_coding", "18"],
    ["chr5", "1400299", "1400300", "frag_F", "0", ".", "1400000", "1400600",
     "chr5", "1414229", "1414230", "ENSG00000900001.3", "0", "-", "GENEX", "protein_coding", "13930"],
    ["chr21", "149", "150", "lonely", "0", ".", "100", "200", ".", "-1", "-1", ".", "-1", ".", ".", ".", "-1"],
]
LINKS = [
    ["element_id", "tissue", "biosample", "accession", "re2g_element", "re2g_class", "target_gene",
     "target_gene_ensembl_id", "target_gene_tss", "is_self_promoter", "re2g_score", "overlap_bp"],
    ["TERT_promoter_fixture", "liver", "HepG2", "ENCFFTEST01", "chr5:1295000-1295500", "promoter", "TERT",
     "ENSG00000164362", "1295068", "TRUE", "1", "300"],
    ["frag_F", "liver", "HepG2", "ENCFFTEST01", "chr5:1400100-1400400", "intergenic", "GENEY", "ENSG00000900002",
     "1500000", "FALSE", "0.4", "300"],
    ["frag_F", "liver", "HepG2", "ENCFFTEST01", "chr5:1400100-1400400", "intergenic", "TERT", "ENSG00000164362",
     "1295068", "FALSE", "0.3", "300"],
]
CENSUS = [
    ["GENE_SYMBOL", "NAME", "COSMIC_GENE_ID", "CHROMOSOME", "TIER", "ROLE_IN_CANCER", "SYNONYMS"],
    ["TERT", "telomerase reverse transcriptase", "COSGTEST1", "5", "1", "oncogene", "TRT,ENSG00000164362.22,TP2"],
    ["GENEX", "gene x", "COSGTEST2", "5", "2", "TSG", "GX"],
    ["GENEY", "gene y", "COSGTEST3", "5", "1", "oncogene, fusion", ""],
    ["GENEZ", "gene z", "COSGTEST4", "1", "2", "TSG", "ENSG00000900009"],
    ["GENEW", "gene w", "COSGTEST5", "2", "1", "oncogene", ""],
]


def write(path, rows, delimiter="\t", gz=False, crlf=False):
    opener = gzip.open if gz else open
    with opener(path, "wt", newline="") as handle:
        csv.writer(handle, delimiter=delimiter, lineterminator="\r\n" if crlf else "\n").writerows(rows)
    return path


def run(tmp_path, gene_list, links=True, pc=True, name="test list"):
    nearest = write(tmp_path / "t.nearest_gene.tsv", NEAREST)
    cmd = [sys.executable, str(BIN), "--nearest", str(nearest), "--gene-list", str(gene_list),
           "--list-name", name, "--output", str(tmp_path / "out.tsv")]
    if pc:
        cmd += ["--nearest-pc", str(write(tmp_path / "t.nearest_pc.bed", NEAREST_PC))]
    if links:
        cmd += ["--links", str(write(tmp_path / "t.links.tsv", LINKS))]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    with open(tmp_path / "out.tsv") as handle:
        return {row["element_id"]: row for row in csv.DictReader(handle, delimiter="\t")}


def test_census_layout_labels_nearest_pc_gene_and_re2g_targets(tmp_path):
    rows = run(tmp_path, write(tmp_path / "census.tsv.gz", CENSUS, gz=True, crlf=True))
    tert = rows["TERT_promoter_fixture"]
    assert tert["nearest_in_list"] == "no"                              # the lncRNA is not in the list
    assert tert["nearest_protein_coding_gene"] == "TERT"
    assert (tert["nearest_pc_in_list"], tert["nearest_pc_match_by"]) == ("yes", "ensembl_id")   # ENSG in SYNONYMS
    assert (tert["nearest_pc_list_tier"], tert["nearest_pc_list_role"]) == ("1", "oncogene")
    assert tert["re2g_targets_in_list"] == "TERT" and tert["re2g_distal_targets_in_list"] == ""
    f = rows["frag_F"]
    assert (f["nearest_in_list"], f["nearest_match_by"], f["nearest_list_role"]) == ("yes", "symbol", "TSG")
    assert f["re2g_targets_in_list"] == "GENEY,TERT" and f["re2g_distal_targets_in_list"] == "GENEY,TERT"
    assert f["list_genes"] == "5" and f["list_name"] == "test list" and len(f["list_md5"]) == 32


def test_plain_list_with_header_and_ensembl_ids(tmp_path):
    plain = write(tmp_path / "plain.tsv", [["gene_symbol", "ensembl_gene_id"], ["GENEX", ""],
                                          ["NOT_TERT", "ENSG00000164362"], ["GENEQ", ""]])
    rows = run(tmp_path, plain)
    # matched by Ensembl id although the symbol differs
    assert (rows["TERT_promoter_fixture"]["nearest_pc_in_list"], rows["TERT_promoter_fixture"]["nearest_pc_match_by"]) == ("yes", "ensembl_id")
    assert rows["frag_F"]["nearest_in_list"] == "yes" and rows["frag_F"]["nearest_list_tier"] == ""


def test_plain_list_without_header(tmp_path):
    rows = run(tmp_path, write(tmp_path / "genes.txt", [["TERT"], ["GENEY"], ["genex"]]))
    assert rows["frag_F"]["nearest_in_list"] == "no"                    # symbols are case-sensitive
    assert rows["TERT_promoter_fixture"]["nearest_pc_in_list"] == "yes"
    assert rows["frag_F"]["re2g_targets_in_list"] == "GENEY,TERT"


def test_symbol_match_is_rejected_when_the_list_gives_another_ensembl_id(tmp_path):
    rows = run(tmp_path, write(tmp_path / "l.tsv", [["gene_symbol", "ensembl_gene_id"], ["GENEX", "ENSG00000000001"]]))
    assert rows["frag_F"]["nearest_in_list"] == "no"


def test_without_rE2G_or_protein_coding_inputs_columns_are_NA(tmp_path):
    rows = run(tmp_path, write(tmp_path / "g.txt", [["TERT"]]), links=False, pc=False)
    assert rows["frag_F"]["re2g_targets_in_list"] == "NA"
    assert rows["frag_F"]["nearest_protein_coding_gene"] == "NA"


def test_unassigned_elements_are_NA(tmp_path):
    rows = run(tmp_path, write(tmp_path / "g.txt", [["TERT"]]))
    assert rows["lonely"]["nearest_in_list"] == "NA" and rows["lonely"]["nearest_pc_in_list"] == "NA"
    assert rows["lonely"]["nearest_protein_coding_gene"] == ""


def test_csv_list_and_bad_list(tmp_path):
    rows = run(tmp_path, write(tmp_path / "g.csv", [["gene_symbol", "tier"], ["GENEX", "2"]], delimiter=","))
    assert rows["frag_F"]["nearest_list_tier"] == "2"
    bad = write(tmp_path / "bad.tsv", [["symbol", "x"], ["A", "B"]])
    result = subprocess.run([sys.executable, str(BIN), "--nearest", str(write(tmp_path / "n.tsv", NEAREST)),
                             "--gene-list", str(bad), "--output", str(tmp_path / "o.tsv")], capture_output=True, text=True)
    assert result.returncode != 0 and "gene_symbol" in result.stderr


@pytest.mark.parametrize("content", [b"", b"# only a comment\n"])
def test_empty_list_is_rejected(tmp_path, content):
    path = tmp_path / "empty.txt"
    path.write_bytes(content)
    result = subprocess.run([sys.executable, str(BIN), "--nearest", str(write(tmp_path / "n.tsv", NEAREST)),
                             "--gene-list", str(path), "--output", str(tmp_path / "o.tsv")], capture_output=True, text=True)
    assert result.returncode != 0 and "empty" in result.stderr
