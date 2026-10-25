import importlib.util
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "rheinbay_concordance", ROOT / "rheinbay_concordance.py")
MOD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MOD)


class FakeLifter:
    def convert_coordinate(self, chrom, pos):
        return [(chrom, pos + 1000, "+", 1)]


def test_bed12_blocks_and_lift_preserve_half_open_intervals():
    fields = ["chr1", "100", "250", "element", "0", "+", "100", "250", "0",
              "2", "10,20,", "0,100,"]
    blocks = MOD.bed12_blocks(fields)
    assert blocks == [("chr1", 100, 110), ("chr1", 200, 220)]
    lifted, status = MOD.lift_blocks(blocks, FakeLifter())
    assert status == "lifted"
    assert lifted == [("chr1", 1100, 1110), ("chr1", 1200, 1220)]


def test_overlap_uses_bed_half_open_boundaries():
    assert MOD.overlaps(("chr1", 10, 20), ("chr1", 19, 21))
    assert not MOD.overlaps(("chr1", 10, 20), ("chr1", 20, 21))
    assert not MOD.overlaps(("chr1", 10, 20), ("chr2", 10, 20))


def test_published_status_separates_filtering_from_significance():
    assert MOD.published_status({"Post-filter q-value": "0.01",
                                 "Final filter judgement": "1"}) == "retained_significant"
    assert MOD.published_status({"Post-filter q-value": "0.5",
                                 "Final filter judgement": "1"}) == "retained_not_significant"
    assert MOD.published_status({"Post-filter q-value": "0.01",
                                 "Final filter judgement": "0"}) == "filtered"


def test_read_table5_without_openpyxl(tmp_path):
    shared = """<?xml version='1.0' encoding='UTF-8'?>
    <sst xmlns='http://schemas.openxmlformats.org/spreadsheetml/2006/main'>
      <si><t>ID</t></si><si><t>tissue</t></si><si><t>gene</t></si>
      <si><t>elt1</t></si><si><t>Liver-HCC</t></si><si><t>TERT</t></si>
      <si><t>Other</t></si>
    </sst>"""
    sheet = """<?xml version='1.0' encoding='UTF-8'?>
    <worksheet xmlns='http://schemas.openxmlformats.org/spreadsheetml/2006/main'><sheetData>
      <row r='1'><c r='A1' t='s'><v>0</v></c><c r='C1' t='s'><v>1</v></c><c r='D1' t='s'><v>2</v></c></row>
      <row r='2'><c r='A2' t='s'><v>3</v></c><c r='C2' t='s'><v>4</v></c><c r='D2' t='s'><v>5</v></c></row>
      <row r='3'><c r='A3' t='s'><v>3</v></c><c r='C3' t='s'><v>6</v></c><c r='D3' t='s'><v>5</v></c></row>
    </sheetData></worksheet>"""
    book = tmp_path / "table.xlsx"
    with zipfile.ZipFile(book, "w") as archive:
        archive.writestr("xl/sharedStrings.xml", shared)
        archive.writestr("xl/worksheets/sheet5.xml", sheet)
    assert MOD.read_table5(book, "Liver-HCC") == [
        {"ID": "elt1", "tissue": "Liver-HCC", "gene": "TERT"}
    ]


def test_compare_matches_locus_not_identifier():
    table = [{"ID": "published", "target": "promCore", "tissue": "Liver-HCC",
              "gene": "TERT", "Post-filter q-value": "1e-8",
              "Final filter judgement": "1"}]
    elements = {"published": {"chrom": "chr5", "start": 100, "end": 110,
                               "blocks": [("chr5", 100, 110)]}}
    consensus = [{"element_id": "different_hg38_id", "chrom": "chr5", "start": "1105",
                  "end": "1120", "consensus_rank": "5", "acat_q": "1e-5",
                  "n_snv": "11", "is_consensus": "True"}]
    candidates = [{"element_id": "different_hg38_id", "sample": "P1", "chrom": "chr5",
                   "pos": "1107", "ref": "G", "alt": "A", "spectrum_status": "clean",
                   "motif_status": "scored", "nearest_gene": "GENE", "nearest_gene_status": "assigned"}]
    rows = MOD.compare(table, elements, FakeLifter(), consensus, candidates)
    assert rows[0]["pipeline_recovered"] == "true"
    assert rows[0]["pipeline_element_ids"] == "different_hg38_id"
    assert rows[0]["pipeline_intervals"] == "chr5:1105-1120"
    assert rows[0]["pipeline_distinct_patients"] == "1"
