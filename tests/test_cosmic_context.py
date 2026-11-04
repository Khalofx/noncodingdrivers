"""bin/cosmic_context.py: exact GRCh38 matching of candidate mutations against a COSMIC NCV file.

Every COSMIC-format row here is synthetic (invented identifiers such as COSVTEST1 and sample
names such as SYN_OTHER); no COSMIC record is used anywhere in the tests.
"""
import csv
import gzip
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import pytest

BIN = Path(__file__).resolve().parents[1] / "bin" / "cosmic_context.py"

# Column order of the v104 file (MUTATION_NC_ID, SAMPLE_NAME before COSMIC_SAMPLE_ID)
COSMIC_HEADER = [
    "GENE_SYMBOL", "COSMIC_GENE_ID", "TRANSCRIPT_ACCESSION", "MUTATION_NC_ID",
    "SAMPLE_NAME", "COSMIC_SAMPLE_ID", "COSMIC_PHENOTYPE_ID",
    "GENOMIC_MUTATION_ID", "LEGACY_MUTATION_ID", "ZYGOSITY", "CHROMOSOME",
    "GENOME_START", "GENOME_STOP", "GENOMIC_WT_ALLELE", "GENOMIC_MUT_ALLELE",
    "COSMIC_STUDY_ID", "PUBMED_PMID", "HGVSG", "MUTATION_SOMATIC_STATUS",
]
CANDIDATE_HEADER = ["sample", "chrom", "pos", "ref", "alt", "vtype", "element_id", "fragment_start", "fragment_end"]
NAME = "Cosmic_NonCodingVariants_v999_GRCh38.tsv.gz"


def cosmic_row(gene, transcript, nc_id, sample_name, coss, coso, cosv, chrom, start, stop, ref, alt,
               study, pmid, status, legacy="", hgvsg=""):
    return [gene, "COSGTEST1", transcript, nc_id, sample_name, coss, coso, cosv, legacy, "het", chrom,
            start, stop, ref, alt, study, pmid, hgvsg, status]


COSMIC_ROWS = [
    cosmic_row("G1", "ENSTTEST1", "COSNTEST1", "S1", "COSSTEST1", "COSOTEST1", "COSVTEST1", "1", 101, 101,
               "C", "G", "COSUTEST1", "111", "Confirmed somatic variant", hgvsg="1:g.101C>G"),
    # same genomic observation on a second transcript: one observation
    cosmic_row("G2", "ENSTTEST2", "COSNTEST2", "S1", "COSSTEST1", "COSOTEST1", "COSVTEST1", "1", 101, 101,
               "C", "G", "COSUTEST1", "111", "Confirmed somatic variant", hgvsg="1:g.101C>G"),
    # MAF '-' must match COSMIC's empty mutant allele
    cosmic_row("", "", "COSNTEST3", "SYN_OTHER", "COSSTEST2", "COSOTEST2", "COSVTEST2", "1", 110, 110,
               "A", "", "", "222", "Variant of unknown origin", hgvsg="1:g.110del"),
    # same position, different allele: not an exact match
    cosmic_row("", "", "COSNTEST4", "SYN_OTHER", "COSSTEST3", "COSOTEST2", "COSVTEST3", "1", 105, 105,
               "G", "T", "", "333", "Confirmed somatic variant"),
    # a donor-prefixed alias of S2 is NOT an exact name: stays not_assessed
    cosmic_row("", "", "COSNTEST5", "S2_T01", "COSSTEST4", "COSOTEST3", "COSVTEST4", "1", 107, 107,
               "T", "C", "COSUTEST2", "444", "Confirmed somatic variant"),
    # COSMIC_SAMPLE_ID equal to a cohort sample id: possible_self_overlap
    cosmic_row("", "", "COSNTEST6", "SYN_X", "S3", "COSOTEST3", "COSVTEST5", "2", 50, 50,
               "A", "G", "COSUTEST3", "", "Confirmed somatic variant"),
    # MT is chrM
    cosmic_row("", "", "COSNTEST7", "SYN_Y", "COSSTEST5", "COSOTEST4", "COSVTEST6", "MT", 300, 300,
               "A", "G", "", "", "Variant of unknown origin"),
    # unrelated row
    cosmic_row("", "", "COSNTEST8", "SYN_Z", "COSSTEST6", "COSOTEST4", "COSVTEST7", "X", 999, 999,
               "A", "G", "", "", "Confirmed somatic variant"),
]

CANDIDATE_ROWS = [
    ["S1", "chr1", 101, "C", "G", "SNP", "E1", 100, 110],
    ["S2", "chr1", 105, "G", "A", "SNP", "E1", 100, 110],
    ["S2", "chr1", 105, "G", "A", "SNP", "E2", 104, 108],   # one mutation in two elements
    ["S3", "chr1", 110, "A", "-", "DEL", "E1", 100, 110],
    ["S2", "chr1", 107, "T", "C", "SNP", "E2", 104, 108],
    ["S3", "chr2", 50, "A", "G", "SNP", "E3", 40, 60],
    ["S4", "chrM", 300, "A", "G", "SNP", "E4", 290, 310],
    ["S5", "chr2", 55, "C", "T", "SNP", "E3", 40, 60],       # absent from COSMIC
]


def write_tsv(path, header, rows, gz=False, crlf=False):
    terminator = "\r\n" if crlf else "\n"
    opener = gzip.open if gz else open
    with opener(path, "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator=terminator)
        writer.writerow(header)
        writer.writerows(rows)
    return path


def md5(path):
    return hashlib.md5(Path(path).read_bytes()).hexdigest()


def run(tmp_path, name=NAME, release="v999", header=COSMIC_HEADER, rows=COSMIC_ROWS, crlf=True,
        expected_md5=None, consensus=True, check=True):
    cosmic = write_tsv(tmp_path / name, header, rows, gz=name.endswith(".gz"), crlf=crlf)
    candidates = write_tsv(tmp_path / "TEST.candidates.tsv", CANDIDATE_HEADER, CANDIDATE_ROWS)
    cmd = [sys.executable, str(BIN), "--cosmic", str(cosmic), "--release", release,
           "--candidates", str(candidates), "--cohort", "TEST", "--outdir", str(tmp_path / "out"),
           "--progress-every", "0"]
    if consensus:
        cons = write_tsv(tmp_path / "TEST.consensus.tsv",
                         ["element_id", "chrom", "start", "end", "consensus_q", "consensus_rank", "is_consensus"],
                         [["E2", "chr1", 104, 108, 0.01, 1, "True"], ["E1", "chr1", 100, 110, 0.02, 2, "True"],
                          ["E3", "chr2", 40, 60, 0.03, 3, "True"], ["E4", "chrM", 290, 310, 0.04, 4, "True"],
                          ["E5", "chr3", 1, 50, 0.05, 5, "True"], ["NS", "chr4", 1, 50, 0.9, 6, "False"]])
        cmd += ["--consensus", str(cons)]
    if expected_md5 is not None:
        cmd += ["--expected-md5", expected_md5]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if check:
        assert result.returncode == 0, result.stderr
    return result, tmp_path / "out", cosmic


def read_tsv(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def by_pos(rows):
    return {(row["chrom"], int(row["pos"])): row for row in rows}


# --- ported from benchmarks/tests/test_cosmic_compare.py (candidate-table input) -------------

def test_candidate_rows_define_the_distinct_variants_and_their_elements(tmp_path):
    _, out, _ = run(tmp_path)
    rows = by_pos(read_tsv(out / "TEST.cosmic.mutation_summary.tsv"))
    assert set(rows) == {("chr1", 101), ("chr1", 105), ("chr1", 107), ("chr1", 110), ("chr2", 50),
                         ("chr2", 55), ("chrM", 300)}
    assert rows[("chr1", 105)]["element_ids"] == "E1,E2" and rows[("chr1", 105)]["n_elements"] == "2"


def test_exact_matching_and_unmatched_are_explicit(tmp_path):
    _, out, _ = run(tmp_path)
    rows = by_pos(read_tsv(out / "TEST.cosmic.mutation_summary.tsv"))
    assert rows[("chr1", 101)]["cosmic_match"] == "yes"
    assert rows[("chr1", 101)]["n_confirmed_somatic_samples"] == "1"
    assert rows[("chr1", 105)]["cosmic_match"] == "no"            # wrong ALT in COSMIC
    assert rows[("chr1", 105)]["independence_status"] == "no_match"
    assert rows[("chr1", 110)]["cosmic_match"] == "yes"           # '-' == empty COSMIC allele
    assert rows[("chrM", 300)]["cosmic_match"] == "yes"           # COSMIC MT == chrM
    assert rows[("chr2", 55)]["cosmic_match"] == "no"
    exact = by_pos(read_tsv(out / "TEST.cosmic.exact_matches.tsv"))
    assert exact[("chr1", 110)]["somatic_status"] == "Variant of unknown origin"


def test_transcript_rows_are_collapsed_but_annotations_preserved(tmp_path):
    _, out, _ = run(tmp_path)
    exact = read_tsv(out / "TEST.cosmic.exact_matches.tsv")
    at_101 = [row for row in exact if row["pos"] == "101"]
    assert len(at_101) == 1
    assert at_101[0]["source_transcript_rows"] == "2"
    assert at_101[0]["gene_symbols"] == "G1,G2"
    assert at_101[0]["mutation_nc_ids"] == "COSNTEST1,COSNTEST2"


def test_provenance_totals_distinguish_rows_observations_and_variants(tmp_path):
    _, out, cosmic = run(tmp_path)
    p = json.loads((out / "TEST.cosmic.provenance.json").read_text())
    assert p["cosmic_release"] == "v999" and p["cosmic_assembly"] == "GRCh38" and p["cosmic_file"] == NAME
    assert p["candidate_elements"] == 5                    # E1-E4 with candidates + E5 from the consensus
    assert p["candidate_sample_events"] == 7               # (sample, variant, vtype) events
    assert p["candidate_distinct_mutations"] == 7
    assert p["cosmic_rows_scanned"] == 8
    assert p["cosmic_raw_exact_rows"] == 6
    assert p["cosmic_collapsed_observations"] == 5
    assert p["cosmic_matched_distinct_mutations"] == 5
    assert p["cosmic_unmatched_distinct_mutations"] == 2
    assert p["cosmic_file_md5"] == md5(cosmic) and p["md5_check"] == "not_requested"


def test_exact_sample_name_overlap_only(tmp_path):
    _, out, _ = run(tmp_path)
    rows = by_pos(read_tsv(out / "TEST.cosmic.mutation_summary.tsv"))
    # SAMPLE_NAME == cohort sample (S1)
    assert rows[("chr1", 101)]["independence_status"] == "possible_self_overlap"
    assert rows[("chr1", 101)]["n_exact_sample_name_overlaps"] == "1"
    # COSMIC_SAMPLE_ID == cohort sample (S3)
    assert rows[("chr2", 50)]["independence_status"] == "possible_self_overlap"
    exact = by_pos(read_tsv(out / "TEST.cosmic.exact_matches.tsv"))
    assert exact[("chr1", 101)]["exact_sample_name_overlap"] == "yes"


def test_aliases_are_not_matched_and_stay_not_assessed(tmp_path):
    _, out, _ = run(tmp_path)
    rows = by_pos(read_tsv(out / "TEST.cosmic.mutation_summary.tsv"))
    assert rows[("chr1", 107)]["cosmic_match"] == "yes"
    assert rows[("chr1", 107)]["independence_status"] == "not_assessed"     # S2_T01 is not S2
    assert rows[("chr1", 110)]["independence_status"] == "not_assessed"


def test_element_summary_counts_in_consensus_order(tmp_path):
    _, out, _ = run(tmp_path)
    rows = read_tsv(out / "TEST.cosmic.element_summary.tsv")
    assert [row["element_id"] for row in rows] == ["E2", "E1", "E3", "E4", "E5"]
    e1 = rows[1]
    assert (e1["n_candidate_sample_events"], e1["n_distinct_candidate_mutations"], e1["n_cosmic_matched_mutations"]) == ("3", "3", "2")
    assert e1["n_exact_sample_name_overlaps"] == "1" and e1["consensus_q"] == "0.02"
    assert rows[4]["n_candidate_sample_events"] == "0" and rows[4]["cosmic_match"] == "no"


def test_element_summary_without_consensus_uses_candidate_elements(tmp_path):
    _, out, _ = run(tmp_path, consensus=False)
    rows = read_tsv(out / "TEST.cosmic.element_summary.tsv")
    assert [row["element_id"] for row in rows] == ["E1", "E2", "E3", "E4"]
    assert rows[0]["consensus_rank"] == "NA"


# --- licence layout: count tables carry no COSMIC identifiers -----------------------------------

def test_count_tables_hold_no_cosmic_identifiers_and_notice_is_written(tmp_path):
    _, out, _ = run(tmp_path)
    for name in ("TEST.cosmic.mutation_summary.tsv", "TEST.cosmic.element_summary.tsv", "TEST.cosmic.provenance.json"):
        text = (out / name).read_text()
        assert "TEST" in text
        for token in ("COSV", "COSS", "COSN", "COSU", "COSO", "SYN_OTHER", "111"):
            assert token + "TEST" not in text and "SYN_OTHER" not in text, (name, token)
    assert "COSVTEST1" in (out / "TEST.cosmic.exact_matches.tsv").read_text()
    notice = (out / "NOTICE.txt").read_text()
    assert "COSMIC v999" in notice and "Do not redistribute" in notice


# --- new checks of the pipeline module ----------------------------------------------------------

def test_crlf_and_lf_files_give_identical_results(tmp_path):
    a = tmp_path / "crlf"
    b = tmp_path / "lf"
    a.mkdir()
    b.mkdir()
    _, out_a, _ = run(a, crlf=True)
    _, out_b, _ = run(b, crlf=False)
    for name in ("mutation_summary", "element_summary", "exact_matches"):
        assert (out_a / f"TEST.cosmic.{name}.tsv").read_text() == (out_b / f"TEST.cosmic.{name}.tsv").read_text()
    assert json.loads((out_a / "TEST.cosmic.provenance.json").read_text())["cosmic_line_endings"] == "CRLF"
    assert json.loads((out_b / "TEST.cosmic.provenance.json").read_text())["cosmic_line_endings"] == "LF"


def test_cosmic_ncv_id_is_accepted_as_the_id_column(tmp_path):
    header = [("COSMIC_NCV_ID" if h == "MUTATION_NC_ID" else h) for h in COSMIC_HEADER]
    _, out, _ = run(tmp_path, header=header)
    p = json.loads((out / "TEST.cosmic.provenance.json").read_text())
    assert p["cosmic_mutation_id_column"] == "COSMIC_NCV_ID" and p["cosmic_matched_distinct_mutations"] == 5


def test_missing_id_column_is_rejected(tmp_path):
    keep = [i for i, h in enumerate(COSMIC_HEADER) if h != "MUTATION_NC_ID"]
    header = [COSMIC_HEADER[i] for i in keep]
    rows = [[r[i] for i in keep] for r in COSMIC_ROWS]
    result, _, _ = run(tmp_path, header=header, rows=rows, check=False)
    assert result.returncode != 0 and "mutation id column" in result.stderr


def test_columns_are_read_by_name_not_position(tmp_path):
    order = list(reversed(range(len(COSMIC_HEADER))))
    header = [COSMIC_HEADER[i] for i in order]
    rows = [[r[i] for i in order] for r in COSMIC_ROWS]
    _, out, _ = run(tmp_path, header=header, rows=rows)
    assert json.loads((out / "TEST.cosmic.provenance.json").read_text())["cosmic_matched_distinct_mutations"] == 5


@pytest.mark.parametrize("name, message", [
    ("Cosmic_NonCodingVariants_v999_GRCh37.tsv.gz", "GRCh37"),
    ("Cosmic_NonCodingVariants_v998_GRCh38.tsv.gz", "release v998"),
    ("cosmic_ncv.tsv.gz", "not COSMIC's own"),
])
def test_file_name_release_and_assembly_are_checked(tmp_path, name, message):
    result, out, _ = run(tmp_path, name=name, check=False)
    assert result.returncode != 0 and message in result.stderr
    assert not out.exists() or not any(out.iterdir())


def test_release_parameter_format_is_checked(tmp_path):
    result, _, _ = run(tmp_path, release="104", check=False)
    assert result.returncode != 0 and "v<number>" in result.stderr


def test_md5_match_is_recorded_and_mismatch_stops(tmp_path):
    good = tmp_path / "good"
    bad = tmp_path / "bad"
    good.mkdir()
    bad.mkdir()
    cosmic = write_tsv(good / "probe.tsv.gz", COSMIC_HEADER, COSMIC_ROWS, gz=True, crlf=True)
    # gzip output is deterministic only with a fixed mtime; write once, then reuse the bytes
    data = cosmic.read_bytes()
    (good / NAME).write_bytes(data)
    expected = hashlib.md5(data).hexdigest()
    cand = write_tsv(good / "TEST.candidates.tsv", CANDIDATE_HEADER, CANDIDATE_ROWS)
    base = [sys.executable, str(BIN), "--cosmic", str(good / NAME), "--release", "v999",
            "--candidates", str(cand), "--cohort", "TEST", "--progress-every", "0"]
    ok = subprocess.run(base + ["--expected-md5", expected, "--outdir", str(good / "out")], capture_output=True, text=True)
    assert ok.returncode == 0, ok.stderr
    p = json.loads((good / "out" / "TEST.cosmic.provenance.json").read_text())
    assert p["md5_check"] == "match" and p["cosmic_file_md5"] == expected == p["cosmic_file_md5_expected"]
    wrong = subprocess.run(base + ["--expected-md5", "0" * 32, "--outdir", str(bad / "out")], capture_output=True, text=True)
    assert wrong.returncode != 0 and "expected 00000000" in wrong.stderr
    assert not (bad / "out").exists()


def test_plain_tsv_md5_covers_the_whole_file(tmp_path):
    name = NAME[:-3]
    _, out, cosmic = run(tmp_path, name=name)
    assert json.loads((out / "TEST.cosmic.provenance.json").read_text())["cosmic_file_md5"] == md5(cosmic)


@pytest.mark.parametrize("value", ["chr1", "M", "GL000220.1", ""])
def test_non_cosmic_chromosome_names_are_rejected(tmp_path, value):
    rows = [list(r) for r in COSMIC_ROWS]
    rows[-1][COSMIC_HEADER.index("CHROMOSOME")] = value
    result, _, _ = run(tmp_path, rows=rows, check=False)
    assert result.returncode != 0 and "not a COSMIC GRCh38 chromosome" in result.stderr


def test_outputs_are_deterministic(tmp_path):
    a = tmp_path / "a"
    b = tmp_path / "b"
    a.mkdir()
    b.mkdir()
    _, out_a, _ = run(a)
    _, out_b, _ = run(b)
    for f in sorted(p.name for p in out_a.iterdir()):
        if f.endswith(".json"):
            ja, jb = (json.loads((o / f).read_text()) for o in (out_a, out_b))
            ja.pop("cosmic_file_md5")
            jb.pop("cosmic_file_md5")     # gzip header mtime differs between the two fixtures
            assert ja == jb
        else:
            assert (out_a / f).read_text() == (out_b / f).read_text(), f
