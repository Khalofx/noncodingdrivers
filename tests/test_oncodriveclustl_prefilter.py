"""bin/oncodriveclustl_prefilter.py: the rows of the ids with >= -emut SNVs, as bedtools intersect -c counts them."""
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(ROOT, "bin", "oncodriveclustl_prefilter.py")
D = os.path.join(ROOT, "tests", "data")


def _run(tmp_path, mutations, elements, min_snvs):
    out = tmp_path / "kept.bed"
    proc = subprocess.run([sys.executable, SCRIPT, "--mutations", str(mutations), "--elements", str(elements),
                           "--min-snvs", str(min_snvs), "--out", str(out)], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    return out.read_text().splitlines(), proc.stderr


def test_fixture_rows_are_the_prefiltered_fixture(tmp_path):
    kept, err = _run(tmp_path, os.path.join(D, "activedriverwgs", "mutations.tsv"),
                     os.path.join(D, "oncodriveclustl", "elements.bed"), 2)
    assert kept == open(os.path.join(D, "oncodriveclustl", "elements.prefiltered.bed")).read().splitlines()
    assert "4 of 9 element ids" in err


def test_bed_boundaries_vtype_recurrence_and_fragments(tmp_path):
    m = tmp_path / "m.tsv"
    m.write_text("sample\tchrom\tpos\tref\talt\tvtype\n"
                 "s1\tchr1\t100\tA\tC\tSNP\n"   # BED start 100: first base is 101, not in e1
                 "s1\tchr1\t101\tA\tC\tSNP\n"   # first base of e1
                 "s2\tchr1\t101\tA\tC\tSNP\n"   # recurrent: counted per sample
                 "s1\tchr1\t200\tA\tC\tSNP\n"   # last base of e1
                 "s1\tchr1\t201\tA\tC\tSNP\n"   # first base after e1
                 "s1\tchr1\t150\tA\tAC\tINS\n"  # not an SNV
                 "s1\tchr2\t11\tA\tC\tSNP\n"    # e2 fragment 1
                 "s1\tchr3\t11\tA\tC\tSNP\n")   # e2 fragment 2
    e = tmp_path / "e.bed"
    e.write_text("chr1\t100\t200\te1\t.\t+\n"
                 "chr1\t300\t400\te3\n"
                 "chr2\t10\t20\te2\n"
                 "chr3\t10\t20\te2\n"
                 "chr1\t200\t201\te4\n"  # the base 201 only
                 "chrUn_x\tNA\t5\te5\n")
    kept, _ = _run(tmp_path, m, e, 2)
    assert kept == ["chr1\t100\t200\te1", "chr2\t10\t20\te2", "chr3\t10\t20\te2"]
    kept, _ = _run(tmp_path, m, e, 4)
    assert kept == []
    kept, _ = _run(tmp_path, m, e, 3)
    assert kept == ["chr1\t100\t200\te1"]
    kept, _ = _run(tmp_path, m, e, 1)
    assert [r.split("\t")[3] for r in kept] == ["e1", "e2", "e2", "e4"]


def test_no_snv_keeps_no_row(tmp_path):
    m = tmp_path / "m.tsv"
    m.write_text("sample\tchrom\tpos\tref\talt\tvtype\n")
    e = tmp_path / "e.bed"
    e.write_text("chr1\t100\t200\te1\n")
    kept, err = _run(tmp_path, m, e, 2)
    assert kept == [] and "0 of 1 element ids" in err
