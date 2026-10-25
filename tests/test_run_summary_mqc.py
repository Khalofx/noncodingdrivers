"""Tests for bin/run_summary_mqc.py.

The QC report is the one output a reader opens FIRST, so the properties under test
are about not misrepresenting the run:

  * significance is the combined q below --fdr, and the method is named
  * per-caller coverage reflects actual p-value presence, so a truncated caller is
    visible rather than averaged away
  * the report contains no candidate ranking and no driver language
"""
import gzip
import os
import subprocess
import sys

import pytest

_PATH = os.path.join(os.path.dirname(__file__), "..", "bin", "run_summary_mqc.py")

HDR = ["element_id", "chrom", "start", "end", "n_tested_callers", "vote_count",
       "consensus_q", "combination_method", "p_oncodrivefml", "p_activedriverwgs"]


def _consensus(tmp_path, rows, gz=False):
    p = tmp_path / ("c.tsv.gz" if gz else "c.tsv")
    op = gzip.open(p, "wt") if gz else open(p, "w")
    with op as fh:
        fh.write("\t".join(HDR) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(h, "")) for h in HDR) + "\n")
    return str(p)


def _run(tmp_path, cons, extra=None):
    out = tmp_path / "sections"
    cmd = [sys.executable, _PATH, "--consensus", cons, "--sample", "S1",
           "--outdir", str(out)] + (extra or [])
    proc = subprocess.run(cmd, capture_output=True, text=True)
    files = {}
    if out.exists():
        for f in sorted(os.listdir(out)):
            files[f] = (out / f).read_text()
    return proc, files


def _sec(files, key):
    for name, text in files.items():
        if key in name:
            return text
    return None


def _val(text, field):
    """Read a value from the `data:` block only.

    The `headers:` block declares the same keys (as nested config), so a naive
    first-match scan returns the header entry's empty value instead of the datum.
    """
    lines = text.split("\n")
    try:
        start = next(i for i, ln in enumerate(lines) if ln.rstrip() == "data:")
    except StopIteration:
        return None
    for line in lines[start + 1:]:
        s = line.strip()
        if s.startswith(field + ":"):
            return s.split(":", 1)[1].strip()
    return None


# ------------------------------------------------------- core counts

def test_overview_counts_match_the_table(tmp_path):
    rows = [
        {"element_id": "e1", "consensus_q": 0.001, "vote_count": 2, "n_tested_callers": 2},
        {"element_id": "e2", "consensus_q": 0.05,  "vote_count": 1, "n_tested_callers": 2},
        {"element_id": "e3", "consensus_q": 0.5,   "vote_count": 1, "n_tested_callers": 2},
        {"element_id": "e4", "consensus_q": 0.9,   "vote_count": 0, "n_tested_callers": 1},
    ]
    proc, files = _run(tmp_path, _consensus(tmp_path, rows))
    assert proc.returncode == 0, proc.stderr
    ov = _sec(files, "overview")
    assert _val(ov, "elements_tested") == "4"
    assert _val(ov, "significant") == "2"      # 0.001 and 0.05 (default --fdr 0.1)
    assert _val(ov, "multi_caller_agreement") == "1"
    assert _val(ov, "callers_run") == "2"


def test_gzipped_consensus_is_read(tmp_path):
    rows = [{"element_id": "e1", "consensus_q": 0.01, "vote_count": 1, "n_tested_callers": 1}]
    proc, files = _run(tmp_path, _consensus(tmp_path, rows, gz=True))
    assert proc.returncode == 0, proc.stderr
    assert _val(_sec(files, "overview"), "elements_tested") == "1"


# ------------------------------------------------------- threshold and method

def test_fdr_threshold_is_the_one_passed(tmp_path):
    """--fdr is the pipeline's --fdr_threshold; the count must follow it."""
    rows = [
        {"element_id": "e1", "consensus_q": 0.001, "vote_count": 2, "n_tested_callers": 2},
        {"element_id": "e2", "consensus_q": 0.05,  "vote_count": 1, "n_tested_callers": 2},
    ]
    proc, files = _run(tmp_path, _consensus(tmp_path, rows), extra=["--fdr", "0.01"])
    assert proc.returncode == 0, proc.stderr
    assert _val(_sec(files, "overview"), "significant") == "1"


def test_overview_names_the_combination_method(tmp_path):
    rows = [{"element_id": "e1", "consensus_q": 0.01, "combination_method": "ebm",
             "vote_count": 1, "n_tested_callers": 2}]
    proc, files = _run(tmp_path, _consensus(tmp_path, rows))
    assert proc.returncode == 0, proc.stderr
    assert "combined q (ebm)" in _sec(files, "overview")


# ------------------------------------------------------- caller coverage

def test_caller_coverage_counts_only_present_p_values(tmp_path):
    """A caller that scored fewer elements must show fewer -- not be averaged away."""
    rows = [
        {"element_id": "e1", "consensus_q": 0.5, "p_oncodrivefml": 0.2, "p_activedriverwgs": 0.3},
        {"element_id": "e2", "consensus_q": 0.6, "p_oncodrivefml": "NA", "p_activedriverwgs": 0.4},
        {"element_id": "e3", "consensus_q": 0.7, "p_oncodrivefml": "",   "p_activedriverwgs": 0.5},
    ]
    proc, files = _run(tmp_path, _consensus(tmp_path, rows))
    cov = _sec(files, "caller_coverage")
    assert _val(cov, "oncodrivefml") == "1"
    assert _val(cov, "activedriverwgs") == "3"


def test_absent_caller_file_is_omitted_with_a_message(tmp_path):
    rows = [{"element_id": "e1", "consensus_q": 0.01, "p_oncodrivefml": 0.001}]
    proc, files = _run(tmp_path, _consensus(tmp_path, rows),
                       extra=["--caller", "ghost=/nonexistent/x.tsv"])
    assert proc.returncode == 0, proc.stderr
    assert "absent, section omitted" in proc.stderr
    assert _sec(files, "caller_rows") is None, "fabricated a row count for a missing file"


def test_caller_rows_counts_data_lines_not_header(tmp_path):
    ct = tmp_path / "S1.oncodrivefml.tsv"
    ct.write_text("ELEMENT\tP_VALUE\na\t0.1\nb\t0.2\nc\t0.3\n")
    rows = [{"element_id": "e1", "consensus_q": 0.01, "p_oncodrivefml": 0.001}]
    proc, files = _run(tmp_path, _consensus(tmp_path, rows),
                       extra=["--caller", f"oncodrivefml={ct}"])
    assert proc.returncode == 0, proc.stderr
    assert _val(_sec(files, "caller_rows"), "oncodrivefml") == "3"


def test_malformed_caller_arg_is_a_hard_error(tmp_path):
    rows = [{"element_id": "e1", "consensus_q": 0.01}]
    proc, _ = _run(tmp_path, _consensus(tmp_path, rows),
                   extra=["--caller", "no_equals_sign"])
    assert proc.returncode != 0
    assert "expects name=path" in proc.stderr


# ------------------------------------------------------- scope discipline

def test_report_contains_no_candidate_ranking_or_driver_language(tmp_path):
    """A QC report answers 'did the run behave', not 'what did we find'.

    If element ids or driver language leak in, the report starts being read as a
    results summary -- which is exactly what the explanation layer is for.
    """
    rows = [
        {"element_id": "CTCF_binding_site_ENSR00000290858", "consensus_q": 1e-9,
         "vote_count": 2, "n_tested_callers": 2},
        {"element_id": "enhancer_ENSR00001267691", "consensus_q": 0.004,
         "vote_count": 1, "n_tested_callers": 2},
    ]
    proc, files = _run(tmp_path, _consensus(tmp_path, rows))
    assert proc.returncode == 0, proc.stderr
    blob = "\n".join(files.values())
    assert "ENSR" not in blob, "an element id leaked into the QC report"
    # Check the PROSE fields only. "driver" is an unavoidable substring elsewhere:
    # the pipeline is named noncodingdrivers, the section ids derive from that, and
    # `activedriverwgs` is a caller's actual name — none of which is a claim about
    # results. What matters is the language in section names and descriptions.
    prose = "\n".join(ln.split(":", 1)[1].lower()
                      for ln in blob.split("\n")
                      if ln.strip().startswith(("section_name:", "description:", "title:")))
    prose = prose.replace("noncodingdrivers", "")
    # "driver" appears legitimately in the DISCLAIMER ("a statistical threshold, not
    # a driver call"), which is the opposite of overclaiming — so this cannot be a
    # bare substring ban. Ban the ASSERTIVE phrasings instead, and separately require
    # the disclaimer to be present.
    for banned in ("top hit", "candidate ranking", "strongest", "most likely",
                   "is a driver", "are drivers", "driver identified", "driver found"):
        assert banned not in prose, f"results language in a QC report: {banned!r}"
    assert "not a driver call" in prose, \
        "the overview lost its disclaimer that significance is not a driver call"


def test_empty_consensus_is_a_hard_error_not_an_empty_report(tmp_path):
    """A header-only table means something upstream failed -- say so, loudly."""
    p = tmp_path / "empty.tsv"
    p.write_text("\t".join(HDR) + "\n")
    proc, _ = _run(tmp_path, str(p))
    assert proc.returncode != 0
    assert "no data rows" in proc.stderr


# ------------------------------------------------------- caller-name derivation
# Regression tests for a SILENT bug: two filename conventions reach this script and
# they have different field counts. Taking a fixed position (e.g. parts[-2]) yields
# "harmonized" for EVERY caller on the 4-field name, collapsing all callers into one
# map key and reporting a single row count under a meaningless label.
#
#   raw caller output   LIRI-JP.oncodrivefml.tsv              -> 3 fields
#   harmonised output   LIRI-JP.mutenricher.harmonized.tsv    -> 4 fields
#
# Both conventions are real: the caller modules emit `*.<caller>.tsv`, while
# modules/local/harmonizecalls emits `${meta.id}.${tool}.harmonized.tsv`.

import importlib.util as _ilu  # noqa: E402  (needs _PATH, defined above)

_SPEC = _ilu.spec_from_file_location("run_summary_mqc", _PATH)
rsm = _ilu.module_from_spec(_SPEC)
_SPEC.loader.exec_module(rsm)


@pytest.mark.parametrize("fname,expected", [
    # raw caller convention: <sample>.<caller>.tsv
    ("LIRI-JP.oncodrivefml.tsv", "oncodrivefml"),
    ("LIRI-JP.activedriverwgs.tsv", "activedriverwgs"),
    ("PACA-CA.mutenricher.tsv", "mutenricher"),
    # harmonised convention: <sample>.<caller>.harmonized.tsv  <- the bug case
    ("LIRI-JP.oncodrivefml.harmonized.tsv", "oncodrivefml"),
    ("LIRI-JP.mutenricher.harmonized.tsv", "mutenricher"),
    ("LIRI-JP.activedriverwgs.harmonised.tsv", "activedriverwgs"),
    # compressed variants
    ("LIRI-JP.oncodrivefml.tsv.gz", "oncodrivefml"),
    ("LIRI-JP.mutenricher.harmonized.tsv.gz", "mutenricher"),
])
def test_caller_name_derivation(fname, expected):
    assert rsm.caller_name_from_path(fname) == expected


def test_harmonised_names_do_not_all_collapse_to_harmonized():
    """The exact bug: three callers must give three distinct names, not 'harmonized'."""
    names = [rsm.caller_name_from_path(f) for f in [
        "LIRI-JP.oncodrivefml.harmonized.tsv",
        "LIRI-JP.activedriverwgs.harmonized.tsv",
        "LIRI-JP.mutenricher.harmonized.tsv",
    ]]
    assert names == ["oncodrivefml", "activedriverwgs", "mutenricher"]
    assert len(set(names)) == 3, "callers collapsed into one map key"
    assert "harmonized" not in names


def test_derivation_handles_a_full_path():
    assert rsm.caller_name_from_path(
        "/work/ab/cd12/LIRI-JP.mutenricher.harmonized.tsv") == "mutenricher"


def test_name_with_no_caller_field_is_visible_not_mislabelled(tmp_path):
    assert rsm.caller_name_from_path("weirdfile.tsv") == "weirdfile"


def test_caller_table_option_derives_names_end_to_end(tmp_path):
    """--caller-table takes paths only and labels each caller correctly."""
    for nm in ("oncodrivefml", "activedriverwgs", "mutenricher"):
        f = tmp_path / f"S1.{nm}.harmonized.tsv"
        f.write_text("element_id\tpvalue\n" + "".join(f"e{i}\t0.1\n" for i in range(nm.count("e") + 2)))
    rows = [{"element_id": "e1", "consensus_q": 0.01, "p_oncodrivefml": 0.001}]
    tables = [str(tmp_path / f"S1.{nm}.harmonized.tsv")
              for nm in ("oncodrivefml", "activedriverwgs", "mutenricher")]
    extra_args = []
    for t in tables:
        extra_args += ["--caller-table", t]
    proc, files = _run(tmp_path, _consensus(tmp_path, rows), extra=extra_args)
    assert proc.returncode == 0, proc.stderr
    sec = _sec(files, "caller_rows")
    assert sec is not None, "caller_rows section missing"
    for nm in ("oncodrivefml", "activedriverwgs", "mutenricher"):
        assert _val(sec, nm) is not None, f"{nm} missing from the report"
    assert _val(sec, "harmonized") is None, "callers were labelled 'harmonized'"


def test_colliding_caller_names_are_a_hard_error(tmp_path):
    """A convention change that merges two callers must FAIL, not silently merge."""
    a_ = tmp_path / "S1.samecaller.tsv"
    b_ = tmp_path / "other" / "S2.samecaller.tsv"
    b_.parent.mkdir()
    for f in (a_, b_):
        f.write_text("element_id\tpvalue\ne1\t0.1\n")
    rows = [{"element_id": "e1", "consensus_q": 0.01}]
    proc, _ = _run(tmp_path, _consensus(tmp_path, rows),
                   extra=["--caller-table", str(a_), "--caller-table", str(b_)])
    assert proc.returncode != 0
    assert "same name" in proc.stderr
