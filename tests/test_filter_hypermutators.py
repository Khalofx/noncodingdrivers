"""Tests for bin/filter_hypermutators.py (D1, D2: sample filters applied once before every caller).

Properties under test:
  * quartiles are the linear-interpolation quartiles (R type 7, numpy default)
  * the intogen rule needs BOTH conditions (count above the cap and above Q3 + k x IQR)
  * kept and removed rows partition the input, header and order preserved
  * a filter that would remove every sample stops the run
  * with donor information, one sample per donor (the first in alphabetical order) is kept,
    before the hypermutator rule; a sample given two donors stops the run
"""
import importlib.util
import os
import statistics
import subprocess
import sys

import pytest

_PATH = os.path.join(os.path.dirname(__file__), "..", "bin", "filter_hypermutators.py")
_spec = importlib.util.spec_from_file_location("filter_hypermutators", _PATH)
fh = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fh)

COUNTS = dict(zip([f"s{i:02d}" for i in range(1, 11)], [10, 12, 14, 16, 18, 20, 22, 24, 26, 120]))


def _table(tmp_path, counts):
    p = tmp_path / "in.tsv"
    with open(p, "w") as out:
        out.write("sample\tchrom\tpos\tref\talt\tvtype\n")
        pos = 0
        for s, n in counts.items():
            for _ in range(n):
                pos += 3
                out.write(f"{s}\tchr1\t{pos}\tC\tT\tSNP\n")
    return p


def _run(tmp_path, counts, *args):
    inp = _table(tmp_path, counts)
    prefix = str(tmp_path / "out")
    r = subprocess.run([sys.executable, _PATH, "--input", str(inp), "--prefix", prefix, *args], capture_output=True, text=True)
    return r, prefix, inp


def _rows(path):
    return open(path).read().splitlines()


def test_quantile_matches_statistics_inclusive():
    v = sorted(COUNTS.values())
    q = statistics.quantiles(v, n=4, method="inclusive")
    assert fh.quantile(v, 0.25) == pytest.approx(q[0])
    assert fh.quantile(v, 0.75) == pytest.approx(q[2])


def test_intogen_needs_both_conditions():
    removed, q1, q3, thr = fh.decide(COUNTS, "intogen", 15, 1.5)
    assert (q1, q3, thr) == (14.5, 23.5, 37.0)
    assert sorted(removed) == ["s10"]          # s04-s09 exceed 15 but not 37
    removed, *_ = fh.decide(COUNTS, "intogen", 200, 1.5)
    assert removed == {}                       # s10 is an outlier but below the count


def test_cap_and_none():
    assert sorted(fh.decide(COUNTS, "cap", 15, 1.5)[0]) == ["s04", "s05", "s06", "s07", "s08", "s09", "s10"]
    assert fh.decide(COUNTS, "none", 15, 1.5)[0] == {}


def test_outputs_partition_the_input(tmp_path):
    r, prefix, inp = _run(tmp_path, COUNTS, "--rule", "intogen", "--min-mutations", "15", "--cohort", "C1")
    assert r.returncode == 0, r.stderr
    src, kept, rem = _rows(inp), _rows(prefix + ".mutations.tsv"), _rows(prefix + ".removed_mutations.tsv")
    assert kept[0] == rem[0] == src[0]
    assert len(kept) - 1 == 162 and len(rem) - 1 == 120
    assert [row for row in src[1:] if not row.startswith("s10\t")] == kept[1:]
    sf = _rows(prefix + ".sample_filter.tsv")
    assert sf[0].split("\t")[:6] == ["sample", "donor", "n_mutations", "kept", "removed_by", "reason"]
    assert sf[1].split("\t")[:5] == ["s10", "", "120", "false", "hypermutator"]
    mqc = _rows(prefix + ".sample_filter_mqc.tsv")
    assert any(row.startswith("# plot_type: 'table'") for row in mqc)
    assert mqc[-1].split("\t")[:4] == ["C1", "10", "0", "1"]


def test_removing_every_sample_stops(tmp_path):
    r, _, _ = _run(tmp_path, {"a": 50, "b": 60}, "--rule", "cap", "--min-mutations", "10")
    assert r.returncode != 0 and "would remove all" in r.stderr


def _map(tmp_path, name, pairs):
    p = tmp_path / name
    p.write_text("sample\tdonor\n" + "".join(f"{s}\t{d}\n" for s, d in pairs))
    return str(p)


def _filter_rows(prefix):
    lines = _rows(prefix + ".sample_filter.tsv")
    cols = lines[0].split("\t")
    return {r.split("\t")[0]: dict(zip(cols, r.split("\t"))) for r in lines[1:]}


def test_donor_keeps_first_sample_in_alphabetical_order(tmp_path):
    dm = _map(tmp_path, "d.tsv", [("s03", "P1"), ("s01", "P1"), ("s05", "P2"), ("absent", "P3")])
    r, prefix, _ = _run(tmp_path, COUNTS, "--rule", "none", "--donor-map", dm)
    assert r.returncode == 0, r.stderr
    sf = _filter_rows(prefix)
    assert sf["s01"]["kept"] == "true" and sf["s01"]["donor"] == "P1"
    assert sf["s03"]["kept"] == "false" and sf["s03"]["removed_by"] == "donor" and "(s01)" in sf["s03"]["reason"]
    assert sf["s05"]["kept"] == "true" and sf["s02"]["donor"] == ""
    assert "absent" not in sf
    assert not any(row.startswith("s03\t") for row in _rows(prefix + ".mutations.tsv"))


def test_donor_step_runs_before_the_hypermutator_rule(tmp_path):
    dm = _map(tmp_path, "d.tsv", [("s09", "P9"), ("s10", "P9")])
    r, prefix, _ = _run(tmp_path, COUNTS, "--rule", "intogen", "--min-mutations", "15", "--donor-map", dm)
    assert r.returncode == 0, r.stderr
    sf = _filter_rows(prefix)
    assert sf["s10"]["removed_by"] == "donor"                       # s09 sorts first
    assert all(v["removed_by"] != "hypermutator" for v in sf.values())  # quartiles over the nine left
    assert sf["s01"]["iqr_threshold"] == "34.0"


def test_two_maps_are_merged_and_conflicts_stop(tmp_path):
    a = _map(tmp_path, "a.tsv", [("s01", "P1")])
    b = _map(tmp_path, "b.tsv", [("s02", "P1")])
    r, prefix, _ = _run(tmp_path, COUNTS, "--rule", "none", "--donor-map", a, "--donor-map", b)
    assert r.returncode == 0 and _filter_rows(prefix)["s02"]["removed_by"] == "donor"
    c = _map(tmp_path, "c.tsv", [("s01", "P2")])
    r, _, _ = _run(tmp_path, COUNTS, "--rule", "none", "--donor-map", a, "--donor-map", c)
    assert r.returncode != 0 and "two donors" in r.stderr
