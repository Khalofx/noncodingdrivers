"""Tests for bin/filter_hypermutators.py (D2: one hypermutator filter before every caller).

Properties under test:
  * quartiles are the linear-interpolation quartiles (R type 7, numpy default)
  * the intogen rule needs BOTH conditions (count above the cap and above Q3 + k x IQR)
  * kept and removed rows partition the input, header and order preserved
  * a filter that would remove every sample stops the run
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
    assert sf[0].split("\t")[:4] == ["sample", "n_mutations", "kept", "reason"]
    assert sf[1].split("\t")[:3] == ["s10", "120", "false"]
    mqc = _rows(prefix + ".sample_filter_mqc.tsv")
    assert any(row.startswith("# plot_type: 'table'") for row in mqc)
    assert mqc[-1].split("\t")[:3] == ["C1", "10", "1"]


def test_removing_every_sample_stops(tmp_path):
    r, _, _ = _run(tmp_path, {"a": 50, "b": 60}, "--rule", "cap", "--min-mutations", "10")
    assert r.returncode != 0 and "would remove all" in r.stderr
