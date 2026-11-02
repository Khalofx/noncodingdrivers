"""bin/merge_oncodriveclustl_shards.py: Benjamini-Hochberg over all shards, OncodriveCLUSTL's order and ranks.

Fixtures (tests/data/oncodriveclustl): the 36 tested elements of the SYNTH-NC OncodriveCLUSTL 1.2.0
run (design prototype, one run, --seed 1) and its 39 clusters, as OncodriveCLUSTL wrote them
(expected.*), and the same rows split into two shuffled shards with shard-local q-values and
ranks overwritten (shard*.*). Merging the shards must give OncodriveCLUSTL's own tables back.
"""
import os
import subprocess
import sys

import pytest

pytest.importorskip("statsmodels")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, "tests", "data", "oncodriveclustl")
SCRIPT = os.path.join(ROOT, "bin", "merge_oncodriveclustl_shards.py")


def _merge(tmp_path, elements, clusters):
    oe, oc = tmp_path / "m.oncodriveclustl.tsv", tmp_path / "m.oncodriveclustl.clusters.tsv"
    proc = subprocess.run([sys.executable, SCRIPT, "--elements", *elements, "--clusters", *clusters,
                           "--out-elements", str(oe), "--out-clusters", str(oc)], capture_output=True, text=True)
    return proc, oe, oc


def test_two_shards_give_the_single_run_tables(tmp_path):
    proc, oe, oc = _merge(tmp_path, [os.path.join(D, "shard1.elements.tsv"), os.path.join(D, "shard2.elements.tsv")],
                          [os.path.join(D, "shard1.clusters.tsv"), os.path.join(D, "shard2.clusters.tsv")])
    assert proc.returncode == 0, proc.stderr
    assert oe.read_bytes() == open(os.path.join(D, "expected.elements.tsv"), "rb").read()
    assert oc.read_bytes() == open(os.path.join(D, "expected.clusters.tsv"), "rb").read()


def test_an_element_in_two_shards_is_an_error(tmp_path):
    e1 = os.path.join(D, "shard1.elements.tsv")
    c1 = os.path.join(D, "shard1.clusters.tsv")
    proc, _, _ = _merge(tmp_path, [e1, e1], [c1, c1])
    assert proc.returncode != 0 and "more than one shard" in proc.stderr


def test_shards_without_a_tested_element_keep_empty_q(tmp_path):
    hdr = open(os.path.join(D, "expected.elements.tsv")).readline()
    chdr = open(os.path.join(D, "expected.clusters.tsv")).readline()
    rows = ["eA\teA\tNon Available\t1\t+\t500\t2\t0\t0\t\t0\t\t\t\t\t\t\n",
            "eB\teB\tNon Available\t2\t+\t600\t3\t0\t0\t\t0\t\t\t\t\t\t\n"]
    paths = []
    for i, r in enumerate(rows):
        p = tmp_path / f"s{i}.tsv"
        p.write_text(hdr + r)
        c = tmp_path / f"c{i}.tsv"
        c.write_text(chdr)
        paths.append((str(p), str(c)))
    proc, oe, oc = _merge(tmp_path, [p for p, _ in paths], [c for _, c in paths])
    assert proc.returncode == 0, proc.stderr
    lines = oe.read_text().splitlines()
    assert len(lines) == 3 and all(line.split("\t")[14] == "" for line in lines[1:])
    assert oc.read_text() == chdr
