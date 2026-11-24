"""Unit tests of benchmarks/score_synth_nc.py on a hand-made answer key and consensus table."""
import json
import os
import subprocess
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(os.path.dirname(HERE), "score_synth_nc.py")
sys.path.insert(0, os.path.dirname(HERE))
from score_synth_nc import compare, score  # noqa: E402

TRUTH = {"donors": 60, "planted": [
    {"planted_id": "E1", "element_id": "promCore__TERT__X", "chrom": "chr5", "start0": 100, "end": 400,
     "donors_mutated": 18, "design": "TERT promoter hotspots"},
    {"planted_id": "E2", "element_id": "enhancer_A", "chrom": "chr1", "start0": 1000, "end": 1800,
     "donors_mutated": 3, "design": "3 donors, random positions"},
    {"planted_id": "E6", "element_id": "enhancer_D", "chrom": "chr2", "start0": 5000, "end": 5600,
     "donors_mutated": 12, "design": "12 donors, A/T only (spectrum decoy)"},
]}


def _cons():
    return pd.DataFrame([
        ("promCore__TERT__X", "chr5", 100, 400, 1e-8, 1e-9),
        ("5utr__TERT__X", "chr5", 350, 600, 0.01, 0.02),     # non-planted, overlapping E1
        ("enhancer_A", "chr1", 1000, 1800, 0.3, 0.05),       # recovered by ACAT only
        ("enhancer_D", "chr2", 5000, 5600, 0.04, 0.2),       # decoy significant under EBM only
        ("enhancer_Z", "chr9", 1, 500, 0.02, 0.5),           # non-planted hit under EBM
    ], columns=["element_id", "chrom", "start", "end", "ebm_q", "acat_q"])


def test_score_counts():
    r = score(TRUTH, _cons(), 0.1)
    assert (r["planted_drivers"], r["planted_decoys"]) == (2, 1)
    assert r["ebm"] == {"recovered": 1, "decoy_significant": True, "non_planted_significant": 2,
                        "overlapping_planted_significant": 1, "missed": ["E2 enhancer_A (3 donors)"]}
    assert r["acat"] == {"recovered": 2, "decoy_significant": False, "non_planted_significant": 1,
                         "overlapping_planted_significant": 1, "missed": []}


def test_compare_and_cli_exit_status(tmp_path):
    r = score(TRUTH, _cons(), 0.1)
    ref = {"ebm": {"recovered": 1, "decoy_significant": True, "non_planted_significant": 2},
           "acat": {"recovered": 2, "decoy_significant": False, "non_planted_significant": 1}}
    assert compare(r, ref) == []
    assert compare(r, {"ebm": {"recovered": 2}}) == ["ebm.recovered: expected 2, got 1"]
    t, c, ok, bad = tmp_path / "t.json", tmp_path / "c.tsv", tmp_path / "ok.json", tmp_path / "bad.json"
    json.dump(TRUTH, open(t, "w"))
    _cons().to_csv(c, sep="\t", index=False)
    json.dump({"synth_nc": ref}, open(ok, "w"))
    json.dump({"synth_nc": {"ebm": {"recovered": 2}}}, open(bad, "w"))

    def run(ref_path):
        return subprocess.run([sys.executable, SCRIPT, "--truth", str(t), "--consensus", str(c),
                               "--reference", str(ref_path)], capture_output=True, text=True)

    assert run(ok).returncode == 0
    p = run(bad)
    assert p.returncode == 1 and "FAIL ebm.recovered" in p.stdout
