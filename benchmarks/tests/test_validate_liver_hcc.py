"""Tests for benchmarks/validate_liver_hcc.py (core layer)."""

import csv
import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("validate_liver_hcc", ROOT / "validate_liver_hcc.py")
MOD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MOD)

COLS = ["element_id", "chrom", "start", "end", "consensus_q", "consensus_rank", "is_consensus", "combination_method"]


EXP = MOD.EXPECTED["acat"]


def _results(tmp_path, n_sig=EXP["significant"], tert_rank=EXP["tert_rank"], tert_q=EXP["tert_q"], method="acat"):
    out = tmp_path / "results" / "consensus"
    out.mkdir(parents=True)
    rows = []
    # EXP["loci"] loci: n_pairs pairs of overlapping elements on chr1, the others single (TERT on chr5 is one locus)
    n_pairs = EXP["significant"] - EXP["loci"]
    pos = 1000
    for i in range(n_sig - 1):
        overlap = i < 2 * n_pairs and i % 2 == 1
        start = pos - 50 if overlap else pos + 10_000
        pos = start + 100
        rows.append({"element_id": f"el{i}", "chrom": "chr1", "start": start, "end": pos,
                     "consensus_q": 0.01, "consensus_rank": i + 10, "is_consensus": "True", "combination_method": method})
    rows.append({"element_id": MOD.TERT_ELEMENT, "chrom": "chr5", "start": 1294666, "end": 1295601,
                 "consensus_q": tert_q, "consensus_rank": tert_rank, "is_consensus": "True", "combination_method": method})
    rows.append({"element_id": "elx", "chrom": "chr2", "start": 1, "end": 50, "consensus_q": 0.5,
                 "consensus_rank": 900, "is_consensus": "False", "combination_method": method})
    with open(out / "Liver-HCC-public.consensus.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    return tmp_path / "results"


def test_core_passes_on_release_values(tmp_path):
    res = _results(tmp_path)
    assert MOD.main(["--results", str(res), "--method", "acat"]) == 0


def test_core_fails_on_tert_rank(tmp_path):
    res = _results(tmp_path, tert_rank=EXP["tert_rank"] + 1)
    with pytest.raises(SystemExit, match="TERT rank"):
        MOD.main(["--results", str(res), "--method", "acat"])


def test_core_fails_on_wrong_method(tmp_path):
    res = _results(tmp_path, method="ebm")
    with pytest.raises(SystemExit, match="combined with"):
        MOD.main(["--results", str(res), "--method", "acat"])


def test_locus_count_merges_overlaps():
    rows = [{"chrom": "chr1", "start": 0, "end": 10}, {"chrom": "chr1", "start": 5, "end": 20},
            {"chrom": "chr1", "start": 20, "end": 30}, {"chrom": "chr2", "start": 0, "end": 10}]
    assert MOD.locus_count(rows) == 3
