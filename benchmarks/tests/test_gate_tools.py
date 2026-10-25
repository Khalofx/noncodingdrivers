import json
import os
import subprocess
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "benchmarks"))

import compare_to_baseline  # noqa: E402
import table_to_vcfs  # noqa: E402

BASELINE = os.path.join(ROOT, "tests", "regression", "baseline.json")


def _v1_table(path, tert_rank_shift=0):
    base = json.load(open(BASELINE))
    tert = base["liver_hcc_public"]["tert"]
    rows = []
    for i in range(1, 41):
        rows.append({"element_id": f"el{i}", "p_activedriverwgs": 1e-3 * i, "p_oncodrivefml": 0.5,
                     "ebm_p": 1e-8 * i, "ebm_q": 0.001 if i <= 30 else 0.5,
                     "acat_p": 1e-8 * i, "acat_q": 0.001 if i <= 27 else 0.5})
    rows.insert(4, {"element_id": tert["element"], "p_activedriverwgs": tert["adwgs_p"],
                    "p_oncodrivefml": tert["ofml_p"], "ebm_p": 4.5e-8 + tert_rank_shift,
                    "ebm_q": base["combination_check_2026_10_06"]["tert_liver_hcc"]["ebm_q"],
                    "acat_p": 4.5e-8 + tert_rank_shift, "acat_q": tert["acat_q"]})
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def test_compare_passes_on_matching_table(tmp_path):
    t = tmp_path / "v1.tsv"
    _v1_table(t)
    out = tmp_path / "report.tsv"
    rc = compare_to_baseline.main(["--baseline", BASELINE, "--cohort", "liver_hcc_public",
                                   "--v1", str(t), "--out", str(out)])
    rep = pd.read_csv(out, sep="\t")
    assert rc == 0, rep.to_string()
    assert set(rep["check"]) >= {"n_significant_ebm", "tert_rank_ebm", "tert_ebm_q", "tert_p_activedriverwgs"}


def test_compare_fails_on_rank_change(tmp_path):
    t = tmp_path / "v1.tsv"
    _v1_table(t, tert_rank_shift=1.0)
    rc = compare_to_baseline.main(["--baseline", BASELINE, "--cohort", "liver_hcc_public",
                                   "--v1", str(t), "--out", str(tmp_path / "r.tsv")])
    assert rc == 1


def test_compare_frozen_set(tmp_path):
    t = tmp_path / "v1.tsv"
    _v1_table(t)
    fr = pd.read_csv(t, sep="\t").rename(columns={"ebm_q": "consensus_q"})
    fr.to_csv(tmp_path / "frozen.tsv", sep="\t", index=False)
    rc = compare_to_baseline.main(["--baseline", BASELINE, "--cohort", "liver_hcc_public", "--v1", str(t),
                                   "--frozen", str(tmp_path / "frozen.tsv"), "--out", str(tmp_path / "r.tsv")])
    rep = pd.read_csv(tmp_path / "r.tsv", sep="\t")
    assert rc == 0, rep.to_string()
    assert "frozen_significant_set_equal" in set(rep["check"])


def _combination_csvs(tmp_path):
    """Comparison and discordant CSVs consistent with _v1_table: EBM 31 (30 + TERT), ACAT 28 (27 + TERT)."""
    pd.DataFrame([{"cohort": "Liver-HCC-public", "elements": 41, "shared_ebm": 28, "acat_only": 0, "ebm_only": 3}]) \
        .to_csv(tmp_path / "comparison.csv", index=False)
    pd.DataFrame([{"cohort": "Liver-HCC-public", "status": "EBM_only", "element_id": f"el{i}"} for i in (28, 29, 30)]) \
        .to_csv(tmp_path / "discordant.csv", index=False)


def test_compare_combination_csvs_pass(tmp_path):
    t = tmp_path / "v1.tsv"
    _v1_table(t)
    d = pd.read_csv(t, sep="\t")
    d.loc[d["element_id"].isin(["el28", "el29", "el30"]), "p_oncodrivefml"] = 0.01   # EBM-only: both callers < 0.05
    d.to_csv(t, sep="\t", index=False)
    _combination_csvs(tmp_path)
    rc = compare_to_baseline.main(["--baseline", BASELINE, "--cohort", "liver_hcc_public", "--v1", str(t),
                                   "--comparison", str(tmp_path / "comparison.csv"),
                                   "--discordant", str(tmp_path / "discordant.csv"), "--out", str(tmp_path / "r.tsv")])
    rep = pd.read_csv(tmp_path / "r.tsv", sep="\t")
    assert rc == 0, rep.to_string()
    assert {"n_elements_tested", "discordant_ebm_only_set_equal", "ebm_only_all_both_callers_p_lt_0.05"} <= set(rep["check"])


def test_compare_discordant_pattern_fails(tmp_path):
    """An EBM-only element with OncodriveFML p >= 0.05 breaks the release pattern."""
    t = tmp_path / "v1.tsv"
    _v1_table(t)                       # el28-el30 keep p_oncodrivefml = 0.5
    _combination_csvs(tmp_path)
    rc = compare_to_baseline.main(["--baseline", BASELINE, "--cohort", "liver_hcc_public", "--v1", str(t),
                                   "--discordant", str(tmp_path / "discordant.csv"), "--out", str(tmp_path / "r.tsv")])
    rep = pd.read_csv(tmp_path / "r.tsv", sep="\t").set_index("check")
    assert rc == 1
    assert rep.at["ebm_only_all_both_callers_p_lt_0.05", "status"] == "FAIL"


def test_compare_acat_q_tolerance(tmp_path):
    # The package ACAT reproduces the release in-house ACAT q within --acat-q-rtol, not exactly
    for factor, expected_rc in ((1 - 3.1e-7, 0), (1.001, 1)):
        t = tmp_path / f"v1_{factor}.tsv"
        _v1_table(t)
        d = pd.read_csv(t, sep="\t")
        tert = json.load(open(BASELINE))["liver_hcc_public"]["tert"]
        d.loc[d["element_id"] == tert["element"], "acat_q"] = tert["acat_q"] * factor
        d.to_csv(t, sep="\t", index=False)
        rc = compare_to_baseline.main(["--baseline", BASELINE, "--cohort", "liver_hcc_public", "--method", "acat",
                                       "--v1", str(t), "--out", str(tmp_path / f"r_{factor}.tsv")])
        assert rc == expected_rc, pd.read_csv(tmp_path / f"r_{factor}.tsv", sep="\t").to_string()


def test_compare_acat_q_uses_recorded_cohort_bound(tmp_path):
    # A 4% package-vs-in-house ACAT difference passes only within the bound recorded for the cohort
    t = tmp_path / "v1.tsv"
    _v1_table(t)
    d = pd.read_csv(t, sep="\t")
    tert = json.load(open(BASELINE))["liver_hcc_public"]["tert"]
    d.loc[d["element_id"] == tert["element"], "acat_q"] = tert["acat_q"] / 1.04
    d.loc[d["element_id"].isin(["el28", "el29", "el30"]), "p_oncodrivefml"] = 0.01
    d.to_csv(t, sep="\t", index=False)
    for dl, expected_rc in ((0.0181, 0), (0.001, 1)):
        pd.DataFrame([{"cohort": "Liver-HCC-public", "elements": 41, "shared_ebm": 28, "acat_only": 0, "ebm_only": 3,
                       "acat_pkg_vs_custom_max_dlog10": dl}]).to_csv(tmp_path / "comparison.csv", index=False)
        rc = compare_to_baseline.main(["--baseline", BASELINE, "--cohort", "liver_hcc_public", "--method", "acat",
                                       "--v1", str(t), "--comparison", str(tmp_path / "comparison.csv"),
                                       "--out", str(tmp_path / "r.tsv")])
        assert rc == expected_rc, pd.read_csv(tmp_path / "r.tsv", sep="\t").to_string()


def test_table_to_vcfs(tmp_path):
    src = os.path.join(ROOT, "tests", "data", "cohort_mutations.tsv")
    fai = os.path.join(ROOT, "tests", "data", "tiny.fa.fai")
    out = tmp_path / "vcfs"
    d = pd.read_csv(src, sep="\t")
    fai_arg = []
    contigs = {line.split("\t")[0] for line in open(fai)}
    if set(d["chrom"]) <= contigs:
        fai_arg = ["--fai", fai]
    assert table_to_vcfs.main([src, str(out)] + fai_arg) == 0
    files = sorted(os.listdir(out))
    assert len(files) == d["sample"].nunique()
    body = [line for f in files for line in open(out / f) if not line.startswith("#")]
    assert len(body) == len(d.drop_duplicates(["sample", "chrom", "pos", "ref", "alt"]))
    assert all(line.split("\t")[6] == "PASS" for line in body)


def test_cli_help():
    for script in ("compare_to_baseline.py", "table_to_vcfs.py"):
        r = subprocess.run([sys.executable, os.path.join(ROOT, "benchmarks", script), "--help"],
                           capture_output=True, text=True)
        assert r.returncode == 0
