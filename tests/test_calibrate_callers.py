"""bin/calibrate_callers.py (CALIBRATECALLERS, D32): inflation factor and band status."""
import os
import subprocess
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(ROOT, "bin", "calibrate_callers.py")


def _bh(p):
    p = np.asarray(p)
    n = p.size
    order = np.argsort(p)
    q = np.empty(n)
    q[order] = np.minimum.accumulate((p[order] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1)


def _tables(d, name, rng, n=4000, power_a=1.0, power_b=1.0):
    ids = [f"e{i}" for i in range(n)]
    pa = rng.uniform(size=n) ** power_a          # ActiveDriverWGS-like column: power > 1 shrinks p (inflation)
    pb = rng.uniform(size=n) ** power_b          # OncodriveFML-like column: power < 1 enlarges p (deflation)
    pe = rng.uniform(size=n)
    cons = pd.DataFrame({"element_id": ids, "p_activedriverwgs": pa, "q_activedriverwgs": _bh(pa),
                         "p_oncodrivefml": pb, "q_oncodrivefml": _bh(pb), "ebm_p": pe, "ebm_q": _bh(pe),
                         "acat_p": pe, "acat_q": _bh(pe)})
    cons.to_csv(d / f"{name}.consensus.tsv", sep="\t", index=False)
    exp = np.where(np.arange(n) < n // 2, 2.0, 0.5)        # half the elements have expected >= 1
    pd.DataFrame({"id": ids, "element_muts_exp": exp}).to_csv(d / f"{name}.activedriverwgs.tsv", sep="\t", index=False)


def test_band_status(tmp_path):
    rng = np.random.default_rng(3)
    _tables(tmp_path, "C", rng)
    for r in (1, 2, 3):
        _tables(tmp_path, f"C.null{r}", rng, power_a=1.6, power_b=0.6)
    out, qq = tmp_path / "cal.tsv", tmp_path / "qq.tsv"
    subprocess.run([sys.executable, SCRIPT, "--cohort", "C", "--observed-consensus", str(tmp_path / "C.consensus.tsv"),
                    "--observed-adwgs", str(tmp_path / "C.activedriverwgs.tsv"),
                    "--null-consensus", *[str(tmp_path / f"C.null{r}.consensus.tsv") for r in (3, 1, 2)],
                    "--null-adwgs", *[str(tmp_path / f"C.null{r}.activedriverwgs.tsv") for r in (1, 2, 3)],
                    "--band", "0.8,1.2", "--calibration", str(out), "--qq", str(qq)], check=True)
    cal = pd.read_csv(out, sep="\t", dtype=str)
    obs = cal[(cal.run == "observed")].set_index("method")
    assert set(obs.index) == {"activedriverwgs", "oncodrivefml", "ebm", "acat"}
    assert (obs["n_elements"] == "2000").all()                               # expected >= 1 only
    assert 0.85 < float(obs.loc["activedriverwgs", "lambda"]) < 1.15
    assert list(cal[cal.run.str.startswith("null") & (cal.run != "null_summary")].run.unique()) == ["null1", "null2", "null3"]
    summ = cal[cal.run == "null_summary"].set_index("method")
    assert summ.loc["activedriverwgs", "status"] == "excluded_for_inflation"   # p^1.6 inflates lambda
    assert summ.loc["oncodrivefml", "status"] == "conservative"                # p^0.6 deflates it
    assert summ.loc["ebm", "status"] == "within_band"
    assert float(summ.loc["ebm", "lambda_min"]) <= float(summ.loc["ebm", "lambda"]) <= float(summ.loc["ebm", "lambda_max"])
    q = pd.read_csv(qq, sep="\t")
    assert set(q.run) == {"observed", "null1", "null2", "null3"} and len(q) > 0


def _add_clustl(d, name, rng, n_tested, n_low, n_q, n_analysed):
    # OncodriveCLUSTL columns: p only for n_tested elements (its clustered ones), n_low of them at p <= 0.05,
    # n_q of those at q < 0.1; the native table lists n_analysed elements (tested ones first)
    cons = pd.read_csv(d / f"{name}.consensus.tsv", sep="\t")
    p = np.full(len(cons), np.nan)
    q = np.full(len(cons), np.nan)
    p[:n_tested] = rng.uniform(0.06, 1.0, size=n_tested)
    q[:n_tested] = np.maximum(p[:n_tested], 0.5)
    p[:n_low] = np.linspace(0.001, 0.05, n_low) if n_low else p[:0]
    q[:n_low] = 0.3
    q[:n_q] = 0.05
    cons.insert(3, "p_oncodriveclustl", p)
    cons.insert(4, "q_oncodriveclustl", q)
    (d / "cl").mkdir(exist_ok=True)
    cons.to_csv(d / "cl" / f"{name}.consensus.tsv", sep="\t", index=False, na_rep="NA")
    pd.DataFrame({"SYMBOL": cons.element_id[:n_analysed], "ENSID": cons.element_id[:n_analysed]}).to_csv(
        d / f"{name}.oncodriveclustl.tsv", sep="\t", index=False)


def _calibrate(d, out, qq, clustl=False, fpr=None):
    c = d / "cl" if clustl else d                                         # consensus tables with or without the columns
    cmd = [sys.executable, SCRIPT, "--cohort", "C", "--observed-consensus", str(c / "C.consensus.tsv"),
           "--observed-adwgs", str(d / "C.activedriverwgs.tsv"),
           "--null-consensus", *[str(c / f"C.null{r}.consensus.tsv") for r in (2, 1, 3)],
           "--null-adwgs", *[str(d / f"C.null{r}.activedriverwgs.tsv") for r in (1, 2, 3)],
           "--band", "0.8,1.2", "--calibration", str(out), "--qq", str(qq)]
    if clustl:
        cmd += ["--observed-oncodriveclustl", str(d / "C.oncodriveclustl.tsv"),
                "--null-oncodriveclustl", *[str(d / f"C.null{r}.oncodriveclustl.tsv") for r in (3, 2, 1)]]
    if fpr:
        cmd += ["--fpr", str(fpr)]
    subprocess.run(cmd, check=True)


def test_oncodriveclustl_fpr_check_replaces_lambda(tmp_path):
    rng = np.random.default_rng(7)
    _tables(tmp_path, "C", rng)
    for r in (1, 2, 3):
        _tables(tmp_path, f"C.null{r}", rng)
    # SYNTH-NC-like counts (docs/validation.md): tested 36 / 30, 30, 34; at p <= 0.05: 3 / 3, 0, 0; q < 0.1: 2 / 0, 0, 0
    for name, tested, low, nq, analysed in (("C", 36, 3, 2, 3088), ("C.null1", 30, 3, 0, 3195), ("C.null2", 30, 0, 0, 3227), ("C.null3", 34, 0, 0, 3219)):
        _add_clustl(tmp_path, name, rng, tested, low, nq, min(analysed, 4000))
    _calibrate(tmp_path, tmp_path / "plain.tsv", tmp_path / "plain_qq.tsv", fpr=tmp_path / "plain_fpr.tsv")
    _calibrate(tmp_path, tmp_path / "cl.tsv", tmp_path / "cl_qq.tsv", clustl=True, fpr=tmp_path / "fpr.tsv")
    # no lambda and no QQ points for OncodriveCLUSTL; every other method's rows are unchanged, byte for byte
    assert open(tmp_path / "cl.tsv").read() == open(tmp_path / "plain.tsv").read()
    assert open(tmp_path / "cl_qq.tsv").read() == open(tmp_path / "plain_qq.tsv").read()
    assert not (tmp_path / "plain_fpr.tsv").exists()                       # written only with OncodriveCLUSTL
    f = pd.read_csv(tmp_path / "fpr.tsv", sep="\t", dtype=str).set_index("run")
    assert list(f.index) == ["observed", "null1", "null2", "null3", "null_summary"]
    assert (f.loc["observed", "n_tested"], f.loc["observed", "n_p_le_0.05"], f.loc["observed", "n_q_lt_fdr"]) == ("36", "3", "2")
    assert [f.loc[f"null{r}", "n_analysed"] for r in (1, 2, 3)] == ["3195", "3227", "3219"]
    s = f.loc["null_summary"]
    assert (s["replicate"], s["n_analysed"], s["n_tested"], s["n_p_le_0.05"], s["n_q_lt_fdr"], s["status"]) == ("3", "9641", "94", "3", "0", "pass")
    assert abs(float(s["fpr_p_le_0.05"]) - 3 / 94) < 1e-6
    assert float(s["min_p"]) == 0.001


def test_oncodriveclustl_fpr_fails_with_two_null_discoveries(tmp_path):
    rng = np.random.default_rng(11)
    _tables(tmp_path, "C", rng)
    for r in (1, 2, 3):
        _tables(tmp_path, f"C.null{r}", rng)
    for name, tested, low, nq in (("C", 10, 1, 1), ("C.null1", 10, 2, 2), ("C.null2", 10, 0, 0), ("C.null3", 10, 1, 1)):
        _add_clustl(tmp_path, name, rng, tested, low, nq, 50)
    _calibrate(tmp_path, tmp_path / "cl.tsv", tmp_path / "cl_qq.tsv", clustl=True, fpr=tmp_path / "fpr.tsv")
    s = pd.read_csv(tmp_path / "fpr.tsv", sep="\t", dtype=str).set_index("run").loc["null_summary"]
    assert (s["n_q_lt_fdr"], s["status"]) == ("2", "fail")                # the largest per-replicate count
