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
