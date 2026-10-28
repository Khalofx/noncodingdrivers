#!/usr/bin/env python3
"""CALIBRATECALLERS: per-caller calibration on neutral (SIMULATENULL) cohorts (in-house, D32).

For the observed run and every null replicate, and for every method (each combined caller, the
EBM and ACAT combinations):
  lambda = chi2_1 quantile of the median p / 0.454936 (genomic inflation factor), over the
           elements with an ActiveDriverWGS expected count element_muts_exp >= --min-expected
           in that run (all elements with a p-value when no ActiveDriverWGS table is given);
  n_q_lt_fdr = elements with q < --fdr (all elements).
QQ points: expected and observed -log10 p over the same element set, all points with observed
-log10 p >= 2 plus up to 500 evenly spaced ranks.
Null summary per method: lambda over the replicates against --band. A method with a replicate
above the band is excluded_for_inflation; one with a replicate below it (and none above) is
conservative; otherwise within_band. Reported only: no p- or q-value is changed.
"""
import argparse
import math
import re
import sys
from statistics import NormalDist

import pandas as pd

CHI2_MEDIAN_1DF = 0.454936423119572


def chi2_isf_1df(p):
    p = min(max(p, 1e-300), 1.0)
    if p >= 1.0:
        return 0.0
    z = NormalDist().inv_cdf(1.0 - p / 2.0)
    return z * z


def inflation(pvals):
    s = sorted(pvals)
    n = len(s)
    if n == 0:
        return float("nan"), float("nan")
    if n % 2:
        med_stat = chi2_isf_1df(s[n // 2])
        med_p = s[n // 2]
    else:
        med_stat = (chi2_isf_1df(s[n // 2 - 1]) + chi2_isf_1df(s[n // 2])) / 2.0
        med_p = (s[n // 2 - 1] + s[n // 2]) / 2.0
    return med_stat / CHI2_MEDIAN_1DF, med_p


def methods(cons):
    out = []
    for col in cons.columns:
        m = re.fullmatch(r"p_(.+)", col)
        if m and f"q_{m.group(1)}" in cons.columns:
            out.append((m.group(1), col, f"q_{m.group(1)}"))
    for name in ("ebm", "acat"):
        if f"{name}_p" in cons.columns and f"{name}_q" in cons.columns:
            out.append((name, f"{name}_p", f"{name}_q"))
    return out


def element_set(cons, adwgs, min_expected):
    if adwgs is None:
        return None, "tested"
    tab = pd.read_csv(adwgs, sep="\t", usecols=["id", "element_muts_exp"])
    keep = set(tab.loc[pd.to_numeric(tab["element_muts_exp"], errors="coerce") >= min_expected, "id"])
    return keep, f"activedriverwgs_expected_ge_{min_expected:g}"


def qq_points(p, run, method, cohort):
    s = sorted(p)
    n = len(s)
    if n == 0:
        return []
    step = max(1, n // 500)
    rows = []
    for i, v in enumerate(s):
        obs = -math.log10(max(v, 1e-300))
        if obs >= 2 or i % step == 0 or i == n - 1:
            rows.append((cohort, run, method, round(-math.log10((i + 0.5) / n), 6), round(obs, 6)))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--observed-consensus", required=True)
    ap.add_argument("--observed-adwgs")
    ap.add_argument("--null-consensus", nargs="+", required=True, help="<cohort>.null<r>.consensus.tsv files")
    ap.add_argument("--null-adwgs", nargs="*", default=[], help="<cohort>.null<r>.activedriverwgs.tsv files")
    ap.add_argument("--band", default="0.8,1.2")
    ap.add_argument("--fdr", type=float, default=0.1)
    ap.add_argument("--min-expected", type=float, default=1.0)
    ap.add_argument("--calibration", required=True)
    ap.add_argument("--qq", required=True)
    args = ap.parse_args()
    lo, hi = (float(x) for x in args.band.split(","))

    def replicate(path, kind):
        m = re.search(r"\.null(\d+)\." + kind + r"\.tsv$", path)
        if not m:
            sys.exit(f"CALIBRATECALLERS: cannot read the replicate number from {path}")
        return int(m.group(1))

    runs = [("observed", None, args.observed_consensus, args.observed_adwgs)]
    adwgs_by_rep = {replicate(f, "activedriverwgs"): f for f in args.null_adwgs}
    for f in sorted(args.null_consensus, key=lambda x: replicate(x, "consensus")):
        r = replicate(f, "consensus")
        runs.append((f"null{r}", r, f, adwgs_by_rep.get(r)))

    cal, qq = [], []
    for run, rep, cons_path, adwgs in runs:
        cons = pd.read_csv(cons_path, sep="\t")
        keep, set_name = element_set(cons, adwgs, args.min_expected)
        for method, pcol, qcol in methods(cons):
            p = pd.to_numeric(cons[pcol], errors="coerce")
            q = pd.to_numeric(cons[qcol], errors="coerce")
            sel = p.notna() & (cons["element_id"].isin(keep) if keep is not None else True)
            lam, med_p = inflation(p[sel].tolist())
            cal.append({"cohort": args.cohort, "run": run, "replicate": "NA" if rep is None else rep, "method": method,
                        "element_set": set_name, "n_elements": int(sel.sum()), "median_p": med_p, "lambda": lam,
                        "n_q_lt_fdr": int((q < args.fdr).sum()), "fdr": args.fdr,
                        "band_low": lo, "band_high": hi, "status": "NA"})
            qq.extend(qq_points(p[sel].tolist(), run, method, args.cohort))

    cal = pd.DataFrame(cal)
    summary = []
    nulls = cal[cal["run"] != "observed"]
    for method, grp in nulls.groupby("method", sort=False):
        lams = grp["lambda"].dropna()
        if lams.empty:
            status = "not_assessable"
        elif (lams > hi).any():
            status = "excluded_for_inflation"
        elif (lams < lo).any():
            status = "conservative"
        else:
            status = "within_band"
        summary.append({"cohort": args.cohort, "run": "null_summary", "replicate": len(grp), "method": method,
                        "element_set": grp["element_set"].iloc[0], "n_elements": int(grp["n_elements"].median()),
                        "median_p": float(grp["median_p"].median()), "lambda": float(lams.median()) if not lams.empty else float("nan"),
                        "n_q_lt_fdr": int(grp["n_q_lt_fdr"].max()), "fdr": args.fdr, "band_low": lo, "band_high": hi,
                        "status": status})
    out = pd.concat([cal, pd.DataFrame(summary)], ignore_index=True)
    out["lambda_min"] = "NA"
    out["lambda_max"] = "NA"
    for i, row in out[out["run"] == "null_summary"].iterrows():
        lams = nulls.loc[nulls["method"] == row["method"], "lambda"].dropna()
        if not lams.empty:
            out.at[i, "lambda_min"] = f"{lams.min():.6g}"
            out.at[i, "lambda_max"] = f"{lams.max():.6g}"
    for col in ("median_p", "lambda"):
        out[col] = out[col].map(lambda v: "NA" if pd.isna(v) else f"{v:.6g}")
    out.to_csv(args.calibration, sep="\t", index=False)
    pd.DataFrame(qq, columns=["cohort", "run", "method", "expected_mlog10p", "observed_mlog10p"]).to_csv(args.qq, sep="\t", index=False)
    for row in summary:
        print(f"CALIBRATECALLERS {args.cohort} {row['method']}: null lambda {row['lambda']:.4g} ({row['status']}), max q<{args.fdr} {row['n_q_lt_fdr']}", file=sys.stderr)


if __name__ == "__main__":
    main()
