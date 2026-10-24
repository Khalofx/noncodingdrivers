#!/usr/bin/env python3
"""Compare a v1 consensus table with the regression baseline (tests/regression/).

Gate tooling only: Nextflow never runs this script.

Checks, one row each in the report:
  * number of significant elements (method q < --fdr) against baseline.json;
  * the named sentinel elements of baseline.json (rank, q, per-caller p) for the cohort;
  * optionally, with --frozen, the significant set and per-caller p of every shared element
    against a frozen consensus table (column names mapped with --map v1col=frozencol);
  * optionally, with --comparison, the number of tested elements and the EBM/ACAT overlap counts
    of combination_method_comparison.csv;
  * optionally, with --discordant, the ACAT-only and EBM-only element sets of
    combination_discordant_elements.csv, and their per-caller pattern (every ACAT-only element has
    OncodriveFML p > 0.05; every EBM-only element has both callers at p < 0.05).

Significance is q < --fdr (default 0.1: the release threshold, and the pipeline's --fdr_threshold).
The release's ACAT values came from an in-house ACAT; the package ACAT matches them to --acat-q-rtol.

Exit status 0 if every check passes, 1 otherwise. The report is tab-separated.
"""

import argparse
import json
import math
import sys

import pandas as pd

COHORT_KEYS = {
    "liver_hcc_public": "significant_annotations_acat",
    "liri_jp": "candidates_acat",
    "paca_ca": "candidates_acat",
}
# Cohort names as written in the combination-check CSVs
COHORT_NAMES = {"liver_hcc_public": "Liver-HCC-public", "liri_jp": "LIRI-JP", "paca_ca": "PACA-CA"}


def close(observed, expected, rtol):
    if expected is None or observed is None:
        return False
    if isinstance(observed, float) and math.isnan(observed):
        return False
    return math.isclose(float(observed), float(expected), rel_tol=rtol, abs_tol=0.0)


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--baseline", required=True, help="tests/regression/baseline.json")
    ap.add_argument("--cohort", required=True, choices=sorted(COHORT_KEYS))
    ap.add_argument("--v1", required=True, help="v1 consensus table (COMBINE_PVALUES output, TSV)")
    ap.add_argument("--method", default="ebm", choices=["ebm", "acat"])
    ap.add_argument("--fdr", type=float, default=0.1, help="q threshold of the significant sets (release: 0.1)")
    ap.add_argument("--p-rtol", type=float, default=1e-9, help="relative tolerance for deterministic p-values")
    ap.add_argument("--ofml-rtol", type=float, default=0.0,
                    help="relative tolerance for OncodriveFML p (0 = exact; raise only if the seed is not fixed)")
    ap.add_argument("--acat-q-rtol", type=float, default=1e-5,
                    help="relative tolerance for the TERT ACAT q of baseline.json, which the release computed with its "
                         "in-house ACAT; the package ACAT (p = 1 inputs removed, not clamped) differs from it by at most "
                         "3.81e-06 in log10 on Liver-HCC-public (combination_method_comparison.csv)")
    ap.add_argument("--q-rtol", type=float, default=1e-2,
                    help="relative tolerance for q values recorded rounded in baseline.json")
    ap.add_argument("--frozen", help="frozen consensus table (TSV) for set and per-element comparison")
    ap.add_argument("--map", action="append", default=[],
                    help="v1col=frozencol column mapping for --frozen (repeatable); identity if omitted")
    ap.add_argument("--frozen-q", default="consensus_q", help="q column of the frozen table for its significant set")
    ap.add_argument("--comparison", help="combination_method_comparison.csv (tested elements, EBM/ACAT overlap counts)")
    ap.add_argument("--discordant", help="combination_discordant_elements.csv (ACAT-only and EBM-only elements)")
    ap.add_argument("--out", default="-", help="report path (default stdout)")
    return ap.parse_args(argv)


def main(argv=None):
    a = parse_args(argv)
    base = json.load(open(a.baseline))
    v1 = pd.read_csv(a.v1, sep="\t", dtype={"element_id": str})
    qcol, pcol = f"{a.method}_q", f"{a.method}_p"
    for col in ("element_id", qcol, pcol):
        if col not in v1.columns:
            sys.exit(f"column {col} missing from {a.v1}")
    v1 = v1.assign(_rank=v1[pcol].rank(method="min"))
    rows = []

    def check(name, expected, observed, ok):
        rows.append({"check": name, "expected": expected, "observed": observed, "status": "PASS" if ok else "FAIL"})

    sig = set(v1.loc[v1[qcol] < a.fdr, "element_id"])
    cc = base.get("combination_check_2026_10_06", {})
    expected_n = (cc.get("ebm_significant", {}) if a.method == "ebm" else cc.get("acat_package_significant", {})).get(a.cohort)
    if expected_n is None and a.method == "acat":
        expected_n = base[a.cohort].get(COHORT_KEYS[a.cohort])
    check(f"n_significant_{a.method}", expected_n, len(sig), expected_n == len(sig))

    tert = base[a.cohort].get("tert")
    if tert:
        hit = v1[v1["element_id"] == tert["element"]]
        if hit.empty:
            check("tert_present", tert["element"], "absent", False)
        else:
            r = hit.iloc[0]
            exp_rank = tert["rank"] if a.method == "acat" else cc.get("tert_liver_hcc", {}).get("ebm_rank", tert["rank"])
            check(f"tert_rank_{a.method}", exp_rank, int(r["_rank"]), int(r["_rank"]) == exp_rank)
            if a.method == "acat":
                check("tert_acat_q", tert["acat_q"], r[qcol], close(r[qcol], tert["acat_q"], a.acat_q_rtol))
            else:
                exp_q = cc.get("tert_liver_hcc", {}).get("ebm_q")
                check("tert_ebm_q", exp_q, r[qcol], close(r[qcol], exp_q, a.q_rtol))
            if "p_activedriverwgs" in v1.columns:
                check("tert_p_activedriverwgs", tert["adwgs_p"], r["p_activedriverwgs"],
                      close(r["p_activedriverwgs"], tert["adwgs_p"], a.p_rtol))
            if "p_oncodrivefml" in v1.columns:
                ok = (float(r["p_oncodrivefml"]) == float(tert["ofml_p"])) if a.ofml_rtol == 0 else \
                    close(r["p_oncodrivefml"], tert["ofml_p"], a.ofml_rtol)
                check("tert_p_oncodrivefml", tert["ofml_p"], r["p_oncodrivefml"], ok)
    if "min_acat_q" in base[a.cohort] and a.method == "acat":
        exp = base[a.cohort]["min_acat_q"]
        obs = float(v1[qcol].min())
        check("min_acat_q", exp, obs, close(obs, exp, a.q_rtol))

    if a.frozen:
        fr = pd.read_csv(a.frozen, sep="\t", dtype=str)
        mapping = dict(m.split("=", 1) for m in a.map)
        fid = mapping.get("element_id", "element_id")
        fsig = set(fr.loc[pd.to_numeric(fr[a.frozen_q], errors="coerce") < a.fdr, fid])
        check("frozen_significant_set_equal", len(fsig), len(sig), fsig == sig)
        check("frozen_significant_only_in_v1", 0, len(sig - fsig), not (sig - fsig))
        check("frozen_significant_only_in_frozen", 0, len(fsig - sig), not (fsig - sig))
        joined = v1.merge(fr, left_on="element_id", right_on=fid, how="inner", suffixes=("", "_frozen"))
        check("frozen_shared_elements", len(fr), len(joined), len(joined) == len(fr))
        for v1col in [c for c in v1.columns if c.startswith("p_")]:
            fcol = mapping.get(v1col, v1col)
            if fcol not in fr.columns:
                continue
            fcol_j = fcol if fcol in joined.columns and fcol != v1col else f"{fcol}_frozen"
            fv = pd.to_numeric(joined[fcol_j], errors="coerce")
            rtol = a.ofml_rtol if v1col == "p_oncodrivefml" else a.p_rtol
            bad = sum(not (pd.isna(x) and pd.isna(y)) and not close(x, y, rtol) and x != y
                      for x, y in zip(joined[v1col], fv))
            check(f"frozen_{v1col}_mismatches", 0, int(bad), bad == 0)

    if a.comparison or a.discordant:
        for col in ("ebm_q", "acat_q"):
            if col not in v1.columns:
                sys.exit(f"column {col} missing from {a.v1} (needed for --comparison/--discordant)")
        ebm_sig = set(v1.loc[v1["ebm_q"] < a.fdr, "element_id"])
        acat_sig = set(v1.loc[v1["acat_q"] < a.fdr, "element_id"])
        name = COHORT_NAMES[a.cohort]
    if a.comparison:
        cmp_rows = pd.read_csv(a.comparison)
        cmp_rows = cmp_rows[cmp_rows["cohort"] == name]
        if cmp_rows.empty:
            check("comparison_cohort_present", name, "absent", False)
        else:
            c = cmp_rows.iloc[0]
            check("n_elements_tested", int(c["elements"]), len(v1), int(c["elements"]) == len(v1))
            check("n_shared_ebm_acat", int(c["shared_ebm"]), len(ebm_sig & acat_sig), int(c["shared_ebm"]) == len(ebm_sig & acat_sig))
            check("n_acat_only", int(c["acat_only"]), len(acat_sig - ebm_sig), int(c["acat_only"]) == len(acat_sig - ebm_sig))
            check("n_ebm_only", int(c["ebm_only"]), len(ebm_sig - acat_sig), int(c["ebm_only"]) == len(ebm_sig - acat_sig))
    if a.discordant:
        dis = pd.read_csv(a.discordant)
        dis = dis[dis["cohort"] == name]
        for status, observed in (("ACAT_only", acat_sig - ebm_sig), ("EBM_only", ebm_sig - acat_sig)):
            expected = set(dis.loc[dis["status"] == status, "element_id"])
            check(f"discordant_{status.lower()}_set_equal", len(expected), len(observed), expected == observed)
        by_id = v1.set_index("element_id")
        acat_only = sorted(acat_sig - ebm_sig)
        ebm_only = sorted(ebm_sig - acat_sig)
        if "p_oncodrivefml" in v1.columns:
            n_bad = sum(not (pd.to_numeric(by_id.at[e, "p_oncodrivefml"], errors="coerce") > 0.05) for e in acat_only)
            check("acat_only_all_oncodrivefml_p_gt_0.05", 0, n_bad, n_bad == 0)
        if {"p_oncodrivefml", "p_activedriverwgs"} <= set(v1.columns):
            n_bad = sum(not all(pd.to_numeric(by_id.at[e, col], errors="coerce") < 0.05
                                for col in ("p_activedriverwgs", "p_oncodrivefml")) for e in ebm_only)
            check("ebm_only_all_both_callers_p_lt_0.05", 0, n_bad, n_bad == 0)

    report = pd.DataFrame(rows, columns=["check", "expected", "observed", "status"])
    report.to_csv(sys.stdout if a.out == "-" else a.out, sep="\t", index=False)
    return 0 if (report["status"] == "PASS").all() else 1


if __name__ == "__main__":
    sys.exit(main())
