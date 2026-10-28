#!/usr/bin/env python3
"""POSTFILTERFLAGS: PCAWG-style post-filter flags per element (in-house, D33; flag only).

For every element of the consensus table, a value, a threshold and pass | fail | not_applied:
  min_support          n_mutations >= --min-mutations and n_patients >= --min-patients
  mappability          share of the element's mutations outside GIAB low-mappability regions
                       (1 - frac_in_low_mappability) > --min-mappable-fraction
  max_sample_fraction  largest share of the element's mutations from one patient <= --max-sample-fraction
The counts and fractions come from VARIANT_FLAGS (element_flags.tsv); without it, or for an
element without mutations, the flags are not_applied. postfilter_status is fail when any flag
fails, pass when all three pass, otherwise not_applied.
--mode pcawg adds q_postfilter: the p-value of the combination that sets consensus_q (column
combination_method) set to 1 for failing elements, then Benjamini-Hochberg over the elements with
a p-value, as PCAWG did. consensus_q is never changed.
"""
import argparse
import sys

import numpy as np
import pandas as pd


def bh(p):
    p = np.asarray(p, dtype=float)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p, kind="mergesort")
    ranked = p[order] * n / np.arange(1, n + 1)
    q = np.empty(n)
    q[order] = np.minimum.accumulate(ranked[::-1])[::-1]
    return np.minimum(q, 1.0)


def fmt(v):
    return "NA" if pd.isna(v) else f"{v:.6g}"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--consensus", required=True)
    ap.add_argument("--element-flags")
    ap.add_argument("--mode", choices=["flag_only", "pcawg"], default="flag_only")
    ap.add_argument("--min-mutations", type=int, default=3)
    ap.add_argument("--min-patients", type=int, default=3)
    ap.add_argument("--min-mappable-fraction", type=float, default=0.5)
    ap.add_argument("--max-sample-fraction", type=float, default=0.5)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    cons = pd.read_csv(args.consensus, sep="\t", dtype={"element_id": str})
    out = cons[["element_id"]].copy()
    if args.element_flags:
        ef = pd.read_csv(args.element_flags, sep="\t", dtype={"element_id": str}, na_values=["NA"])
        ef = ef[["element_id", "n_mutations", "n_patients", "max_sample_fraction", "frac_in_low_mappability"]]
        out = out.merge(ef, on="element_id", how="left")
    else:
        for col in ("n_mutations", "n_patients", "max_sample_fraction", "frac_in_low_mappability"):
            out[col] = np.nan
    has = out["n_mutations"].fillna(0) > 0
    nm, npat = out["n_mutations"], out["n_patients"]

    out["min_support_threshold"] = f"{args.min_mutations};{args.min_patients}"
    out["min_support_status"] = np.where(nm.isna(), "not_applied",
                                         np.where((nm >= args.min_mutations) & (npat >= args.min_patients), "pass", "fail"))
    out["mappable_fraction"] = 1.0 - out["frac_in_low_mappability"]
    out["mappability_threshold"] = args.min_mappable_fraction
    out["mappability_status"] = np.where(~has | out["mappable_fraction"].isna(), "not_applied",
                                         np.where(out["mappable_fraction"] > args.min_mappable_fraction, "pass", "fail"))
    out["max_sample_fraction_threshold"] = args.max_sample_fraction
    out["max_sample_fraction_status"] = np.where(~has | out["max_sample_fraction"].isna(), "not_applied",
                                                 np.where(out["max_sample_fraction"] <= args.max_sample_fraction, "pass", "fail"))
    st = out[["min_support_status", "mappability_status", "max_sample_fraction_status"]]
    out["postfilter_status"] = np.where((st == "fail").any(axis=1), "fail",
                                        np.where((st == "pass").all(axis=1), "pass", "not_applied"))

    cols = ["element_id", "n_mutations", "n_patients", "min_support_threshold", "min_support_status",
            "mappable_fraction", "mappability_threshold", "mappability_status",
            "max_sample_fraction", "max_sample_fraction_threshold", "max_sample_fraction_status", "postfilter_status"]
    if args.mode == "pcawg":
        methods = cons["combination_method"].dropna().unique()
        if len(methods) != 1 or methods[0] not in ("ebm", "acat"):
            sys.exit(f"POSTFILTERFLAGS: consensus table has combination_method {list(methods)}; expected one of ebm, acat")
        p = pd.to_numeric(cons[f"{methods[0]}_p"], errors="coerce").to_numpy()
        p = np.where(out["postfilter_status"].to_numpy() == "fail", np.where(np.isnan(p), np.nan, 1.0), p)
        q = np.full(p.size, np.nan)
        ok = ~np.isnan(p)
        q[ok] = bh(p[ok])
        out["q_postfilter"] = q
        cols.append("q_postfilter")
    for col in ("n_mutations", "n_patients"):
        out[col] = out[col].map(lambda v: "NA" if pd.isna(v) else str(int(v)))
    for col in ("mappable_fraction", "max_sample_fraction", "q_postfilter"):
        if col in out:
            out[col] = out[col].map(fmt)
    out[cols].to_csv(args.output, sep="\t", index=False)
    counts = out["postfilter_status"].value_counts().to_dict()
    print(f"POSTFILTERFLAGS mode={args.mode} " + " ".join(f"{k}={v}" for k, v in sorted(counts.items())), file=sys.stderr)


if __name__ == "__main__":
    main()
