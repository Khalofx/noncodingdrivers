#!/usr/bin/env python3
"""WRITEDRIVEREVIDENCE: long-form evidence table, one row per cohort x element x evidence item (D33).

Columns: cohort, element_id, evidence_item, value, unit, status, threshold, source, source_version,
native_or_lifted (validated by assets/schema_driver_evidence.json). Sources, units and versions come
from assets/evidence_sources.tsv; "pipeline" versions are replaced by --pipeline-version.
native_or_lifted is the cohort's coordinate provenance (--coordinates) for every item computed from
its mutations, and native for items computed from the GRCh38 resources alone (callable_fraction).

Elements written (--scope): nominal = any caller or combination p < 0.01 or consensus_q < 0.25
(the entry rule of the lowest downstream tier); all = every element of the consensus table.
Statuses: significant | not_significant (consensus_q and caller q against --fdr), observed,
not_tested (NA p), flagged | not_flagged (mask fractions above --mask-fraction, germline share
above it), pass | fail | not_applied (post-filters), and the null-calibration band status. No
evidence item is a tier: tiers are assigned downstream (reports repository).
"""
import argparse
import sys

import numpy as np
import pandas as pd

CALLER_ITEMS = ("p_activedriverwgs", "q_activedriverwgs", "p_oncodrivefml", "q_oncodrivefml")
CONS_ITEMS = ("consensus_q", "consensus_rank", "ebm_p", "ebm_q", "acat_p", "acat_q", "n_tested_callers", "vote_count")
FLAG_ITEMS = ("n_mutations", "n_patients", "max_sample_fraction", "frac_in_blacklist_v2", "frac_in_segdup",
              "frac_in_low_mappability", "mean_umap_unique_fraction", "frac_germline_pop_common", "callable_fraction")
MASK_ITEMS = ("frac_in_blacklist_v2", "frac_in_segdup", "frac_in_low_mappability", "frac_germline_pop_common")


def fmt(v):
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return "NA"
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    try:
        f = float(v)
    except (TypeError, ValueError):
        return str(v)
    if np.isnan(f):
        return "NA"
    return str(int(f)) if f.is_integer() and abs(f) < 1e15 else f"{f:.6g}"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--consensus", required=True)
    ap.add_argument("--element-flags")
    ap.add_argument("--postfilter")
    ap.add_argument("--calibration")
    ap.add_argument("--sources", required=True)
    ap.add_argument("--coordinates", choices=["native", "lifted", "mixed"], required=True)
    ap.add_argument("--pipeline-version", required=True)
    ap.add_argument("--scope", choices=["nominal", "all"], default="nominal")
    ap.add_argument("--fdr", type=float, default=0.1)
    ap.add_argument("--mask-fraction", type=float, default=0.5)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    src = pd.read_csv(args.sources, sep="\t", dtype=str).set_index("evidence_item")
    src["source_version"] = src["source_version"].replace("pipeline", f"noncodingdrivers {args.pipeline_version}")
    cons = pd.read_csv(args.consensus, sep="\t", dtype={"element_id": str})
    pcols = [c for c in cons.columns if c in CALLER_ITEMS[0::2] or c in ("ebm_p", "acat_p")]
    minp = cons[pcols].apply(pd.to_numeric, errors="coerce").min(axis=1)
    cq = pd.to_numeric(cons["consensus_q"], errors="coerce")
    keep = (minp < 0.01) | (cq < 0.25) if args.scope == "nominal" else pd.Series(True, index=cons.index)
    cons = cons[keep]
    ids = cons["element_id"]
    ef = pd.read_csv(args.element_flags, sep="\t", dtype={"element_id": str}).set_index("element_id") if args.element_flags else None
    pf = pd.read_csv(args.postfilter, sep="\t", dtype=str).set_index("element_id") if args.postfilter else None

    rows = []

    def add(eid, item, value, status, threshold="NA"):
        s = src.loc[item]
        coords = "native" if s["coordinates"] == "native" else args.coordinates
        rows.append((args.cohort, eid, item, fmt(value), s["unit"], status, threshold, s["source"], s["source_version"], coords))

    for _, r in cons.iterrows():
        eid = r["element_id"]
        for item in CONS_ITEMS + CALLER_ITEMS:
            if item not in cons.columns:
                continue
            v = pd.to_numeric(r[item], errors="coerce")
            if item in ("consensus_q",) or item.startswith("q_") or item.endswith("_q"):
                status = "not_tested" if pd.isna(v) else ("significant" if v < args.fdr else "not_significant")
                add(eid, item, v, status, fmt(args.fdr))
            elif item.startswith("p_") or item.endswith("_p"):
                add(eid, item, v, "not_tested" if pd.isna(v) else "observed")
            else:
                add(eid, item, v, "not_tested" if pd.isna(v) else "observed")
        if ef is not None:
            fr = ef.loc[eid] if eid in ef.index else None
            for item in FLAG_ITEMS:
                v = pd.to_numeric(fr[item], errors="coerce") if fr is not None else np.nan
                if pd.isna(v):
                    add(eid, item, v, "not_applied")
                elif item in MASK_ITEMS:
                    add(eid, item, v, "flagged" if v > args.mask_fraction else "not_flagged", fmt(args.mask_fraction))
                else:
                    add(eid, item, v, "observed")
        if pf is not None and eid in pf.index:
            p = pf.loc[eid]
            ms = min(pd.to_numeric(p["n_mutations"], errors="coerce"), pd.to_numeric(p["n_patients"], errors="coerce"))
            add(eid, "postfilter_min_support", ms, p["min_support_status"], p["min_support_threshold"])
            add(eid, "postfilter_mappability", p["mappable_fraction"], p["mappability_status"], fmt(p["mappability_threshold"]))
            add(eid, "postfilter_max_sample_fraction", p["max_sample_fraction"], p["max_sample_fraction_status"], fmt(p["max_sample_fraction_threshold"]))
            add(eid, "postfilter_status", "NA", p["postfilter_status"])
            if "q_postfilter" in pf.columns:
                v = pd.to_numeric(p["q_postfilter"], errors="coerce")
                add(eid, "q_postfilter", v, "not_tested" if pd.isna(v) else ("significant" if v < args.fdr else "not_significant"), fmt(args.fdr))
    if args.calibration:
        cal = pd.read_csv(args.calibration, sep="\t", dtype=str)
        for _, r in cal[cal["run"] == "null_summary"].iterrows():
            s = src.loc["null_lambda"]
            rows.append((args.cohort, ".", f"null_lambda_{r['method']}", fmt(r["lambda"]), s["unit"], r["status"],
                         f"{r['band_low']},{r['band_high']}", s["source"], s["source_version"], args.coordinates))

    cols = ["cohort", "element_id", "evidence_item", "value", "unit", "status", "threshold", "source", "source_version", "native_or_lifted"]
    pd.DataFrame(rows, columns=cols).to_csv(args.output, sep="\t", index=False)
    print(f"WRITEDRIVEREVIDENCE {args.cohort}: {len(ids)} elements ({args.scope}), {len(rows)} rows", file=sys.stderr)


if __name__ == "__main__":
    main()
