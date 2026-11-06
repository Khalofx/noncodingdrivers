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

OncodriveCLUSTL (--oncodriveclustl, its native table): oncodriveclustl_n_clusters (CLUSTERS) and
oncodriveclustl_clustered_mutations (CLUSTERED_MUT) per element; not_tested for an element absent from
the table (below the mutation cutoff).

Callers: every p_<caller>/q_<caller> pair of the consensus table is written, combined callers and
callers reported alongside alike (known callers first, in CALLER_ORDER). --consensus-callers names
the combined callers (the pipeline's --consensus_callers); when at least one caller of the table is
not among them, one cohort-level row per caller (element_id '.', caller_role_<caller>) records its
role: value 1 status combined, or value 0 status reported_alongside. When every caller is combined
(the default) no role rows are written, so the default table is unchanged.

Annotation and context layers (S7A.33): for every written element with is_consensus True (the
elements the layers annotate), a fixed set of summary items per layer whose table is given; a
layer that is off gives no table and no row, so the default table is unchanged. The items are
element-level summaries of the layer's own table (no per-variant or per-link rows); evidence never
filters or re-ranks. A significant element missing from a given layer table gets the layer's items
with value NA, status not_applied.
  --avi-summary      avi_max_phred (score), avi_n_scored, avi_n_unscored (candidate variants)
  --opencravat       opencravat_n_<annotator> for ccre_screen, repeat, clinvar: candidate variants
                     with at least one non-NA column of that annotator
  --cosmic-elements  cosmic_n_matched_mutations, cosmic_n_samples, cosmic_n_sample_name_overlaps
                     (counts only; no COSMIC identifier or record is read or written)
  --cancer-genes     cancer_gene_nearest_in_list, cancer_gene_nearest_pc_in_list (1/0, in_list or
                     not_in_list), cancer_gene_re2g_targets_in_list (count; NA without ENCODE-rE2G)
  --re2g-summary     re2g_n_distal_targets (distinct distal target genes over the requested tissues),
                     re2g_top_distal_target (gene label; status = nearest_gene_status over the tissues)
  --ccre-summary     ccre_v4_dominant_class (class label), ccre_v4_covered_fraction
Null false-positive rate (--null-fpr, CALIBRATECALLERS <cohort>.null_fpr.tsv): one cohort-level row
per method of its null_summary (null_fpr_<method>, OncodriveCLUSTL): the share of tested null
elements at p <= 0.05 over the replicates.
"""
import argparse
import sys

import numpy as np
import pandas as pd

CALLER_ORDER = ("activedriverwgs", "oncodrivefml", "oncodriveclustl")
# Per-element counts from the native OncodriveCLUSTL table (--oncodriveclustl): item -> column
CLUSTL_ITEMS = (("oncodriveclustl_n_clusters", "CLUSTERS"), ("oncodriveclustl_clustered_mutations", "CLUSTERED_MUT"))
CONS_ITEMS = ("consensus_q", "consensus_rank", "ebm_p", "ebm_q", "acat_p", "acat_q", "n_tested_callers", "vote_count")
FLAG_ITEMS = ("n_mutations", "n_patients", "max_sample_fraction", "frac_in_blacklist_v2", "frac_in_segdup",
              "frac_in_low_mappability", "mean_umap_unique_fraction", "frac_germline_pop_common", "callable_fraction")
MASK_ITEMS = ("frac_in_blacklist_v2", "frac_in_segdup", "frac_in_low_mappability", "frac_germline_pop_common")
OPENCRAVAT_ANNOTATORS = ("ccre_screen", "repeat", "clinvar")
RE2G_LINKED = ("top_is_nearest", "nearest_in_targets")


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


class _Fill(dict):
    """source_version placeholders ({release}, {list}, {md5}, {tissues}); NA when the table had none."""

    def __missing__(self, key):
        return "NA"


def read_table(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def missing(v):
    return v is None or v in ("", "NA", "nan")


def layer_rows(args, sig_ids):
    """Per significant element: [(item, value, status, threshold)] of every layer whose table is given,
    and the source_version placeholders of each layer."""
    out = {eid: [] for eid in sig_ids}
    fills = {}
    if args.avi_summary:
        t = read_table(args.avi_summary).set_index("element_id")
        for eid in sig_ids:
            if eid not in t.index:
                out[eid] += [(i, "NA", "not_applied", "NA") for i in ("avi_max_phred", "avi_n_scored", "avi_n_unscored")]
                continue
            r = t.loc[eid]
            n_cand, n_scored = int(r["n_candidate_variants"]), int(r["n_scored"])
            mp = r["max_avi_phred"]
            out[eid] += [("avi_max_phred", "NA" if missing(mp) else float(mp), "not_applied" if missing(mp) else "observed", "NA"),
                         ("avi_n_scored", n_scored, "observed", "NA"),
                         ("avi_n_unscored", n_cand - n_scored, "observed", "NA")]
    if args.opencravat:
        t = read_table(args.opencravat)
        hits = {eid: None for eid in sig_ids}
        cols = {a: [c for c in t.columns if c.startswith(f"{a}__")] for a in OPENCRAVAT_ANNOTATORS}
        for _, r in t.iterrows():
            for eid in r["element_ids"].split(","):
                if eid not in hits:
                    continue
                if hits[eid] is None:
                    hits[eid] = {a: 0 for a in OPENCRAVAT_ANNOTATORS}
                for a in OPENCRAVAT_ANNOTATORS:
                    if any(not missing(r[c]) for c in cols[a]):
                        hits[eid][a] += 1
        for eid in sig_ids:
            for a in OPENCRAVAT_ANNOTATORS:
                if hits[eid] is None or not cols[a]:
                    out[eid].append((f"opencravat_n_{a}", "NA", "not_applied", "NA"))
                else:
                    out[eid].append((f"opencravat_n_{a}", hits[eid][a], "observed", "NA"))
    if args.cosmic_elements:
        t = read_table(args.cosmic_elements).set_index("element_id")
        releases = sorted(set(t["cosmic_release"])) if len(t) else []
        fills["release"] = ",".join(releases) if releases else "NA"
        items = (("cosmic_n_matched_mutations", "n_cosmic_matched_mutations"), ("cosmic_n_samples", "n_cosmic_samples"),
                 ("cosmic_n_sample_name_overlaps", "n_exact_sample_name_overlaps"))
        for eid in sig_ids:
            for item, col in items:
                if eid in t.index and not missing(t.loc[eid, col]):
                    out[eid].append((item, int(t.loc[eid, col]), "observed", "NA"))
                else:
                    out[eid].append((item, "NA", "not_applied", "NA"))
    if args.cancer_genes:
        t = read_table(args.cancer_genes).set_index("element_id")
        if len(t):
            fills["list"] = t["list_name"].iloc[0]
            fills["md5"] = t["list_md5"].iloc[0]
        for eid in sig_ids:
            r = t.loc[eid] if eid in t.index else None
            for item, col in (("cancer_gene_nearest_in_list", "nearest_in_list"), ("cancer_gene_nearest_pc_in_list", "nearest_pc_in_list")):
                v = r[col] if r is not None else "NA"
                if v == "yes":
                    out[eid].append((item, 1, "in_list", "NA"))
                elif v == "no":
                    out[eid].append((item, 0, "not_in_list", "NA"))
                else:
                    out[eid].append((item, "NA", "not_applied", "NA"))
            v = r["re2g_targets_in_list"] if r is not None else "NA"
            if r is None or v == "NA":
                out[eid].append(("cancer_gene_re2g_targets_in_list", "NA", "not_applied", "NA"))
            else:
                out[eid].append(("cancer_gene_re2g_targets_in_list", len([g for g in v.split(",") if g]), "observed", "NA"))
    if args.re2g_summary:
        t = read_table(args.re2g_summary)
        tissues = list(dict.fromkeys(t["tissue"]))
        fills["tissues"] = ",".join(tissues) if tissues else "NA"
        for eid in sig_ids:
            rows = t[t["element_id"] == eid]
            if rows.empty:
                out[eid] += [("re2g_n_distal_targets", "NA", "not_applied", "NA"), ("re2g_top_distal_target", "NA", "not_applied", "NA")]
                continue
            genes, best = set(), None
            for _, r in rows.iterrows():
                genes |= {x.rsplit(":", 1)[0] for x in r["distal_targets_by_max_score"].split(",") if x}
                if not missing(r["top_distal_target"]):
                    cand = (-float(r["top_distal_score"]), r["top_distal_target"], r["nearest_gene_status"])
                    best = cand if best is None or cand < best else best
            statuses = set(rows["nearest_gene_status"])
            if "no_nearest_gene" in statuses:
                status = "no_nearest_gene"
            elif not genes:
                status = "no_link"
            elif best[2] == "top_is_nearest":
                status = "top_is_nearest"
            elif statuses & set(RE2G_LINKED):
                status = "nearest_in_targets"
            else:
                status = "nearest_not_in_targets"
            out[eid] += [("re2g_n_distal_targets", len(genes), "observed", "NA"),
                         ("re2g_top_distal_target", best[1] if best else "NA", status, "NA")]
    if args.ccre_summary:
        t = read_table(args.ccre_summary).set_index("element_id")
        for eid in sig_ids:
            if eid not in t.index:
                out[eid] += [("ccre_v4_dominant_class", "NA", "not_applied", "NA"), ("ccre_v4_covered_fraction", "NA", "not_applied", "NA")]
                continue
            r = t.loc[eid]
            out[eid] += [("ccre_v4_dominant_class", r["dominant_class"], "observed", "NA"),
                         ("ccre_v4_covered_fraction", float(r["covered_frac"]), "observed", "NA")]
    return out, fills


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
    ap.add_argument("--consensus-callers", help="comma-separated combined callers; default: every caller of the table")
    ap.add_argument("--oncodriveclustl", help="native OncodriveCLUSTL table: cluster counts per element")
    ap.add_argument("--avi-summary", help="EXPLAIN <cohort>.avi.element_summary.tsv (--run_alphagenome_avi)")
    ap.add_argument("--opencravat", help="EXPLAIN <cohort>.opencravat.tsv (--run_opencravat)")
    ap.add_argument("--cosmic-elements", help="CONTEXT_LAYERS <cohort>.cosmic.element_summary.tsv (--run_cosmic_context; counts only)")
    ap.add_argument("--cancer-genes", help="CONTEXT_LAYERS <cohort>.cancer_gene_context.tsv (--cancer_gene_list)")
    ap.add_argument("--re2g-summary", help="CONTEXT_LAYERS <cohort>.enhancer_gene_summary.tsv (--run_enhancer_gene_links)")
    ap.add_argument("--ccre-summary", help="CONTEXT_LAYERS <cohort>.ccre_v4_summary.tsv (--run_ccre_overlay)")
    ap.add_argument("--null-fpr", help="CALIBRATECALLERS <cohort>.null_fpr.tsv (null false-positive rate of OncodriveCLUSTL)")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    src = pd.read_csv(args.sources, sep="\t", dtype=str).set_index("evidence_item")
    src["source_version"] = src["source_version"].replace("pipeline", f"noncodingdrivers {args.pipeline_version}")
    cons = pd.read_csv(args.consensus, sep="\t", dtype={"element_id": str})
    callers = [c[2:] for c in cons.columns if c.startswith("p_")]
    callers = [c for c in CALLER_ORDER if c in callers] + sorted(c for c in callers if c not in CALLER_ORDER)
    caller_items = tuple(item for c in callers for item in (f"p_{c}", f"q_{c}"))
    pcols = [c for c in cons.columns if c in caller_items[0::2] or c in ("ebm_p", "acat_p")]
    minp = cons[pcols].apply(pd.to_numeric, errors="coerce").min(axis=1)
    cq = pd.to_numeric(cons["consensus_q"], errors="coerce")
    keep = (minp < 0.01) | (cq < 0.25) if args.scope == "nominal" else pd.Series(True, index=cons.index)
    cons = cons[keep]
    ids = cons["element_id"]
    ef = pd.read_csv(args.element_flags, sep="\t", dtype={"element_id": str}).set_index("element_id") if args.element_flags else None
    pf = pd.read_csv(args.postfilter, sep="\t", dtype=str).set_index("element_id") if args.postfilter else None
    cl = (pd.read_csv(args.oncodriveclustl, sep="\t", dtype=str, keep_default_na=False).set_index("ENSID")
          if args.oncodriveclustl else None)

    if "is_consensus" in cons.columns:
        sig = cons["is_consensus"].astype(str) == "True"
    else:
        sig = pd.to_numeric(cons["consensus_q"], errors="coerce") < args.fdr
    layers, fills = layer_rows(args, list(cons.loc[sig, "element_id"]))

    rows = []

    def add(eid, item, value, status, threshold="NA"):
        s = src.loc[item]
        coords = "native" if s["coordinates"] == "native" else args.coordinates
        rows.append((args.cohort, eid, item, fmt(value), s["unit"], status, threshold, s["source"], s["source_version"], coords))

    for _, r in cons.iterrows():
        eid = r["element_id"]
        for item in CONS_ITEMS + caller_items:
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
        if cl is not None:
            # elements below OncodriveCLUSTL's mutation cutoff are not in its table: not tested
            cr = cl.loc[eid] if eid in cl.index else None
            for item, col in CLUSTL_ITEMS:
                v = pd.to_numeric(cr[col], errors="coerce") if cr is not None and cr[col] != "" else np.nan
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
        for item, value, status, threshold in layers.get(eid, []):
            s = src.loc[item]
            coords = "native" if s["coordinates"] == "native" else args.coordinates
            rows.append((args.cohort, eid, item, fmt(value), s["unit"], status, threshold, s["source"],
                         s["source_version"].format_map(_Fill(fills)), coords))
    named = callers if args.consensus_callers is None else [c for c in args.consensus_callers.split(",") if c]
    if any(c not in named for c in callers):
        s = src.loc["caller_role"]
        for c in callers:
            role = "combined" if c in named else "reported_alongside"
            rows.append((args.cohort, ".", f"caller_role_{c}", "1" if role == "combined" else "0", s["unit"], role,
                         "NA", s["source"], s["source_version"], args.coordinates))
    if args.calibration:
        cal = pd.read_csv(args.calibration, sep="\t", dtype=str)
        for _, r in cal[cal["run"] == "null_summary"].iterrows():
            s = src.loc["null_lambda"]
            rows.append((args.cohort, ".", f"null_lambda_{r['method']}", fmt(r["lambda"]), s["unit"], r["status"],
                         f"{r['band_low']},{r['band_high']}", s["source"], s["source_version"], args.coordinates))
    if args.null_fpr:
        fpr = pd.read_csv(args.null_fpr, sep="\t", dtype=str, keep_default_na=False)
        s = src.loc["null_fpr"]
        for _, r in fpr[fpr["run"] == "null_summary"].iterrows():
            v = pd.to_numeric(r["fpr_p_le_0.05"], errors="coerce")
            rows.append((args.cohort, ".", f"null_fpr_{r['method']}", fmt(v), s["unit"], "not_tested" if pd.isna(v) else "observed",
                         "0.05", s["source"], s["source_version"], args.coordinates))

    cols = ["cohort", "element_id", "evidence_item", "value", "unit", "status", "threshold", "source", "source_version", "native_or_lifted"]
    pd.DataFrame(rows, columns=cols).to_csv(args.output, sep="\t", index=False)
    print(f"WRITEDRIVEREVIDENCE {args.cohort}: {len(ids)} elements ({args.scope}), {len(rows)} rows", file=sys.stderr)


if __name__ == "__main__":
    main()
