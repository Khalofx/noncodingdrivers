#!/usr/bin/env python3
"""Merge the OncodriveCLUSTL tables of the shards of one cohort (ONCODRIVECLUSTL, >= 2 shards).

Every shard ran OncodriveCLUSTL 1.2.0 on a disjoint set of element ids with the same cohort
signature. This script only re-does, over the union of the shards, the two steps OncodriveCLUSTL
itself applies to the elements of one run (oncodriveclustl/utils/postprocessing.py):

  elements  Q_ANALYTICAL, Q_EMPIRICAL, Q_TOPCLUSTER = Benjamini-Hochberg of the matching P column
            (statsmodels multipletests, method fdr_bh, the call OncodriveCLUSTL makes) over the
            elements with a finite P_ANALYTICAL; empty elsewhere. Rows sorted by Q_ANALYTICAL,
            P_ANALYTICAL (ascending), SCORE, CGC (descending), then ENSID for a fixed order.
  clusters  RANK = rank of the element in that order; rows sorted by RANK, P (ascending), SCORE
            (descending).

Every other value is copied as the shard wrote it (read as text), so a value keeps its formatting.
Runs in the OncodriveCLUSTL image (pandas, numpy, statsmodels).

Usage: merge_oncodriveclustl_shards.py --elements a/elements_results.txt b/elements_results.txt
           --clusters a/clusters_results.tsv b/clusters_results.tsv
           --out-elements X.oncodriveclustl.tsv --out-clusters X.oncodriveclustl.clusters.tsv
"""
import argparse
import sys

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

ELEMENT_COLUMNS = ["SYMBOL", "ENSID", "CGC", "CHROMOSOME", "STRAND", "LENGTH", "TOTAL_MUT", "CLUSTERED_MUT",
                   "CLUSTERS", "SIM_CLUSTERS", "SCORE", "P_EMPIRICAL", "Q_EMPIRICAL", "P_ANALYTICAL",
                   "Q_ANALYTICAL", "P_TOPCLUSTER", "Q_TOPCLUSTER"]
CLUSTER_COLUMNS = ["RANK", "SYMBOL", "ENSID", "CGC", "CHROMOSOME", "STRAND", "COORDINATES", "MAX_COORD", "WIDTH",
                   "N_MUT", "N_SAMPLES", "FRA_UNIQ_SAMPLES", "SCORE", "P"]


def read_text(path, columns):
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if list(df.columns) != columns:
        sys.exit(f"ERROR: {path}: unexpected columns {list(df.columns)}")
    return df


def num(series):
    # float() per value: exact round trip of the shard's text (pandas' fast parsers are not)
    return series.map(lambda v: float(v) if v != "" else np.nan).astype(float)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--elements", nargs="+", required=True)
    ap.add_argument("--clusters", nargs="+", required=True)
    ap.add_argument("--out-elements", required=True)
    ap.add_argument("--out-clusters", required=True)
    a = ap.parse_args()

    el = pd.concat([read_text(p, ELEMENT_COLUMNS) for p in a.elements], ignore_index=True)
    if el["ENSID"].duplicated().any():
        sys.exit("ERROR: an element id appears in more than one shard: " + ", ".join(el.loc[el["ENSID"].duplicated(), "ENSID"].head(5)))

    tested = np.isfinite(num(el["P_ANALYTICAL"]).to_numpy(dtype=float))
    out = el.copy()
    for kind in ("ANALYTICAL", "EMPIRICAL", "TOPCLUSTER"):
        q = pd.Series(np.nan, index=el.index, dtype=float)
        if tested.any():
            q[tested] = multipletests(num(el.loc[tested, "P_" + kind]).to_numpy(dtype=float), alpha=0.01, method="fdr_bh")[1]
        out["Q_" + kind] = q
    keys = pd.DataFrame({"q": out["Q_ANALYTICAL"], "p": num(el["P_ANALYTICAL"]), "score": num(el["SCORE"]),
                         "cgc": el["CGC"], "id": el["ENSID"]})
    order = keys.sort_values(["q", "p", "score", "cgc", "id"], ascending=[True, True, False, False, True],
                             kind="mergesort", na_position="last").index
    out = out.loc[order]
    out.to_csv(a.out_elements, sep="\t", index=False, na_rep="")

    cl = pd.concat([read_text(p, CLUSTER_COLUMNS) for p in a.clusters], ignore_index=True)
    rank = {sym: i + 1 for i, sym in enumerate(out["SYMBOL"])}
    missing = sorted(set(cl["SYMBOL"]) - set(rank))
    if missing:
        sys.exit("ERROR: clusters of elements absent from the element tables: " + ", ".join(missing[:5]))
    cl["RANK"] = cl["SYMBOL"].map(rank).astype(int)
    ck = pd.DataFrame({"rank": cl["RANK"], "p": num(cl["P"]), "score": num(cl["SCORE"])})
    cl = cl.loc[ck.sort_values(["rank", "p", "score"], ascending=[True, True, False], kind="mergesort").index]
    cl.to_csv(a.out_clusters, sep="\t", index=False, na_rep="")
    print(f"merge_oncodriveclustl_shards: {len(a.elements)} shards, {len(out)} elements ({int(tested.sum())} with a p-value), "
          f"{len(cl)} clusters", file=sys.stderr)


if __name__ == "__main__":
    main()
