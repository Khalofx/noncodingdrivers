#!/usr/bin/env python3
"""Score a SYNTH-NC run against its answer key (data/synth_nc*/SYNTH-NC.truth.json).

For each combination method (EBM: column ebm_q; ACAT: column acat_q) at q < --fdr:
  * recovered: planted true drivers (designs without "decoy") whose element id is significant;
  * missed: the planted true drivers that are not, as "<planted_id> <element_id> (<n> donors)";
  * decoy_significant: whether every planted decoy is significant;
  * non_planted_significant: significant elements that are not planted elements (the definition of the
    regression references);
  * overlapping_planted_significant: the subset of those whose span overlaps a planted element
    (elements can overlap, e.g. the TERT promoter core and the TERT 5'UTR).
With --reference, the counts are compared with the synth_nc block of a regression reference JSON and the
exit status is 1 on any difference.

usage: score_synth_nc.py --truth T --consensus C [--fdr 0.1] [--reference R] [--out O]
"""
import argparse
import json
import sys

import pandas as pd

METHODS = {"ebm": "ebm_q", "acat": "acat_q"}


def score(truth, cons, fdr=0.1):
    planted = truth["planted"]
    drivers = [p for p in planted if "decoy" not in p["design"]]
    decoys = [p for p in planted if "decoy" in p["design"]]
    ids = {p["element_id"] for p in planted}
    spans = [(p["chrom"], int(p["start0"]), int(p["end"])) for p in planted]

    def touches(row):
        return any(row.chrom == c and int(row.start) < e and s < int(row.end) for c, s, e in spans)

    out = {"donors": truth["donors"], "planted_drivers": len(drivers), "planted_decoys": len(decoys)}
    for m, col in METHODS.items():
        if col not in cons.columns:
            continue
        sig = cons[cons[col] < fdr]
        sig_ids = set(sig.element_id)
        rec = [p for p in drivers if p["element_id"] in sig_ids]
        other = sig[~sig.element_id.isin(ids)]
        near = other.apply(touches, axis=1) if len(other) else pd.Series([], dtype=bool)
        out[m] = {
            "recovered": len(rec),
            "decoy_significant": bool(decoys) and all(p["element_id"] in sig_ids for p in decoys),
            "non_planted_significant": int(len(other)),
            "overlapping_planted_significant": int(near.sum()) if len(other) else 0,
            "missed": [f'{p["planted_id"]} {p["element_id"]} ({p["donors_mutated"]} donors)'
                       for p in drivers if p["element_id"] not in sig_ids],
        }
    return out


def compare(res, ref):
    bad = []
    for m in METHODS:
        if m not in ref:
            continue
        for k in ("recovered", "decoy_significant", "non_planted_significant", "missed"):
            if k in ref[m] and res.get(m, {}).get(k) != ref[m][k]:
                bad.append(f"{m}.{k}: expected {ref[m][k]!r}, got {res.get(m, {}).get(k)!r}")
    return bad


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--truth", required=True)
    ap.add_argument("--consensus", required=True)
    ap.add_argument("--fdr", type=float, default=0.1)
    ap.add_argument("--reference", help="regression reference JSON with a synth_nc block")
    ap.add_argument("--out", help="write the result as JSON")
    a = ap.parse_args()
    truth = json.load(open(a.truth))
    cons = pd.read_csv(a.consensus, sep="\t", dtype={"element_id": str})
    res = score(truth, cons, a.fdr)
    print(json.dumps(res, indent=2))
    if a.out:
        json.dump(res, open(a.out, "w"), indent=2)
    if a.reference:
        bad = compare(res, json.load(open(a.reference))["synth_nc"])
        for b in bad:
            print("FAIL", b)
        print("PASS" if not bad else "FAIL")
        sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
