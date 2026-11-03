#!/usr/bin/env python3
"""OncodriveCLUSTL prefilter (layout only; exact).

Keeps the catalogue rows of the element ids whose rows hold at least --min-snvs SNVs of the cohort
in total, in catalogue order, as BED4. OncodriveCLUSTL analyses only those ids (-emut, the same
value), and its run time grows with the size of the regions file, so the other rows are left out:
the analysed elements and their results do not change (docs/validation.md).

An SNV at 1-based position p lies in a BED row [start, end) when start < p <= end (bedtools
intersect -c of the SNV as the BED point [p - 1, p)); every SNV row counts, so recurrent positions
count once per sample, as OncodriveCLUSTL counts them. Rows with a non-numeric start or end count
no SNV (the regions adapter of the module skips them).
"""
import argparse
import sys

import numpy as np
import pandas as pd


def snv_positions(path):
    """Sorted 1-based SNV positions per contig of the canonical table (vtype SNP)."""
    m = pd.read_csv(path, sep="\t", usecols=["chrom", "pos", "vtype"], dtype={"chrom": str, "vtype": str},
                    keep_default_na=False)
    m = m[m["vtype"] == "SNP"]
    pos = pd.to_numeric(m["pos"], errors="raise").astype(np.int64)
    return {c: np.sort(g.to_numpy()) for c, g in pos.groupby(m["chrom"].to_numpy())}


def read_bed(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4 or line.startswith(("#", "track", "browser")):
                continue
            rows.append(f[:4])
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--mutations", required=True, help="canonical cohort table (header: sample chrom pos ref alt vtype ...)")
    ap.add_argument("--elements", required=True, help="element catalogue rows, headerless BED (chrom start end id ...)")
    ap.add_argument("--min-snvs", type=int, required=True, help="OncodriveCLUSTL -emut")
    ap.add_argument("--out", required=True, help="BED4 rows of the ids with >= --min-snvs SNVs")
    a = ap.parse_args(argv)

    pos = snv_positions(a.mutations)
    rows = read_bed(a.elements)
    counts = np.zeros(len(rows), dtype=np.int64)
    by_chrom = {}
    for i, r in enumerate(rows):
        if r[1].isdigit() and r[2].isdigit():
            by_chrom.setdefault(r[0], []).append(i)
    for chrom, idx in by_chrom.items():
        p = pos.get(chrom)
        if p is None:
            continue
        idx = np.asarray(idx)
        start = np.array([int(rows[i][1]) for i in idx], dtype=np.int64)
        end = np.array([int(rows[i][2]) for i in idx], dtype=np.int64)
        counts[idx] = np.maximum(np.searchsorted(p, end, side="right") - np.searchsorted(p, start, side="right"), 0)
    total = {}
    for r, n in zip(rows, counts):
        total[r[3]] = total.get(r[3], 0) + int(n)
    keep = {e for e, n in total.items() if n >= a.min_snvs}
    n_rows = 0
    with open(a.out, "w") as out:
        for r in rows:
            if r[3] in keep:
                out.write("\t".join(r) + "\n")
                n_rows += 1
    print(f"oncodriveclustl_prefilter: {len(keep)} of {len(total)} element ids hold >= {a.min_snvs} SNVs "
          f"({n_rows} of {len(rows)} rows kept)", file=sys.stderr)


if __name__ == "__main__":
    main()
