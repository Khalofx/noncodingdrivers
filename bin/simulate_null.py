#!/usr/bin/env python3
"""SIMULATENULL: a neutral copy of a cohort mutation table (in-house, D32).

Each SNV moves to a position drawn uniformly from the positions with the same trinucleotide
within +/- --window bases, in the same patient (the PCAWG Sanger null). Trinucleotides are
collapsed with their reverse complement (32 classes); when the drawn position carries the
reverse-complement context, REF and ALT are complemented, so the substitution class is kept.
The draw includes the original position. Non-SNV rows (MNVs, indels) are not carried into the
null table; SNVs whose REF differs from the FASTA, or whose context holds a non-ACGT base, are
dropped. All counts go to the log (stderr) and to --summary.

Deterministic: numpy PCG64 seeded with --seed; contigs in FASTA-index order, classes 0..31, SNVs
in table order within a class. Output rows keep the input order of the SNVs that were moved;
exact duplicates (same sample, position and alleles) are written once.
"""
import argparse
import sys

import numpy as np
import pandas as pd

COMP = str.maketrans("ACGT", "TGCA")
LUT = np.full(256, 4, dtype=np.uint8)
for i, b in enumerate(b"ACGT"):
    LUT[b] = i
    LUT[ord(chr(b).lower())] = i
# strand-collapsed trinucleotide classes: representative min(t, revcomp(t)) -> dense id 0..31
_T = np.arange(64)
_RC = 16 * (3 - _T % 4) + 4 * (3 - (_T // 4) % 4) + (3 - _T // 16)
_REPS = np.unique(np.minimum(_T, _RC))
DENSE = np.full(256, 255, dtype=np.uint8)
DENSE[_REPS] = np.arange(_REPS.size, dtype=np.uint8)


def read_fai(path):
    fai = []
    with open(path) as fh:
        for line in fh:
            name, length, offset, linebases, linewidth = line.split("\t")[:5]
            fai.append((name, int(length), int(offset), int(linebases), int(linewidth)))
    return fai


def read_contig(fasta, length, offset, linebases, linewidth):
    """Base codes 0-3 (A, C, G, T; 4 = other) of one contig, read through its .fai entry."""
    nlines = (length + linebases - 1) // linebases
    fasta.seek(offset)
    raw = fasta.read(nlines * linewidth)
    raw = raw + b"\n" * (nlines * linewidth - len(raw))
    seq = np.frombuffer(raw, dtype=np.uint8).reshape(nlines, linewidth)[:, :linebases].reshape(-1)[:length]
    return LUT[seq]


def context_classes(codes):
    """Class 0..31 of every position (255 where undefined) and whether its + strand context is the representative."""
    n = codes.size
    cls = np.full(n, 255, dtype=np.uint8)
    fwd = np.zeros(n, dtype=bool)
    if n < 3:
        return cls, fwd
    a, b, c = codes[:-2], codes[1:-1], codes[2:]
    ok = (a < 4) & (b < 4) & (c < 4)
    tri = (a << 4) + (b << 2) + c                                    # uint8; garbage where not ok
    rc = ((3 - c) << 4) + ((3 - b) << 2) + (3 - a)
    rep = np.minimum(tri, rc)
    del rc
    mid = DENSE[rep]
    mid[~ok] = 255
    cls[1:-1] = mid
    fwd[1:-1] = ok & (tri == rep)
    return cls, fwd


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mutations", required=True, help="canonical cohort table (sample chrom pos ref alt vtype)")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--fai", required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--window", type=int, default=50000)
    ap.add_argument("--output", required=True)
    ap.add_argument("--summary", required=True)
    args = ap.parse_args()

    tab = pd.read_csv(args.mutations, sep="\t", dtype=str, keep_default_na=False)
    missing = {"sample", "chrom", "pos", "ref", "alt", "vtype"} - set(tab.columns)
    if missing:
        sys.exit(f"SIMULATENULL: table lacks {sorted(missing)}")
    n_in = len(tab)
    acgt = tab["ref"].isin(list("ACGT")) & tab["alt"].isin(list("ACGT")) & (tab["ref"] != tab["alt"])
    snv = tab.loc[acgt, ["sample", "chrom", "pos", "ref", "alt"]].copy()
    snv["pos0"] = snv["pos"].astype(np.int64) - 1
    snv["order"] = np.arange(len(snv))
    counts = {"rows_in": n_in, "snvs_in": len(snv), "non_snv_not_carried": n_in - len(snv),
              "ref_mismatch_dropped": 0, "context_undefined_dropped": 0, "contig_absent_dropped": 0}

    rng = np.random.Generator(np.random.PCG64(args.seed))
    fai = read_fai(args.fai)
    known = {f[0] for f in fai}
    counts["contig_absent_dropped"] = int((~snv["chrom"].isin(known)).sum())
    parts = []
    with open(args.fasta, "rb") as fasta:
        for name, length, offset, linebases, linewidth in fai:
            sub = snv[snv["chrom"] == name]
            if sub.empty:
                continue
            codes = read_contig(fasta, length, offset, linebases, linewidth)
            cls, fwd = context_classes(codes)
            pos0 = sub["pos0"].to_numpy()
            inside = (pos0 >= 0) & (pos0 < length)
            p = np.clip(pos0, 0, length - 1)
            ref_code = LUT[np.frombuffer("".join(sub["ref"]).encode(), dtype=np.uint8)]
            ok_ref = inside & (codes[p] == ref_code)
            counts["ref_mismatch_dropped"] += int((~ok_ref).sum())
            mcls = np.where(ok_ref, cls[p], 255)
            counts["context_undefined_dropped"] += int((ok_ref & (mcls == 255)).sum())
            new_pos = np.full(pos0.size, -1, dtype=np.int64)
            flip = np.zeros(pos0.size, dtype=bool)
            for k in range(32):
                sel = np.flatnonzero(mcls == k)
                if sel.size == 0:
                    continue
                cand = np.flatnonzero(cls == k)
                lo = np.searchsorted(cand, pos0[sel] - args.window, side="left")
                hi = np.searchsorted(cand, pos0[sel] + args.window, side="right")
                pick = cand[lo + np.floor(rng.random(sel.size) * (hi - lo)).astype(np.int64)]
                new_pos[sel] = pick
                flip[sel] = fwd[pos0[sel]] != fwd[pick]
            keep = new_pos >= 0
            moved = sub.loc[keep, ["sample", "chrom", "ref", "alt", "order"]].copy()
            moved["pos"] = new_pos[keep] + 1
            fl = flip[keep]
            moved.loc[fl, "ref"] = moved.loc[fl, "ref"].str.translate(COMP)
            moved.loc[fl, "alt"] = moved.loc[fl, "alt"].str.translate(COMP)
            parts.append(moved)
            del codes, cls, fwd

    out = pd.concat(parts) if parts else pd.DataFrame(columns=["sample", "chrom", "ref", "alt", "order", "pos"])
    out = out.sort_values("order", kind="stable")
    n_moved = len(out)
    out = out.drop_duplicates(subset=["sample", "chrom", "pos", "ref", "alt"], keep="first")
    out["vtype"] = "SNP"
    out[["sample", "chrom", "pos", "ref", "alt", "vtype"]].to_csv(args.output, sep="\t", index=False)
    counts["duplicates_collapsed"] = n_moved - len(out)
    counts["snvs_out"] = len(out)
    with open(args.summary, "w") as fh:
        fh.write("metric\tvalue\n")
        for k, v in counts.items():
            fh.write(f"{k}\t{v}\n")
    print("SIMULATENULL " + " ".join(f"{k}={v}" for k, v in counts.items()), file=sys.stderr)


if __name__ == "__main__":
    main()
