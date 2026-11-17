#!/usr/bin/env python3
"""Remove hypermutated samples from a cohort mutation table, once, before every caller.

Rules (--rule):
  intogen  remove a sample with more than --min-mutations rows AND more than Q3 + --iqr-factor x IQR
           of the per-sample counts of its cohort (IntOGen preprocessing; 10,000 for WGS)
  cap      remove a sample with more than --min-mutations rows
  none     keep every sample

A sample's count is its number of rows in the table (all variant types). Quartiles use linear
interpolation between order statistics (numpy's default). Comparisons are strict.

Outputs (prefix P):
  P.mutations.tsv           the kept rows, header and input order preserved
  P.removed_mutations.tsv   the rows of removed samples (header only when none)
  P.sample_filter.tsv       one row per sample: count, kept, reason and the cohort thresholds
  P.sample_filter_mqc.tsv   MultiQC custom-content table (one row for the cohort)
"""
import argparse
import os
import sys
from collections import Counter


def quantile(sorted_values, p):
    """Linear interpolation between order statistics (numpy 'linear', R type 7)."""
    n = len(sorted_values)
    k = (n - 1) * p
    i = int(k)
    j = min(i + 1, n - 1)
    return sorted_values[i] + (sorted_values[j] - sorted_values[i]) * (k - i)


def decide(counts, rule, min_mutations, iqr_factor):
    """Return ({sample: reason} for removed samples, q1, q3, iqr_threshold)."""
    values = sorted(counts.values())
    q1, q3 = quantile(values, 0.25), quantile(values, 0.75)
    iqr_threshold = q3 + iqr_factor * (q3 - q1)
    removed = {}
    for sample, n in counts.items():
        if rule == "intogen" and n > min_mutations and n > iqr_threshold:
            removed[sample] = f"more than {min_mutations} mutations and above Q3 + {iqr_factor:g} x IQR ({iqr_threshold:.1f})"
        elif rule == "cap" and n > min_mutations:
            removed[sample] = f"more than {min_mutations} mutations"
    return removed, q1, q3, iqr_threshold


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="cohort mutation table with a 'sample' column (TSV, header)")
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--cohort", default=None, help="cohort name for the MultiQC table (default: prefix)")
    ap.add_argument("--rule", choices=["intogen", "cap", "none"], default="intogen")
    ap.add_argument("--min-mutations", type=int, default=10000)
    ap.add_argument("--iqr-factor", type=float, default=1.5)
    a = ap.parse_args(argv)
    if a.min_mutations < 0 or a.iqr_factor < 0:
        sys.exit("--min-mutations and --iqr-factor must be non-negative")
    cohort = a.cohort or os.path.basename(a.prefix)

    with open(a.input) as fh:
        header = fh.readline().rstrip("\n")
        cols = header.split("\t")
        if "sample" not in cols:
            sys.exit(f"{a.input}: no 'sample' column in header: {header}")
        si = cols.index("sample")
        counts = Counter(line.split("\t", si + 1)[si] for line in fh if line.strip())
    if not counts:
        sys.exit(f"{a.input}: no mutation rows")

    removed, q1, q3, thr = decide(counts, a.rule, a.min_mutations, a.iqr_factor)
    if len(removed) == len(counts):
        sys.exit(f"hypermutator filter ({a.rule}) would remove all {len(counts)} samples of {cohort}")

    n_kept_rows = n_removed_rows = 0
    with open(a.input) as fh, open(f"{a.prefix}.mutations.tsv", "w") as keep, open(f"{a.prefix}.removed_mutations.tsv", "w") as drop:
        fh.readline()
        keep.write(header + "\n"); drop.write(header + "\n")
        for line in fh:
            if not line.strip():
                continue
            if line.split("\t", si + 1)[si] in removed:
                drop.write(line); n_removed_rows += 1
            else:
                keep.write(line); n_kept_rows += 1

    with open(f"{a.prefix}.sample_filter.tsv", "w") as out:
        out.write("sample\tn_mutations\tkept\treason\trule\tmin_mutations\tiqr_factor\tq1\tq3\tiqr_threshold\n")
        for sample, n in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])):
            out.write(f"{sample}\t{n}\t{'false' if sample in removed else 'true'}\t{removed.get(sample, '')}\t"
                      f"{a.rule}\t{a.min_mutations}\t{a.iqr_factor:g}\t{q1:.1f}\t{q3:.1f}\t{thr:.1f}\n")

    rm_list = ",".join(sorted(removed)) if removed else "none"
    with open(f"{a.prefix}.sample_filter_mqc.tsv", "w") as out:
        out.write(f"# id: 'hypermutator_filter_{cohort}'\n"
                  f"# section_name: 'Hypermutator filter: {cohort}'\n"
                  f"# description: 'Samples removed once, before every caller (rule: {a.rule}). "
                  f"Per-sample counts and reasons: mutations/{a.prefix}.sample_filter.tsv.'\n"
                  "# plot_type: 'table'\n")
        out.write("cohort\tsamples\tsamples_removed\tmutations_kept\tmutations_removed\trule\tmin_mutations\tiqr_threshold\tremoved_samples\n")
        out.write(f"{cohort}\t{len(counts)}\t{len(removed)}\t{n_kept_rows}\t{n_removed_rows}\t{a.rule}\t{a.min_mutations}\t{thr:.1f}\t{rm_list}\n")

    print(f"{cohort}: {len(counts)} samples, {len(removed)} removed ({a.rule}); rows kept {n_kept_rows}, removed {n_removed_rows}; "
          f"Q1 {q1:.1f} Q3 {q3:.1f} threshold {thr:.1f}", file=sys.stderr)


if __name__ == "__main__":
    main()
