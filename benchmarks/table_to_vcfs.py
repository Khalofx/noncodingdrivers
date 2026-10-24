#!/usr/bin/env python3
"""Project a canonical mutation table into per-sample PASS VCFs.

Gate tooling only (the VCF-parity gate of the build guide): the pipeline itself reads VCFs
with nf-core modules and never runs this script.

Input: tab-separated with columns sample, chrom, pos, ref, alt (1-based, GRCh38, chr-named);
other columns are ignored. Output: <outdir>/<sample>.vcf, sorted by the contig order of the
FASTA index when --fai is given (otherwise by first appearance), FILTER=PASS on every record.
"""

import argparse
import os
import sys

import pandas as pd


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("table")
    ap.add_argument("outdir")
    ap.add_argument("--fai", help="FASTA index: contig order and ##contig lines")
    a = ap.parse_args(argv)

    d = pd.read_csv(a.table, sep="\t", dtype={"sample": str, "chrom": str, "pos": int, "ref": str, "alt": str})
    missing = {"sample", "chrom", "pos", "ref", "alt"} - set(d.columns)
    if missing:
        sys.exit(f"missing columns: {sorted(missing)}")
    contigs = []
    if a.fai:
        contigs = [line.split("\t")[:2] for line in open(a.fai) if line.strip()]
        order = {c: i for i, (c, _) in enumerate(contigs)}
        unknown = set(d["chrom"]) - set(order)
        if unknown:
            sys.exit(f"contigs not in {a.fai}: {sorted(unknown)}")
    else:
        order = {c: i for i, c in enumerate(dict.fromkeys(d["chrom"]))}
    d = d.assign(_o=d["chrom"].map(order)).sort_values(["sample", "_o", "pos", "ref", "alt"])
    d = d.drop_duplicates(["sample", "chrom", "pos", "ref", "alt"])
    os.makedirs(a.outdir, exist_ok=True)
    for sample, g in d.groupby("sample", sort=True):
        with open(os.path.join(a.outdir, f"{sample}.vcf"), "w") as fh:
            fh.write("##fileformat=VCFv4.2\n")
            fh.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
            for name, length in contigs:
                fh.write(f"##contig=<ID={name},length={length}>\n")
            fh.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n")
            for r in g.itertuples(index=False):
                fh.write(f"{r.chrom}\t{r.pos}\t.\t{r.ref}\t{r.alt}\t.\tPASS\t.\tGT\t0/1\n")
    print(f"{d['sample'].nunique()} samples, {len(d)} records", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
