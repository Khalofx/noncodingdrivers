#!/usr/bin/env python3
"""Project a canonical mutation table into per-sample PASS VCFs.

Gate tooling only (the VCF-parity gate of the build guide): the pipeline itself reads VCFs
with nf-core modules and never runs this script.

Input: tab-separated with columns sample, chrom, pos, ref, alt (1-based, GRCh38, chr-named,
MAF allele conventions as in the mutations tier); other columns are ignored.

The canonical table writes indels the MAF way: a deletion starts at its first deleted base
with alt '-', an insertion sits on the base before it with ref '-'. VCF needs an anchor base,
so with --fasta each indel is written VCF-style:
  deletion  pos p, ref R, alt '-'  ->  POS p-1, REF = base(p-1) + R, ALT = base(p-1)
  insertion pos p, ref '-', alt A  ->  POS p,   REF = base(p),       ALT = base(p) + A
The deleted bases are compared with the FASTA; a disagreement stops the script. Without
--fasta, a table that contains indels is refused (the alleles would not be valid VCF).

Output: <outdir>/<sample>.vcf, records sorted by the contig order of the FASTA index and by
position, FILTER=PASS, GT 0/1. Contigs and order come from <fasta>.fai (or --fai). Rows on
contigs absent from the index stop the script, or are dropped and listed with
--drop-unknown-contigs.
"""

import argparse
import os
import sys

import pandas as pd


class Fasta:
    """Random access to an uncompressed, samtools-indexed FASTA (pure Python)."""

    def __init__(self, path, fai=None):
        self.fh = open(path, "rb")
        self.index = {}
        for line in open(fai or f"{path}.fai"):
            name, length, offset, line_bases, line_width = line.rstrip("\n").split("\t")[:5]
            self.index[name] = (int(length), int(offset), int(line_bases), int(line_width))

    def fetch(self, chrom, start, end):
        """Bases start..end (1-based, inclusive), upper case."""
        length, offset, line_bases, line_width = self.index[chrom]
        if start < 1 or end > length or start > end:
            sys.exit(f"{chrom}:{start}-{end} is outside the FASTA contig (length {length})")
        first = offset + (start - 1) // line_bases * line_width + (start - 1) % line_bases
        last = offset + (end - 1) // line_bases * line_width + (end - 1) % line_bases
        self.fh.seek(first)
        return self.fh.read(last - first + 1).replace(b"\n", b"").replace(b"\r", b"").decode().upper()


def to_vcf_alleles(fasta, chrom, pos, ref, alt):
    """Canonical (MAF-convention) alleles -> VCF POS, REF, ALT."""
    if ref != "-" and alt != "-":
        return pos, ref, alt
    if fasta is None:
        sys.exit("the table contains indels ('-' alleles): pass --fasta to anchor them")
    if alt == "-":  # deletion of ref starting at pos
        seq = fasta.fetch(chrom, pos - 1, pos + len(ref) - 1)
        if seq[1:] != ref.upper():
            sys.exit(f"deleted bases {chrom}:{pos} {ref} do not match the FASTA ({seq[1:]})")
        return pos - 1, seq, seq[0]
    anchor = fasta.fetch(chrom, pos, pos)  # insertion after pos
    return pos, anchor, anchor + alt


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("table")
    ap.add_argument("outdir")
    ap.add_argument("--fasta", help="uncompressed FASTA with .fai beside it: indel anchor bases, contigs and order")
    ap.add_argument("--fai", help="FASTA index for contig order and ##contig lines (default: <fasta>.fai)")
    ap.add_argument("--drop-unknown-contigs", action="store_true",
                    help="drop rows on contigs absent from the FASTA index (listed on stderr) instead of stopping")
    a = ap.parse_args(argv)

    d = pd.read_csv(a.table, sep="\t", dtype={"sample": str, "chrom": str, "pos": int, "ref": str, "alt": str})
    missing = {"sample", "chrom", "pos", "ref", "alt"} - set(d.columns)
    if missing:
        sys.exit(f"missing columns: {sorted(missing)}")
    fasta = Fasta(a.fasta, a.fai) if a.fasta else None
    fai = a.fai or (f"{a.fasta}.fai" if a.fasta else None)
    contigs = []
    if fai:
        contigs = [line.split("\t")[:2] for line in open(fai) if line.strip()]
        order = {c: i for i, (c, _) in enumerate(contigs)}
        unknown = set(d["chrom"]) - set(order)
        if unknown and not a.drop_unknown_contigs:
            sys.exit(f"contigs not in {fai}: {sorted(unknown)}")
        if unknown:
            off = d[d["chrom"].isin(unknown)]
            for chrom, n in off["chrom"].value_counts().sort_index().items():
                print(f"dropped {n} rows on {chrom} (not in {fai})", file=sys.stderr)
            d = d[~d["chrom"].isin(unknown)]
    else:
        order = {c: i for i, c in enumerate(dict.fromkeys(d["chrom"]))}
    d = d.drop_duplicates(["sample", "chrom", "pos", "ref", "alt"])
    vcf = [to_vcf_alleles(fasta, r.chrom, r.pos, r.ref, r.alt) for r in d.itertuples(index=False)]
    d = d.assign(vpos=[v[0] for v in vcf], vref=[v[1] for v in vcf], valt=[v[2] for v in vcf])
    d = d.assign(_o=d["chrom"].map(order)).sort_values(["sample", "_o", "vpos", "vref", "valt"])
    n_indel = int(((d["ref"] == "-") | (d["alt"] == "-")).sum())
    os.makedirs(a.outdir, exist_ok=True)
    for sample, g in d.groupby("sample", sort=True):
        with open(os.path.join(a.outdir, f"{sample}.vcf"), "w") as fh:
            fh.write("##fileformat=VCFv4.2\n")
            fh.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
            for name, length in contigs:
                fh.write(f"##contig=<ID={name},length={length}>\n")
            fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            fh.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n")
            for r in g.itertuples(index=False):
                fh.write(f"{r.chrom}\t{r.vpos}\t.\t{r.vref}\t{r.valt}\t.\tPASS\t.\tGT\t0/1\n")
    print(f"{d['sample'].nunique()} samples, {len(d)} records ({n_indel} indels anchored)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
