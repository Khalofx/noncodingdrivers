#!/usr/bin/env python3
"""Build the PCAWG-style native-GRCh38 element catalogue (the pipeline's default element_catalogue).

Gene-anchored layer from GENCODE (comprehensive GTF) following the recipe of Rheinbay et al. 2020, Supplementary
Information section 4, and its subtraction order, with these departures:
  * lncRNAs: only filters 1-5 are applied (1-2 as one test: no protein-coding gene within 5 kb on the sense strand,
    which includes sense overlap; 3: spliced length > 200; 4: GENCODE biotype lncRNA; 5: no exon in a UCSC chainSelf
    chain with normScore >= 90, target span, self-diagonal removed). Rheinbay also required filter 6 (> 20% of bases
    in conserved elements) OR filter 7 (expressed in > 10% of PCAWG RNA-seq samples); that gate is NOT applied (not
    reproducible here), and neither MiTranscriptome nor the literature supplement is used, so this lncRNA set is
    broader than Rheinbay's;
  * small RNAs from GENCODE biotypes; with --trna, tRNA genes are added from the GENCODE predicted-tRNA file
    (tRNAscan-SE; primary chromosomes; Pseudo_tRNA excluded) as class trna. No separate Rfam or snoRNAdb import;
  * miRNA precursors from GENCODE biotype miRNA (class mirna); with --mirbase, mature miRNAs are added from miRBase
    (GFF3 type miRNA, GRCh38), collapsed per mature accession (MIMAT) across loci, as class mirnaMature. The set of
    potential novel miRNAs used by Rheinbay is not included.
  * trna and mirnaMature have the global CDS, pc_SS, 5UTR, 3UTR and promCore bases subtracted, as for smallrna.
Without --trna and --mirbase these two classes are omitted. The recipe and the expected MD5 are in docs/references.md.
Regulatory layer from the current Regulatory Build catalogue with promoters dropped and enhancers trimmed to the
ENCODE cCREs inside them.

Output: headerless BED (chrom start end element_id), sorted, 24 primary chromosomes; element ids
<class>__<gene_name>__<gene_id> for the gene-anchored layer and the Regulatory Build ids for the regulatory layer.
Also writes <out>.summary.tsv (ids, fragments, bases per class) and <out>.lncrna_filters.tsv.

usage: build_pcawg_style_catalogue.py --gtf G --catalogue C --ccre CCRE --chainself CS --out OUT
           [--chain-min-norm 90] [--trna gencode.vNN.tRNAs.gtf.gz] [--mirbase hsa.gff3] [--bedtools bedtools]
"""
import argparse
import gzip
import os
import re
import subprocess
import sys

import pandas as pd

BT = "bedtools"  # set from --bedtools in main()
PRIMARY = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
SMALL_RNA = {"snRNA", "snoRNA", "misc_RNA", "rRNA", "scaRNA", "ribozyme"}
ATTR = re.compile(r'(\w+) "([^"]*)"')


def read_gtf(path):
    """exons, CDS, UTR and transcripts of the primary chromosomes; 0-based half-open."""
    rows = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] not in PRIMARY or f[2] not in ("gene", "transcript", "exon", "CDS", "UTR"):
                continue
            a = dict(ATTR.findall(f[8]))
            rows.append((f[0], int(f[3]) - 1, int(f[4]), f[6], f[2], a.get("gene_id"), a.get("gene_name"), a.get("gene_type"),
                         a.get("transcript_id"), a.get("transcript_type")))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "strand", "feature", "gene_id", "gene_name", "gene_type",
                                       "transcript_id", "transcript_type"])


def read_trna(path):
    """GENCODE predicted tRNA genes (feature tRNA); primary chromosomes; Pseudo_tRNA excluded; 0-based half-open."""
    rows = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] not in PRIMARY or f[2] != "tRNA":
                continue
            a = dict(ATTR.findall(f[8]))
            if a.get("gene_type") == "Pseudo_tRNA":
                continue
            rows.append((f[0], int(f[3]) - 1, int(f[4]), a.get("gene_type"), a.get("gene_id")))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "gene_name", "gene_id"])


def read_mirbase_mature(path):
    """miRBase GFF3 mature miRNAs (type miRNA); primary chromosomes; gene_name = Name, gene_id = MIMAT accession."""
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[0] not in PRIMARY or f[2] != "miRNA":
                continue
            a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            rows.append((f[0], int(f[3]) - 1, int(f[4]), a.get("Name"), a.get("Alias", a.get("ID"))))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "gene_name", "gene_id"])


def safe(s):
    return re.sub(r"[^A-Za-z0-9._-]", "_", str(s))


def write_bed(df, path, cols=("chrom", "start", "end", "name")):
    df = df[df.end > df.start]
    df.loc[:, list(cols)].to_csv(path, sep="\t", header=False, index=False)


def bt(*args, out):
    with open(out, "w") as fh:
        subprocess.run([BT, *args], stdout=fh, check=True)


def merge_per_gene(df, cls, tmp):
    """Collapse bases per gene: one name per gene, intervals merged within the name."""
    if df.empty:
        return pd.DataFrame(columns=["chrom", "start", "end", "name"])
    d = df.assign(name=cls + "__" + df.gene_name.map(safe) + "__" + df.gene_id.map(safe))
    d = d[d.end > d.start].sort_values(["name", "chrom", "start"])
    out = []
    for (name, chrom), g in d.groupby(["name", "chrom"], sort=False):
        s0 = e0 = None
        for s, e in zip(g.start.values, g.end.values):
            if s0 is None:
                s0, e0 = s, e
            elif s <= e0:
                e0 = max(e0, e)
            else:
                out.append((chrom, s0, e0, name))
                s0, e0 = s, e
        out.append((chrom, s0, e0, name))
    return pd.DataFrame(out, columns=["chrom", "start", "end", "name"])


def subtract(el, masks, tmp, tag):
    """el minus the union of masks (each a DataFrame chrom start end ...)."""
    if el.empty or not masks:
        return el
    a, b, o = f"{tmp}/{tag}.a.bed", f"{tmp}/{tag}.b.bed", f"{tmp}/{tag}.o.bed"
    write_bed(el, a)
    m = pd.concat([x[["chrom", "start", "end"]] for x in masks if not x.empty])
    m = m.assign(name=".")
    write_bed(m, b)
    bt("subtract", "-a", a, "-b", b, out=o)
    r = pd.read_csv(o, sep="\t", header=None, names=["chrom", "start", "end", "name"]) if os.path.getsize(o) else el.iloc[0:0]
    return r


def introns(ex):
    """splice-site windows: 6 bp at donors, 20 bp at acceptors (intronic), per transcript."""
    rows = []
    for (tid, strand), g in ex.sort_values("start").groupby(["transcript_id", "strand"], sort=False):
        s, e = g.start.values, g.end.values
        gid, gname, chrom = g.gene_id.iloc[0], g.gene_name.iloc[0], g.chrom.iloc[0]
        for i in range(len(s) - 1):
            i0, i1 = e[i], s[i + 1]          # intron [i0, i1)
            if i1 <= i0:
                continue
            if strand == "+":
                rows += [(chrom, i0, min(i0 + 6, i1), gid, gname), (chrom, max(i1 - 20, i0), i1, gid, gname)]
            else:
                rows += [(chrom, max(i1 - 6, i0), i1, gid, gname), (chrom, i0, min(i0 + 20, i1), gid, gname)]
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "gene_id", "gene_name"])


def tss_windows(tx, flank=200):
    t = tx.copy()
    tss = t.start.where(t.strand == "+", t.end - 1)
    t["start"] = (tss - flank).clip(lower=0)
    t["end"] = tss + flank + 1
    return t


def main():
    ap = argparse.ArgumentParser(description="Build the PCAWG-style GRCh38 element catalogue")
    ap.add_argument("--gtf", required=True, help="GENCODE comprehensive annotation GTF (gz)")
    ap.add_argument("--catalogue", required=True, help="Regulatory Build element BED (chrom start end id)")
    ap.add_argument("--ccre", required=True, help="SCREEN cCRE BED (gz)")
    ap.add_argument("--chainself", required=True, help="UCSC chainSelf.txt.gz")
    ap.add_argument("--chain-min-norm", type=float, default=90.0)
    ap.add_argument("--out", required=True)
    ap.add_argument("--trna", default=None, help="GENCODE predicted tRNA GTF (adds class trna)")
    ap.add_argument("--mirbase", default=None, help="miRBase GFF3 on GRCh38 (adds class mirnaMature)")
    ap.add_argument("--bedtools", default="bedtools", help="bedtools executable")
    a = ap.parse_args()
    global BT
    BT = a.bedtools
    tmp = a.out + ".tmp"
    os.makedirs(tmp, exist_ok=True)
    g = read_gtf(a.gtf)
    print("gtf rows", len(g), file=sys.stderr)
    genes = g[g.feature == "gene"]
    tx = g[g.feature == "transcript"]
    ex = g[g.feature == "exon"]
    cds = g[g.feature == "CDS"]
    # --- protein-coding elements -----------------------------------------------------------------------
    pc_tx = tx[(tx.transcript_type == "protein_coding") & (tx.gene_type == "protein_coding")]
    pc_ex = ex[ex.transcript_id.isin(pc_tx.transcript_id)]
    utr = g[(g.feature == "UTR") & g.transcript_id.isin(pc_tx.transcript_id)]
    cspan = cds[cds.transcript_id.isin(pc_tx.transcript_id)].groupby("transcript_id").agg(cs=("start", "min"), ce=("end", "max"))
    utr = utr.join(cspan, on="transcript_id")
    five = utr[((utr.strand == "+") & (utr.end <= utr.cs)) | ((utr.strand == "-") & (utr.start >= utr.ce))]
    three = utr[((utr.strand == "+") & (utr.start >= utr.ce)) | ((utr.strand == "-") & (utr.end <= utr.cs))]
    CDS = cds[["chrom", "start", "end"]]
    pcSS = subtract(merge_per_gene(introns(pc_ex), "ss", tmp), [CDS], tmp, "ss")
    u5 = subtract(merge_per_gene(five, "5utr", tmp), [CDS, pcSS], tmp, "u5")
    u3 = subtract(merge_per_gene(three, "3utr", tmp), [CDS, pcSS, u5], tmp, "u3")
    prom = subtract(merge_per_gene(tss_windows(pc_tx), "promCore", tmp), [CDS, pcSS], tmp, "prom")
    # --- lncRNA selection: filters 1-5 of SI section 4 only; the required 6-or-7 gate (conservation or PCAWG
    # expression) is not applied, so more lncRNAs are kept than in Rheinbay's set ----------------------------
    lt = tx[(tx.gene_type == "lncRNA") & (tx.transcript_type == "lncRNA")].copy()
    lex = ex[ex.transcript_id.isin(lt.transcript_id)]
    slen = (lex.end - lex.start).groupby(lex.transcript_id).sum()
    lt["spliced"] = lt.transcript_id.map(slen).fillna(0)
    pcg = genes[genes.gene_type == "protein_coding"]
    write_bed(lt.assign(name=lt.transcript_id, score=0)[["chrom", "start", "end", "name", "score", "strand"]], f"{tmp}/lt.bed",
              cols=("chrom", "start", "end", "name", "score", "strand"))
    write_bed(pcg.assign(name=pcg.gene_id, score=0)[["chrom", "start", "end", "name", "score", "strand"]], f"{tmp}/pcg.bed",
              cols=("chrom", "start", "end", "name", "score", "strand"))
    bt("window", "-a", f"{tmp}/lt.bed", "-b", f"{tmp}/pcg.bed", "-w", "5000", "-sm", "-u", out=f"{tmp}/lt_near.bed")
    near = set(pd.read_csv(f"{tmp}/lt_near.bed", sep="\t", header=None)[3]) if os.path.getsize(f"{tmp}/lt_near.bed") else set()
    # self-chain: UCSC chainSelf (bin score tName tSize tStart tEnd qName qSize qStrand qStart qEnd id normScore)
    cs = pd.read_csv(a.chainself, sep="\t", header=None, usecols=[2, 4, 5, 6, 9, 10, 12],
                     names=["tName", "tStart", "tEnd", "qName", "qStart", "qEnd", "normScore"])
    cs = cs[(cs.normScore >= a.chain_min_norm) & cs.tName.isin(PRIMARY)]
    diag = (cs.tName == cs.qName) & (cs.tStart < cs.qEnd) & (cs.qStart < cs.tEnd)
    cs = cs[~diag]
    write_bed(cs.rename(columns={"tName": "chrom", "tStart": "start", "tEnd": "end"}).assign(name="."), f"{tmp}/selfchain.bed")
    write_bed(lex.assign(name=lex.transcript_id), f"{tmp}/lex.bed")
    bt("intersect", "-a", f"{tmp}/lex.bed", "-b", f"{tmp}/selfchain.bed", "-u", out=f"{tmp}/lex_self.bed")
    selfc = set(pd.read_csv(f"{tmp}/lex_self.bed", sep="\t", header=None)[3]) if os.path.getsize(f"{tmp}/lex_self.bed") else set()
    lt["f12_pc_within_5kb_sense"] = lt.transcript_id.isin(near)
    lt["f3_short"] = lt.spliced <= 200
    lt["f5_selfchain"] = lt.transcript_id.isin(selfc)
    keep = lt[~(lt.f12_pc_within_5kb_sense | lt.f3_short | lt.f5_selfchain)]
    lt[["transcript_id", "gene_id", "gene_name", "spliced", "f12_pc_within_5kb_sense", "f3_short", "f5_selfchain"]].to_csv(
        a.out + ".lncrna_filters.tsv", sep="\t", index=False)
    kex = lex[lex.transcript_id.isin(keep.transcript_id)]
    lncSS = subtract(merge_per_gene(introns(kex), "lncrnaSS", tmp), [CDS, pcSS, u5, u3, prom], tmp, "lss")
    lnc = subtract(merge_per_gene(kex, "lncrna", tmp), [CDS, pcSS, u5, u3, prom, lncSS], tmp, "lnc")
    lprom = subtract(merge_per_gene(tss_windows(keep), "lncrnaPromCore", tmp), [CDS, pcSS, u5, u3, prom, lncSS], tmp, "lprom")
    # --- short RNAs and miRNAs ---------------------------------------------------------------------------
    sm = subtract(merge_per_gene(ex[ex.gene_type.isin(SMALL_RNA)], "smallrna", tmp), [CDS, pcSS, u5, u3, prom], tmp, "sm")
    mi = subtract(merge_per_gene(ex[ex.gene_type == "miRNA"], "mirna", tmp), [CDS, pcSS, u5, u3, prom], tmp, "mi")
    extra = {}
    if a.trna:
        extra["trna"] = subtract(merge_per_gene(read_trna(a.trna), "trna", tmp), [CDS, pcSS, u5, u3, prom], tmp, "trna")
    if a.mirbase:
        extra["mirnaMature"] = subtract(merge_per_gene(read_mirbase_mature(a.mirbase), "mirnaMature", tmp),
                                        [CDS, pcSS, u5, u3, prom], tmp, "mm")
    # --- regulatory layer ------------------------------------------------------------------------------
    cat = pd.read_csv(a.catalogue, sep="\t", header=None, usecols=[0, 1, 2, 3], names=["chrom", "start", "end", "name"])
    cat["cls"] = cat.name.str.replace(r"_ENSR\d+$", "", regex=True)
    enh = cat[cat.cls == "enhancer"]
    keep_rb = cat[cat.cls.isin(["CTCF_binding_site", "open_chromatin_region", "TF_binding_site"])]
    write_bed(enh, f"{tmp}/enh.bed")
    cc = pd.read_csv(a.ccre, sep="\t", header=None, usecols=[0, 1, 2]).rename(columns={0: "chrom", 1: "start", 2: "end"}).assign(name=".")
    write_bed(cc, f"{tmp}/ccre.bed")
    bt("intersect", "-a", f"{tmp}/enh.bed", "-b", f"{tmp}/ccre.bed", out=f"{tmp}/enh_trim_raw.bed")
    et = pd.read_csv(f"{tmp}/enh_trim_raw.bed", sep="\t", header=None, names=["chrom", "start", "end", "name"])
    et = et.sort_values(["name", "chrom", "start"])
    # merge overlapping pieces within an enhancer id
    rows = []
    for (name, chrom), gg in et.groupby(["name", "chrom"], sort=False):
        s0 = e0 = None
        for s, e in zip(gg.start.values, gg.end.values):
            if s0 is None:
                s0, e0 = s, e
            elif s <= e0:
                e0 = max(e0, e)
            else:
                rows.append((chrom, s0, e0, name))
                s0, e0 = s, e
        rows.append((chrom, s0, e0, name))
    et = pd.DataFrame(rows, columns=["chrom", "start", "end", "name"])
    # --- assemble ------------------------------------------------------------------------------------------
    layers = {"ss": pcSS, "5utr": u5, "3utr": u3, "promCore": prom, "lncrnaSS": lncSS, "lncrna": lnc, "lncrnaPromCore": lprom,
              "smallrna": sm, "mirna": mi, "enhancer": et, "CTCF_binding_site": keep_rb[keep_rb.cls == "CTCF_binding_site"],
              "open_chromatin_region": keep_rb[keep_rb.cls == "open_chromatin_region"], "TF_binding_site": keep_rb[keep_rb.cls == "TF_binding_site"]}
    layers.update(extra)
    allel = pd.concat([d[["chrom", "start", "end", "name"]] for d in layers.values()])
    allel = allel[allel.chrom.isin(PRIMARY) & (allel.end > allel.start)]
    allel["ci"] = allel.chrom.map({c: i for i, c in enumerate(PRIMARY)})
    allel = allel.sort_values(["ci", "start", "end", "name"]).drop(columns="ci")
    allel.to_csv(a.out, sep="\t", header=False, index=False)
    summ = [dict(cls=k, ids=d.name.nunique(), fragments=len(d), bases=int((d.end - d.start).sum()),
                 median_id_bases=int((d.end - d.start).groupby(d.name).sum().median()) if len(d) else 0) for k, d in layers.items()]
    summ.append(dict(cls="enhancer_before_trim", ids=enh.name.nunique(), fragments=len(enh), bases=int((enh.end - enh.start).sum()),
                     median_id_bases=int((enh.end - enh.start).groupby(enh.name).sum().median())))
    pd.DataFrame(summ).to_csv(a.out + ".summary.tsv", sep="\t", index=False)
    print(pd.DataFrame(summ).to_string(index=False))
    print(f"lncRNA transcripts {len(lt)}: within 5 kb of a sense protein-coding gene {int(lt.f12_pc_within_5kb_sense.sum())}, "
          f"<= 200 bases {int(lt.f3_short.sum())}, self-chain {int(lt.f5_selfchain.sum())}; kept {len(keep)} "
          f"({keep.gene_id.nunique()} genes); self-chain intervals used {len(cs)} (normScore >= {a.chain_min_norm})")


if __name__ == "__main__":
    main()
