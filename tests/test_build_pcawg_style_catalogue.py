"""Unit test of the element-catalogue recipe (assets/recipes/build_pcawg_style_catalogue.py).

Tiny synthetic inputs: one protein-coding gene, one lncRNA, one snRNA, one miRNA, GENCODE-style
predicted tRNAs, miRBase mature miRNAs and a three-element Regulatory Build catalogue. Checks the
classes, the Rheinbay subtraction order (no element overlaps CDS; promCore overlaps the 5'UTR; the
3'UTR does not overlap splice sites), the cCRE trimming of enhancers, the tRNA and mature-miRNA
rules, and that the two optional inputs only add rows. Needs bedtools on PATH.
"""
import gzip
import os
import shutil
import subprocess
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(ROOT, "assets", "recipes", "build_pcawg_style_catalogue.py")
BEDTOOLS = shutil.which("bedtools")

pytestmark = pytest.mark.skipif(BEDTOOLS is None, reason="bedtools not on PATH")


def _gz(path, lines):
    with gzip.open(path, "wt") as fh:
        fh.write("".join(line + "\n" for line in lines))
    return str(path)


def _plain(path, lines):
    with open(path, "w") as fh:
        fh.write("".join(line + "\n" for line in lines))
    return str(path)


def _gtf_row(chrom, feature, start, end, strand, attrs):
    a = " ".join(f'{k} "{v}";' for k, v in attrs.items())
    return "\t".join([chrom, "TEST", feature, str(start), str(end), ".", strand, ".", a])


def _inputs(tmp):
    pc = {"gene_id": "ENSG01.1", "gene_name": "PCG", "gene_type": "protein_coding"}
    pt = dict(pc, transcript_id="ENST01.1", transcript_type="protein_coding")
    ln = {"gene_id": "ENSG02.1", "gene_name": "LNC", "gene_type": "lncRNA"}
    lt = dict(ln, transcript_id="ENST02.1", transcript_type="lncRNA")
    sn = {"gene_id": "ENSG03.1", "gene_name": "SNR", "gene_type": "snRNA"}
    st = dict(sn, transcript_id="ENST03.1", transcript_type="snRNA")
    mi = {"gene_id": "ENSG04.1", "gene_name": "MIR1", "gene_type": "miRNA"}
    mt = dict(mi, transcript_id="ENST04.1", transcript_type="miRNA")
    gtf = [
        # protein-coding gene, + strand: exons 1001-1100 and 1201-1400; CDS 1051-1100, 1201-1300
        _gtf_row("chr1", "gene", 1001, 1400, "+", pc),
        _gtf_row("chr1", "transcript", 1001, 1400, "+", pt),
        _gtf_row("chr1", "exon", 1001, 1100, "+", pt),
        _gtf_row("chr1", "exon", 1201, 1400, "+", pt),
        _gtf_row("chr1", "CDS", 1051, 1100, "+", pt),
        _gtf_row("chr1", "CDS", 1201, 1300, "+", pt),
        _gtf_row("chr1", "UTR", 1001, 1050, "+", pt),
        _gtf_row("chr1", "UTR", 1301, 1400, "+", pt),
        # spliced lncRNA far from the protein-coding gene (350 bases spliced)
        _gtf_row("chr1", "gene", 60001, 60500, "+", ln),
        _gtf_row("chr1", "transcript", 60001, 60500, "+", lt),
        _gtf_row("chr1", "exon", 60001, 60150, "+", lt),
        _gtf_row("chr1", "exon", 60301, 60500, "+", lt),
        _gtf_row("chr1", "gene", 80001, 80150, "+", sn),
        _gtf_row("chr1", "transcript", 80001, 80150, "+", st),
        _gtf_row("chr1", "exon", 80001, 80150, "+", st),
        _gtf_row("chr1", "gene", 90001, 90080, "+", mi),
        _gtf_row("chr1", "transcript", 90001, 90080, "+", mt),
        _gtf_row("chr1", "exon", 90001, 90080, "+", mt),
    ]
    trna = [
        _gtf_row("chr1", "tRNA", 95001, 95072, "-", {"gene_id": "101", "gene_type": "Ala_tRNA", "gene_name": "101"}),
        _gtf_row("chr1", "tRNA", 96001, 96072, "-", {"gene_id": "102", "gene_type": "Pseudo_tRNA", "gene_name": "102"}),
        _gtf_row("GL000251.2", "tRNA", 101, 172, "+", {"gene_id": "103", "gene_type": "Gly_tRNA", "gene_name": "103"}),
    ]
    mirbase = [
        "##gff-version 3",
        "chr1\t.\tmiRNA_primary_transcript\t90001\t90080\t.\t+\t.\tID=MI0000001;Alias=MI0000001;Name=hsa-mir-t1",
        "chr1\t.\tmiRNA\t90011\t90032\t.\t+\t.\tID=MIMAT0000001;Alias=MIMAT0000001;Name=hsa-miR-t1-5p;Derives_from=MI0000001",
        "chr1\t.\tmiRNA\t90051\t90072\t.\t+\t.\tID=MIMAT0000002;Alias=MIMAT0000002;Name=hsa-miR-t1-3p;Derives_from=MI0000001",
        # the same mature accession from a second locus: one element, two fragments
        "chr2\t.\tmiRNA\t5001\t5022\t.\t+\t.\tID=MIMAT0000001;Alias=MIMAT0000001;Name=hsa-miR-t1-5p;Derives_from=MI0000009",
        # a mature miRNA inside the 3'UTR: removed by the subtraction
        "chr1\t.\tmiRNA\t1321\t1342\t.\t+\t.\tID=MIMAT0000003;Alias=MIMAT0000003;Name=hsa-miR-t3-5p;Derives_from=MI0000003",
    ]
    rb = ["chr1\t70000\t71000\tenhancer_ENSR00000000001", "chr1\t72000\t72400\tCTCF_binding_site_ENSR00000000002",
          "chr1\t73000\t73400\tpromoter_ENSR00000000003"]
    ccre = ["chr1\t70100\t70300", "chr1\t70600\t70700", "chr1\t72100\t72200"]
    chain = ["0\t1000\tchr1\t248956422\t500000\t600000\tchr3\t198295559\t+\t1000\t101000\t1\t95.0"]
    return dict(
        gtf=_gz(tmp / "a.gtf.gz", gtf), trna=_gz(tmp / "t.gtf.gz", trna), mirbase=_plain(tmp / "hsa.gff3", mirbase),
        catalogue=_plain(tmp / "rb.bed", rb), ccre=_gz(tmp / "ccre.bed.gz", ccre), chainself=_gz(tmp / "chain.txt.gz", chain))


def _build(tmp, inp, out, extra=True):
    cmd = [sys.executable, SCRIPT, "--gtf", inp["gtf"], "--catalogue", inp["catalogue"], "--ccre", inp["ccre"],
           "--chainself", inp["chainself"], "--bedtools", BEDTOOLS, "--out", str(tmp / out)]
    if extra:
        cmd += ["--trna", inp["trna"], "--mirbase", inp["mirbase"]]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    rows = [line.rstrip("\n").split("\t") for line in open(tmp / out)]
    return [(c, int(s), int(e), n) for c, s, e, n in rows]


def _cls(name):
    return name.split("__")[0] if "__" in name else name.rsplit("_ENSR", 1)[0]


def _overlaps(rows, chrom, s, e):
    return [r for r in rows if r[0] == chrom and r[1] < e and s < r[2]]


def test_classes_and_subtraction_order(tmp_path):
    rows = _build(tmp_path, _inputs(tmp_path), "cat.bed")
    classes = {_cls(r[3]) for r in rows}
    assert classes == {"promCore", "5utr", "3utr", "ss", "lncrna", "lncrnaPromCore", "lncrnaSS", "smallrna", "mirna",
                       "mirnaMature", "trna", "enhancer", "CTCF_binding_site"}
    # Regulatory Build promoters are not carried over
    assert not [r for r in rows if r[3].startswith("promoter_")]
    # no element overlaps coding sequence (0-based CDS [1050, 1100) and [1200, 1300))
    for s, e in ((1050, 1100), (1200, 1300)):
        assert not _overlaps(rows, "chr1", s, e)
    # promCore: TSS +/- 200 minus CDS and splice sites; it keeps the 5'UTR bases (Rheinbay order)
    prom = sorted((r[1], r[2]) for r in rows if r[3] == "promCore__PCG__ENSG01.1")
    assert prom == [(800, 1050), (1106, 1180)]
    assert [(r[1], r[2]) for r in rows if r[3] == "5utr__PCG__ENSG01.1"] == [(1000, 1050)]
    # splice sites: 6 bp at the donor, 20 bp at the acceptor; the 3'UTR does not touch them
    assert sorted((r[1], r[2]) for r in rows if r[3] == "ss__PCG__ENSG01.1") == [(1100, 1106), (1180, 1200)]
    assert [(r[1], r[2]) for r in rows if r[3] == "3utr__PCG__ENSG01.1"] == [(1300, 1400)]
    # enhancers are trimmed to the cCREs inside them; other regulatory classes are unchanged
    assert sorted((r[1], r[2]) for r in rows if r[3] == "enhancer_ENSR00000000001") == [(70100, 70300), (70600, 70700)]
    assert [(r[1], r[2]) for r in rows if r[3] == "CTCF_binding_site_ENSR00000000002"] == [(72000, 72400)]
    # the lncRNA passes filters 1-5 and gets exons, promoter and splice-site elements
    assert sorted((r[1], r[2]) for r in rows if r[3] == "lncrnaSS__LNC__ENSG02.1") == [(60150, 60156), (60280, 60300)]


def test_trna_and_mature_mirna_rules(tmp_path):
    rows = _build(tmp_path, _inputs(tmp_path), "cat.bed")
    trna = [r for r in rows if r[3].startswith("trna__")]
    assert trna == [("chr1", 95000, 95072, "trna__Ala_tRNA__101")]  # no Pseudo_tRNA, no non-primary contig
    mm = sorted(r for r in rows if r[3].startswith("mirnaMature__"))
    assert mm == [("chr1", 90010, 90032, "mirnaMature__hsa-miR-t1-5p__MIMAT0000001"),
                  ("chr1", 90050, 90072, "mirnaMature__hsa-miR-t1-3p__MIMAT0000002"),
                  ("chr2", 5000, 5022, "mirnaMature__hsa-miR-t1-5p__MIMAT0000001")]
    # the precursor element stays the GENCODE miRNA gene
    assert [(r[1], r[2]) for r in rows if r[3] == "mirna__MIR1__ENSG04.1"] == [(90000, 90080)]


def test_optional_inputs_only_add_rows(tmp_path):
    inp = _inputs(tmp_path)
    full = _build(tmp_path, inp, "full.bed")
    base = _build(tmp_path, inp, "base.bed", extra=False)
    assert [r for r in full if not r[3].startswith(("trna__", "mirnaMature__"))] == base
    assert not [r for r in base if r[3].startswith(("trna__", "mirnaMature__"))]
