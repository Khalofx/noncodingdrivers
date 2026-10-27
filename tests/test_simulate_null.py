"""bin/simulate_null.py (SIMULATENULL, D32): determinism and the properties of the null."""
import hashlib
import os
import subprocess
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(ROOT, "bin", "simulate_null.py")
COMP = str.maketrans("ACGT", "TGCA")


def _context_class(tri):
    rc = tri.translate(COMP)[::-1]
    return min(tri, rc)


@pytest.fixture(scope="module")
def cohort(tmp_path_factory):
    d = tmp_path_factory.mktemp("null")
    rng = np.random.default_rng(7)
    seqs = {"chrA": "".join(rng.choice(list("ACGT"), 300_000)), "chrB": "N" * 50 + "".join(rng.choice(list("ACGT"), 120_000))}
    fa, fai = d / "g.fa", d / "g.fa.fai"
    offset = 0
    with open(fa, "w") as fh, open(fai, "w") as fx:
        for name, s in seqs.items():
            head = f">{name}\n"
            fh.write(head)
            offset += len(head)
            lines = [s[i:i + 60] for i in range(0, len(s), 60)]
            fh.write("\n".join(lines) + "\n")
            fx.write(f"{name}\t{len(s)}\t{offset}\t60\t61\n")
            offset += sum(len(x) + 1 for x in lines)
    rows = ["sample\tchrom\tpos\tref\talt\tvtype"]
    for k in range(400):
        name = "chrA" if k % 3 else "chrB"
        s = seqs[name]
        p = int(rng.integers(100, len(s) - 100))
        while s[p - 2:p + 1].count("N"):
            p += 1
        ref = s[p - 1]
        alt = "ACGT".replace(ref, "")[k % 3]
        rows.append(f"S{k % 7}\t{name}\t{p}\t{ref}\t{alt}\tSNP")
    rows.append("S1\tchrA\t5000\tAC\t-\tDEL")                 # not carried
    rows.append("S2\tchrA\t6000\t" + ("A" if seqs["chrA"][5999] != "A" else "C") + "\tG\tSNP")   # REF mismatch (or G>G skipped as non-SNV)
    tab = d / "m.tsv"
    tab.write_text("\n".join(rows) + "\n")
    return d, seqs, tab, fa, fai


def _run(d, tab, fa, fai, seed, name):
    out, summ = d / f"{name}.tsv", d / f"{name}.summary.tsv"
    subprocess.run([sys.executable, SCRIPT, "--mutations", str(tab), "--fasta", str(fa), "--fai", str(fai),
                    "--seed", str(seed), "--window", "50000", "--output", str(out), "--summary", str(summ)], check=True)
    return out, dict(line.rstrip("\n").split("\t") for line in open(summ) if not line.startswith("metric"))


def test_same_seed_same_table(cohort):
    d, _, tab, fa, fai = cohort
    a, _ = _run(d, tab, fa, fai, 11, "a")
    b, _ = _run(d, tab, fa, fai, 11, "b")
    c, _ = _run(d, tab, fa, fai, 12, "c")
    md5 = lambda p: hashlib.md5(open(p, "rb").read()).hexdigest()   # noqa: E731
    assert md5(a) == md5(b)
    assert md5(a) != md5(c)


def test_null_keeps_context_window_and_patient(cohort):
    d, seqs, tab, fa, fai = cohort
    out, summary = _run(d, tab, fa, fai, 11, "props")
    obs = [line.rstrip("\n").split("\t") for line in open(tab)][1:]
    null = [line.rstrip("\n").split("\t") for line in open(out)][1:]
    assert int(summary["non_snv_not_carried"]) >= 1
    assert int(summary["snvs_out"]) == len(null) > 300
    obs_snv = [r for r in obs if r[5] == "SNP" and len(r[3]) == 1 and r[3] != r[4] and seqs[r[1]][int(r[2]) - 1] == r[3]]
    assert int(summary["ref_mismatch_dropped"]) == sum(1 for r in obs if r[5] == "SNP" and r[3] != r[4] and len(r[3]) == 1) - len(obs_snv)
    moved = 0
    for o, n in zip(obs_snv, null):          # same order; no duplicates expected at this density
        assert n[0] == o[0] and n[1] == o[1]                    # patient and contig kept
        p, q = int(o[2]), int(n[2])
        assert abs(p - q) <= 50000                               # window
        s = seqs[n[1]]
        assert s[q - 1] == n[3]                                  # REF matches the genome at the new site
        assert _context_class(s[p - 2:p + 1]) == _context_class(s[q - 2:q + 1])
        same = s[p - 2:p + 1] == s[q - 2:q + 1]
        assert (n[3], n[4]) == ((o[3], o[4]) if same else (o[3].translate(COMP), o[4].translate(COMP)))
        moved += p != q
    assert moved > 0.9 * len(obs_snv)
