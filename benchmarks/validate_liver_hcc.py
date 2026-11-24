#!/usr/bin/env python3
"""Liver-HCC-public acceptance checks for a noncodingdrivers results folder.

Layer "core" (driver core) reads consensus/<cohort>.consensus.tsv and checks the significant set
(is_consensus, which follows --combination_method and --fdr_threshold), its locus count, and the
TERT promoter core element of the PCAWG-style catalogue (promCore__TERT__ENSG00000164362.23):
significant, consensus rank and consensus q. Layer "full" adds the explanation-layer outputs: the
candidate-mutation table (only significant elements; 10 TERT rows from 10 patients) and the MultiQC
report.

Expected values are those of the PCAWG-style catalogue reference
(tests/regression/reference_pcawg_style.json, datasheet 1.3.0): under ACAT, 28 significant
annotations at 23 loci with TERT at rank 3 (q 2.36e-08); under EBM, 38 significant at 28 loci with
TERT at rank 5 (q 1.18e-07, recorded to three significant digits). On the Regulatory Build catalogue
(tests/regression/reference_phase2.json) the values were ACAT 23 at 19 loci and EBM 24 at 19 loci
(TERT promoter_ENSR00000746745, 11 rows); the frozen pre-v1 baseline (tests/regression/baseline.json)
had ACAT 28 at 23 loci and EBM 31, TERT rank 5.
"""

import argparse
import csv
import json
from pathlib import Path


TERT_ELEMENT = "promCore__TERT__ENSG00000164362.23"
TERT_ROWS = 10  # candidate rows (= distinct patients) of the TERT element in the full layer
EXPECTED = {
    "acat": {"significant": 28, "loci": 23, "tert_rank": 3, "tert_q": 2.36051422129596e-08, "q_tol": 1e-9},
    "ebm": {"significant": 38, "loci": 28, "tert_rank": 5, "tert_q": 1.18e-07, "q_tol": 5e-10},
}


def read_tsv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def truth(value):
    return str(value).strip().lower() in {"1", "true", "t", "yes"}


def locus_count(rows):
    """Number of loci: element intervals (half-open) merged where they overlap, per chromosome."""
    intervals = sorted((row["chrom"], int(float(row["start"])),
                        int(float(row["end"]))) for row in rows)
    count = 0
    current_chrom = None
    current_end = -1
    for chrom, start, end in intervals:
        if chrom != current_chrom or start >= current_end:
            count += 1
            current_chrom, current_end = chrom, end
        else:
            current_end = max(current_end, end)
    return count


def require(condition, message, failures):
    if not condition:
        failures.append(message)


def check_core(consensus, method, failures):
    exp = EXPECTED[method]
    significant = [row for row in consensus if truth(row.get("is_consensus"))]
    methods = {row.get("combination_method", "") for row in consensus} - {""}
    require(methods <= {method}, f"table was combined with {sorted(methods)}, not {method}", failures)
    require(len(significant) == exp["significant"],
            f"expected {exp['significant']} significant annotations, found {len(significant)}", failures)
    if exp["loci"] is not None:
        require(locus_count(significant) == exp["loci"],
                f"expected {exp['loci']} significant loci, found {locus_count(significant)}", failures)
    tert = [row for row in significant if row["element_id"] == TERT_ELEMENT]
    require(len(tert) == 1, f"expected one significant TERT element row, found {len(tert)}", failures)
    if len(tert) == 1:
        require(int(float(tert[0]["consensus_rank"])) == exp["tert_rank"],
                f"expected TERT rank {exp['tert_rank']}, found {tert[0]['consensus_rank']}", failures)
        qvalue = float(tert[0]["consensus_q"])
        require(abs(qvalue - exp["tert_q"]) <= exp["q_tol"],
                f"TERT {method} q outside tolerance: {qvalue}", failures)
    return significant


def check_full(root, cohort, significant, failures):
    variants_path = root / "candidates" / f"{cohort}.candidates.tsv"
    multiqc = root / "multiqc" / "multiqc_report.html"
    for path in (variants_path, multiqc):
        require(path.is_file() and path.stat().st_size > 0, f"missing or empty output: {path}", failures)
    if failures:
        return 0
    sig_ids = {row["element_id"] for row in significant}
    variants = read_tsv(variants_path)
    require(all(row.get("element_id") in sig_ids for row in variants),
            "candidate table contains a non-significant element", failures)
    tert_variants = [row for row in variants if row["element_id"] == TERT_ELEMENT]
    require(len(tert_variants) == TERT_ROWS, f"expected {TERT_ROWS} TERT mutation rows, found {len(tert_variants)}", failures)
    require(len({row["sample"] for row in tert_variants}) == TERT_ROWS,
            f"expected {TERT_ROWS} distinct TERT patients", failures)
    require(multiqc.stat().st_size > 100_000, "MultiQC report is unexpectedly small", failures)
    return len(variants)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results", required=True, type=Path, help="pipeline --outdir")
    parser.add_argument("--cohort", default="Liver-HCC-public")
    parser.add_argument("--method", default="ebm", choices=sorted(EXPECTED))
    parser.add_argument("--layer", default="core", choices=["core", "full"])
    args = parser.parse_args(argv)

    consensus_path = args.results / "consensus" / f"{args.cohort}.consensus.tsv"
    failures = []
    require(consensus_path.is_file() and consensus_path.stat().st_size > 0,
            f"missing or empty output: {consensus_path}", failures)
    if failures:
        raise SystemExit("ACCEPTANCE CHECK FAILED\n- " + "\n- ".join(failures))
    significant = check_core(read_tsv(consensus_path), args.method, failures)
    candidate_rows = check_full(args.results, args.cohort, significant, failures) if args.layer == "full" else None
    if failures:
        raise SystemExit("ACCEPTANCE CHECK FAILED\n- " + "\n- ".join(failures))
    print(json.dumps({
        "status": "pass", "layer": args.layer, "method": args.method,
        "significant_annotations": len(significant), "significant_loci": locus_count(significant),
        "candidate_rows": candidate_rows,
    }, indent=2))
    return 0


if __name__ == "__main__":
    main()
