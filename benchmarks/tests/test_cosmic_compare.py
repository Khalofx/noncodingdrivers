"""Tests for exact GRCh38 mutation matching against COSMIC NCV."""

import csv
import gzip
import json
import subprocess
import sys
from pathlib import Path


BIN = Path(__file__).resolve().parents[1] / "cosmic_compare.py"

COSMIC_HEADER = [
    "GENE_SYMBOL", "COSMIC_GENE_ID", "TRANSCRIPT_ACCESSION", "MUTATION_NC_ID",
    "SAMPLE_NAME", "COSMIC_SAMPLE_ID", "COSMIC_PHENOTYPE_ID",
    "GENOMIC_MUTATION_ID", "LEGACY_MUTATION_ID", "ZYGOSITY", "CHROMOSOME",
    "GENOME_START", "GENOME_STOP", "GENOMIC_WT_ALLELE", "GENOMIC_MUT_ALLELE",
    "COSMIC_STUDY_ID", "PUBMED_PMID", "HGVSG", "MUTATION_SOMATIC_STATUS",
]


def write_tsv(path, header, rows, gz=False):
    opener = gzip.open if gz else open
    with opener(path, "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
    return path


def cosmic_row(gene, transcript, nc_id, sample_name, coss, coso, cosv, chrom,
               start, stop, ref, alt, study, pmid, status, legacy="", hgvsg=""):
    return [gene, "COSG1", transcript, nc_id, sample_name, coss, coso, cosv,
            legacy, "Heterozygous", chrom, start, stop, ref, alt, study, pmid,
            hgvsg, status]


def run_comparison(tmp_path):
    consensus = write_tsv(
        tmp_path / "consensus.tsv",
        ["consensus_rank", "element_id", "chrom", "start", "end", "acat_q",
         "vote_count", "is_consensus"],
        [
            [1, "E1", "chr1", 100, 110, 0.01, 2, "True"],
            [2, "E2", "chr1", 104, 108, 0.02, 1, "true"],
            [3, "NOT_SELECTED", "chr1", 200, 210, 0.5, 0, "False"],
        ])
    mutations = write_tsv(
        tmp_path / "mutations.tsv",
        ["sample", "chrom", "pos", "ref", "alt", "vtype"],
        [
            ["START_BOUNDARY", "chr1", 100, "A", "T", "SNP"],
            ["S1", "chr1", 101, "C", "G", "SNP"],
            ["S2", "chr1", 105, "G", "A", "SNP"],
            ["S3", "chr1", 110, "A", "-", "DEL"],
            ["END_PLUS_ONE", "chr1", 111, "T", "C", "SNP"],
            ["NOT_CANDIDATE", "chr1", 205, "A", "C", "SNP"],
        ])
    cosmic = write_tsv(
        tmp_path / "cosmic.tsv.gz", COSMIC_HEADER,
        [
            cosmic_row("G1", "ENST1", "COSN1", "S1", "COSS1", "COSO1",
                       "COSV1", "1", 101, 101, "C", "G", "COSU1", "111",
                       "Confirmed somatic variant", hgvsg="1:g.101C>G"),
            # Same genomic observation on a second transcript: one recurrence.
            cosmic_row("G2", "ENST2", "COSN2", "S1", "COSS1", "COSO1",
                       "COSV1", "1", 101, 101, "C", "G", "COSU1", "111",
                       "Confirmed somatic variant", hgvsg="1:g.101C>G"),
            # MAF '-' must match COSMIC's empty mutant allele.
            cosmic_row("", "", "COSN3", "OTHER", "COSS2", "COSO2",
                       "COSV2", "1", 110, 110, "A", "", "", "222",
                       "Variant of unknown origin", hgvsg="1:g.110del"),
            # Same position but different allele is not an exact match.
            cosmic_row("", "", "COSN4", "OTHER", "COSS3", "COSO2",
                       "COSV3", "1", 105, 105, "G", "T", "", "333",
                       "Confirmed somatic variant"),
        ], gz=True)
    outdir = tmp_path / "out"
    result = subprocess.run(
        [sys.executable, str(BIN), "--cosmic", str(cosmic),
         "--consensus", str(consensus), "--mutations", str(mutations),
         "--cohort", "TEST", "--outdir", str(outdir), "--progress-every", "0"],
        capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    return outdir


def read_tsv(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_bed_boundaries_and_candidate_selection(tmp_path):
    outdir = run_comparison(tmp_path)
    rows = read_tsv(outdir / "TEST.cosmic.mutation_summary.tsv")
    positions = {int(row["pos"]) for row in rows}
    assert positions == {101, 105, 110}
    assert 100 not in positions                 # BED start is excluded
    assert 110 in positions                     # BED end maps to the last 1-based base
    assert all("NOT_SELECTED" not in row["element_ids"] for row in rows)
    overlap = next(row for row in rows if row["pos"] == "105")
    assert set(overlap["element_ids"].split(",")) == {"E1", "E2"}


def test_exact_matching_and_unmatched_are_explicit(tmp_path):
    outdir = run_comparison(tmp_path)
    rows = {int(row["pos"]): row for row in
            read_tsv(outdir / "TEST.cosmic.mutation_summary.tsv")}
    assert rows[101]["cosmic_match"] == "yes"
    assert rows[101]["n_confirmed_somatic_samples"] == "1"
    assert rows[105]["cosmic_match"] == "no"     # wrong ALT in COSMIC
    assert rows[110]["cosmic_match"] == "yes"    # '-' == empty COSMIC allele
    assert rows[110]["somatic_statuses"] == "Variant of unknown origin"


def test_transcript_rows_are_deduplicated_but_annotations_preserved(tmp_path):
    outdir = run_comparison(tmp_path)
    exact = read_tsv(outdir / "TEST.cosmic.exact_matches.tsv")
    at_101 = [row for row in exact if row["pos"] == "101"]
    assert len(at_101) == 1
    assert at_101[0]["source_transcript_rows"] == "2"
    assert set(at_101[0]["gene_symbols"].split(",")) == {"G1", "G2"}
    assert at_101[0]["exact_sample_name_overlap"] == "yes"
    assert at_101[0]["independence_status"] == "possible_self_overlap"


def test_summary_totals_distinguish_rows_observations_and_variants(tmp_path):
    outdir = run_comparison(tmp_path)
    summary = json.loads((outdir / "TEST.cosmic.summary.json").read_text())
    assert summary["candidate_elements"] == 2
    assert summary["candidate_sample_events"] == 3
    assert summary["candidate_distinct_mutations"] == 3
    assert summary["cosmic_rows_scanned"] == 4
    assert summary["cosmic_raw_exact_rows"] == 3
    assert summary["cosmic_deduplicated_observations"] == 2
    assert summary["cosmic_matched_distinct_mutations"] == 2
    assert summary["cosmic_unmatched_distinct_mutations"] == 1
