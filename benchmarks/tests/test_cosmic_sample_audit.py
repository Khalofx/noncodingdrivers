"""Tests for PCAWG-to-COSMIC sample identity auditing."""

import csv
import json
import os
import sys


sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import cosmic_sample_audit as A


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
    return path


def read_tsv(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fixture(tmp_path):
    sample_map = write_tsv(
        tmp_path / "samples.tsv",
        ["dcc_project_code", "aliquot_id", "submitter_specimen_id",
         "icgc_specimen_id", "submitter_sample_id", "icgc_sample_id",
         "submitter_donor_id", "icgc_donor_id", "dcc_specimen_type",
         "library_strategy"],
        [
            ["LIRI-JP", "UUID-T", "RK001_C", "SP1", "RK001_C", "SA1",
             "RK001", "DO1", "Primary tumour - solid tissue", "WGS"],
            ["LIRI-JP", "UUID-N", "RK001_B", "SP2", "RK001_B", "SA2",
             "RK001", "DO1", "Normal - blood derived", "WGS"],
            ["OTHER", "UUID-O", "RK999_C", "SP9", "RK999_C", "SA9",
             "RK999", "DO9", "Primary tumour - solid tissue", "WGS"],
        ])
    exact_header = [
        "cohort", "chrom", "pos", "ref", "alt", "element_ids",
        "cosmic_sample_id", "cosmic_sample_name", "somatic_status",
        "exact_sample_name_overlap", "independence_status",
    ]
    exact = write_tsv(
        tmp_path / "exact.tsv", exact_header,
        [
            ["LIRI-JP", "chr1", 10, "A", "T", "E1", "COSS1", "RK001_C01",
             "Confirmed somatic variant", "no", "not_assessed"],
            ["LIRI-JP", "chr1", 10, "A", "T", "E1", "COSS2", "EXT-1",
             "Confirmed somatic variant", "no", "not_assessed"],
            ["LIRI-JP", "chr2", 20, "C", "G", "E2", "COSS3", "RK001_C",
             "Variant of unknown origin", "no", "not_assessed"],
        ])
    mutation_header = [
        "cohort", "chrom", "pos", "ref", "alt", "cosmic_match",
        "n_exact_sample_name_overlaps", "independence_status",
    ]
    mutations = write_tsv(
        tmp_path / "mut.tsv", mutation_header,
        [
            ["LIRI-JP", "chr1", 10, "A", "T", "yes", 0, "not_assessed"],
            ["LIRI-JP", "chr2", 20, "C", "G", "yes", 0, "not_assessed"],
            ["LIRI-JP", "chr3", 30, "G", "A", "no", 0, "no_match"],
        ])
    elements = write_tsv(
        tmp_path / "elements.tsv", ["element_id", "cosmic_match"],
        [["E1", "yes"], ["E2", "yes"], ["E3", "no"]])
    summary = tmp_path / "summary.json"
    summary.write_text(json.dumps({"cosmic_matched_distinct_mutations": 2}))
    return sample_map, exact, mutations, elements, summary


def test_crosswalk_marks_same_cohort_and_excludes_normal_alias(tmp_path):
    sample_map, exact, mutations, elements, summary = fixture(tmp_path)
    out_prefix = tmp_path / "audited"
    result = A.audit(exact, mutations, elements, summary, sample_map,
                     "LIRI-JP", out_prefix)
    rows = read_tsv(tmp_path / "audited.exact_matches.tsv")
    by_name = {row["cosmic_sample_name"]: row for row in rows}
    assert by_name["RK001_C01"]["cohort_sample_mapping_overlap"] == "yes"
    assert by_name["RK001_C01"]["mapped_cohort_samples"] == "UUID-T"
    assert by_name["RK001_C"]["independence_status"] == "same_cohort"
    assert by_name["EXT-1"]["independence_status"] == "not_assessed"
    assert result["cosmic_observations_mapping_to_input_cohort"] == 2
    assert result["cosmic_observations_not_mapping_to_input_cohort"] == 1


def test_mutation_and_element_summaries_separate_external_observations(tmp_path):
    sample_map, exact, mutations, elements, summary = fixture(tmp_path)
    out_prefix = tmp_path / "audited"
    result = A.audit(exact, mutations, elements, summary, sample_map,
                     "LIRI-JP", out_prefix)
    rows = {row["pos"]: row for row in
            read_tsv(tmp_path / "audited.mutation_summary.tsv")}
    assert rows["10"]["independence_status"] == "external_not_assessed"
    assert rows["10"]["n_cohort_mapped_cosmic_observations"] == "1"
    assert rows["10"]["n_external_cosmic_observations"] == "1"
    assert rows["20"]["independence_status"] == "cohort_only"
    assert rows["30"]["independence_status"] == "no_match"
    element_rows = {row["element_id"]: row for row in
                    read_tsv(tmp_path / "audited.element_summary.tsv")}
    assert element_rows["E1"]["external_cosmic_observation"] == "yes"
    assert element_rows["E2"]["external_cosmic_observation"] == "no"
    assert result["cosmic_matched_mutations_with_external_observation"] == 1
    assert result["cosmic_matched_mutations_cohort_only"] == 1
