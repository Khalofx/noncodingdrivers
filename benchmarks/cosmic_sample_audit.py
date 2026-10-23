#!/usr/bin/env python3
"""Audit COSMIC exact matches for observations originating in the input cohort.

PCAWG mutation tables use aliquot UUIDs while COSMIC commonly uses submitter
sample names (for example RK004_C01). A literal string comparison therefore
misses self-cohort observations. This script crosswalks the exact-match output
through the PCAWG sample sheet and writes audited copies of all result tables.

The audit labels an observation ``same_cohort`` when its COSMIC sample name maps
to a tumour WGS aliquot from the requested project. Unmapped observations remain
``not_assessed`` rather than being called independent: establishing independence
also requires study and tissue provenance.
"""

import argparse
import csv
import json
import os
import sys
from collections import defaultdict
from pathlib import Path


def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader), list(reader.fieldnames or [])


def write_tsv(path, rows, fieldnames):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with open(temporary, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n",
            extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def _require(fieldnames, required, label):
    missing = set(required) - set(fieldnames)
    if missing:
        raise ValueError(f"{label} lacks required columns: {sorted(missing)}")


def load_mapper(path, project):
    """Return exact-alias and donor-prefix maps to tumour WGS aliquot UUIDs."""
    rows, fields = read_tsv(path)
    _require(fields, {"dcc_project_code", "aliquot_id"}, "sample map")
    exact = defaultdict(set)
    donors = defaultdict(set)
    exact_columns = (
        "aliquot_id", "submitter_specimen_id", "icgc_specimen_id",
        "submitter_sample_id", "icgc_sample_id",
    )
    donor_columns = ("submitter_donor_id", "icgc_donor_id")
    for row in rows:
        if row.get("dcc_project_code", "").strip() != project:
            continue
        strategy = row.get("library_strategy", "").strip().upper()
        if strategy and strategy != "WGS":
            continue
        specimen_type = row.get("dcc_specimen_type", "").strip().lower()
        if specimen_type and "normal" in specimen_type:
            continue
        aliquot = row["aliquot_id"].strip()
        if not aliquot:
            continue
        for column in exact_columns:
            value = row.get(column, "").strip()
            if value:
                exact[value].add(aliquot)
        for column in donor_columns:
            value = row.get(column, "").strip()
            if value:
                donors[value].add(aliquot)
    return {"exact": exact, "donor": donors}


def map_cosmic_sample(mapper, cosmic_sample_id, sample_name):
    matches = set()
    for value in (cosmic_sample_id, sample_name):
        if value:
            matches.update(mapper["exact"].get(value, ()))
    if sample_name:
        for donor, aliquots in mapper["donor"].items():
            if sample_name == donor or sample_name.startswith(donor + "_"):
                matches.update(aliquots)
    return matches


def variant_key(row):
    return row["chrom"], row["pos"], row["ref"], row["alt"]


def audit(exact_path, mutation_path, element_path, summary_path,
          sample_map_path, project, out_prefix):
    mapper = load_mapper(sample_map_path, project)
    exact_rows, exact_fields = read_tsv(exact_path)
    mutation_rows, mutation_fields = read_tsv(mutation_path)
    element_rows, element_fields = read_tsv(element_path)
    _require(exact_fields,
             {"chrom", "pos", "ref", "alt", "cosmic_sample_id",
              "cosmic_sample_name", "somatic_status", "element_ids"},
             "exact-match table")
    _require(mutation_fields, {"chrom", "pos", "ref", "alt", "cosmic_match"},
             "mutation-summary table")
    _require(element_fields, {"element_id", "cosmic_match"},
             "element-summary table")

    audited_exact = []
    observations_by_variant = defaultdict(list)
    external_samples_by_element = defaultdict(set)
    external_observations_by_element = defaultdict(int)
    for source in exact_rows:
        row = dict(source)
        mapped = map_cosmic_sample(
            mapper, row.get("cosmic_sample_id", ""), row.get("cosmic_sample_name", ""))
        row.pop("exact_sample_name_overlap", None)
        row["cohort_sample_mapping_overlap"] = "yes" if mapped else "no"
        row["mapped_cohort_samples"] = ",".join(sorted(mapped))
        row["independence_status"] = "same_cohort" if mapped else "not_assessed"
        audited_exact.append(row)
        observations_by_variant[variant_key(row)].append(row)
        if not mapped:
            for element_id in filter(None, row.get("element_ids", "").split(",")):
                external_observations_by_element[element_id] += 1
                cosmic_sample = row.get("cosmic_sample_id", "")
                if cosmic_sample:
                    external_samples_by_element[element_id].add(cosmic_sample)

    audited_mutations = []
    for source in mutation_rows:
        row = dict(source)
        row.pop("n_exact_sample_name_overlaps", None)
        observations = observations_by_variant.get(variant_key(row), [])
        cohort_obs = [r for r in observations
                      if r["cohort_sample_mapping_overlap"] == "yes"]
        external_obs = [r for r in observations
                        if r["cohort_sample_mapping_overlap"] == "no"]
        cohort_samples = {r["cosmic_sample_id"] for r in cohort_obs
                          if r["cosmic_sample_id"]}
        external_samples = {r["cosmic_sample_id"] for r in external_obs
                            if r["cosmic_sample_id"]}
        row["n_cohort_mapped_cosmic_observations"] = len(cohort_obs)
        row["n_cohort_mapped_cosmic_samples"] = len(cohort_samples)
        row["n_external_cosmic_observations"] = len(external_obs)
        row["n_external_cosmic_samples"] = len(external_samples)
        if not observations:
            row["independence_status"] = "no_match"
        elif not external_obs:
            row["independence_status"] = "cohort_only"
        else:
            row["independence_status"] = "external_not_assessed"
        audited_mutations.append(row)

    audited_elements = []
    for source in element_rows:
        row = dict(source)
        element_id = row["element_id"]
        row["n_external_cosmic_observations"] = external_observations_by_element[element_id]
        row["n_external_cosmic_samples"] = len(external_samples_by_element[element_id])
        row["external_cosmic_observation"] = (
            "yes" if external_observations_by_element[element_id] else "no")
        audited_elements.append(row)

    exact_fields = [f for f in exact_fields if f != "exact_sample_name_overlap"]
    for field in ("cohort_sample_mapping_overlap", "mapped_cohort_samples",
                  "independence_status"):
        if field not in exact_fields:
            exact_fields.append(field)
    mutation_fields = [f for f in mutation_fields if f != "n_exact_sample_name_overlaps"]
    for field in ("n_cohort_mapped_cosmic_observations",
                  "n_cohort_mapped_cosmic_samples",
                  "n_external_cosmic_observations", "n_external_cosmic_samples",
                  "independence_status"):
        if field not in mutation_fields:
            mutation_fields.append(field)
    for field in ("n_external_cosmic_observations", "n_external_cosmic_samples",
                  "external_cosmic_observation"):
        if field not in element_fields:
            element_fields.append(field)

    write_tsv(f"{out_prefix}.exact_matches.tsv", audited_exact, exact_fields)
    write_tsv(f"{out_prefix}.mutation_summary.tsv", audited_mutations, mutation_fields)
    write_tsv(f"{out_prefix}.element_summary.tsv", audited_elements, element_fields)

    summary = json.loads(Path(summary_path).read_text())
    same_cohort_observations = sum(
        row["cohort_sample_mapping_overlap"] == "yes" for row in audited_exact)
    external_observations = len(audited_exact) - same_cohort_observations
    mutations_with_external = sum(
        int(row["n_external_cosmic_observations"]) > 0 for row in audited_mutations)
    mutations_cohort_only = sum(
        row["independence_status"] == "cohort_only" for row in audited_mutations)
    confirmed_external_mutations = sum(
        any(obs["cohort_sample_mapping_overlap"] == "no"
            and obs["somatic_status"] == "Confirmed somatic variant"
            for obs in observations_by_variant.get(variant_key(row), []))
        for row in audited_mutations)
    summary.update({
        "sample_identity_audit": "PCAWG submitter/specimen/donor crosswalk",
        "sample_map_file": str(Path(sample_map_path).resolve()),
        "sample_map_project": project,
        "cosmic_observations_mapping_to_input_cohort": same_cohort_observations,
        "cosmic_observations_not_mapping_to_input_cohort": external_observations,
        "cosmic_matched_mutations_cohort_only": mutations_cohort_only,
        "cosmic_matched_mutations_with_external_observation": mutations_with_external,
        "cosmic_confirmed_somatic_mutations_with_external_observation":
            confirmed_external_mutations,
        "interpretation": (
            "COSMIC presence records prior observation; it is not a driver label. "
            "Mapped LIRI-JP observations are same-cohort, while unmapped observations "
            "remain not_assessed for independence until study/tissue provenance is resolved."
        ),
    })
    summary_out = Path(f"{out_prefix}.summary.json")
    temporary = summary_out.with_name(summary_out.name + ".tmp")
    temporary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, summary_out)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exact-matches", required=True)
    parser.add_argument("--mutation-summary", required=True)
    parser.add_argument("--element-summary", required=True)
    parser.add_argument("--summary-json", required=True)
    parser.add_argument("--sample-map", required=True)
    parser.add_argument("--project", required=True)
    parser.add_argument("--out-prefix", required=True)
    args = parser.parse_args()
    summary = audit(
        args.exact_matches, args.mutation_summary, args.element_summary,
        args.summary_json, args.sample_map, args.project, args.out_prefix)
    sys.stderr.write(json.dumps(summary, sort_keys=True) + "\n")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, csv.Error, json.JSONDecodeError) as error:
        sys.exit(f"ERROR: {error}")
