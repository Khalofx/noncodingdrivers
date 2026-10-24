#!/usr/bin/env python3
"""Compare mutations in significant consensus elements with COSMIC NCV.

The comparison is deliberately variant-level. Candidate elements are BED intervals
(0-based, half-open), while the canonical mutation table and COSMIC coordinates are
1-based. Therefore a mutation position ``p`` belongs to ``[start, end)`` when
``start < p <= end``.

An exact COSMIC match requires the same GRCh38 chromosome, start position,
reference allele, and alternate allele. ``-`` in the canonical mutation table is
treated as COSMIC's empty allele for MAF-style insertions/deletions. No fuzzy or
interval-only matches are called exact.

COSMIC NCV is transcript-oriented, so one genomic observation may occupy several
rows. Recurrence is counted after collapsing transcript rows by genomic variant,
COSMIC sample, phenotype, study, PMID, and somatic status.

Outputs under --outdir:
  <cohort>.cosmic.mutation_summary.tsv  one row per distinct candidate variant
  <cohort>.cosmic.exact_matches.tsv     one row per deduplicated COSMIC observation
  <cohort>.cosmic.element_summary.tsv   one row per significant element
  <cohort>.cosmic.summary.json          machine-readable audit totals

COSMIC presence means a mutation was recorded; it is not itself a driver label.
Study independence is not inferred by this script. Exact sample-name overlap is
flagged, while all other matched records remain ``not_assessed`` for independence.
"""

import argparse
import csv
import gzip
import json
import os
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


TRUE_VALUES = {"1", "true", "t", "yes", "y"}

COSMIC_REQUIRED = {
    "GENOMIC_MUTATION_ID",
    "SAMPLE_NAME",
    "COSMIC_SAMPLE_ID",
    "COSMIC_PHENOTYPE_ID",
    "CHROMOSOME",
    "GENOME_START",
    "GENOME_STOP",
    "GENOMIC_WT_ALLELE",
    "GENOMIC_MUT_ALLELE",
    "COSMIC_STUDY_ID",
    "PUBMED_PMID",
    "HGVSG",
    "MUTATION_SOMATIC_STATUS",
}


def open_text(path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") \
        if str(path).endswith(".gz") else open(path, "r", encoding="utf-8", newline="")


def normalize_chrom(value):
    """Return UCSC-style chromosome names, including COSMIC MT -> chrM."""
    chrom = str(value or "").strip()
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    if chrom.upper() in {"M", "MT"}:
        return "chrM"
    return f"chr{chrom}" if chrom else ""


def normalize_allele(value):
    """Map MAF '-' and COSMIC empty alleles onto the same comparison value."""
    allele = str(value or "").strip().upper()
    return "" if allele == "-" else allele


def shown_allele(value):
    return value if value else "-"


def _need_columns(fieldnames, required, label):
    missing = required - set(fieldnames or [])
    if missing:
        raise ValueError(f"{label} lacks required columns: {sorted(missing)}")


def load_candidate_elements(path):
    """Load only rows explicitly marked true in the is_consensus column."""
    elements = {}
    by_chrom = defaultdict(list)
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"element_id", "chrom", "start", "end", "is_consensus"}
        _need_columns(reader.fieldnames, required, "consensus table")
        for row in reader:
            if str(row["is_consensus"]).strip().lower() not in TRUE_VALUES:
                continue
            element_id = row["element_id"].strip()
            chrom = normalize_chrom(row["chrom"])
            start = int(row["start"])
            end = int(row["end"])
            if not element_id or not chrom or end <= start:
                raise ValueError(f"invalid candidate element row: {row}")
            element = {
                "element_id": element_id,
                "chrom": chrom,
                "start": start,
                "end": end,
                "consensus_rank": row.get("consensus_rank", ""),
                "acat_q": row.get("acat_q", row.get("consensus_q", "")),
                "vote_count": row.get("vote_count", ""),
            }
            previous = elements.get(element_id)
            if previous and previous != element:
                raise ValueError(f"candidate element has conflicting rows: {element_id}")
            elements[element_id] = element

    for element in elements.values():
        by_chrom[element["chrom"]].append(element)
    for chrom in by_chrom:
        by_chrom[chrom].sort(key=lambda item: (item["start"], item["end"], item["element_id"]))
    return elements, by_chrom


def overlapping_elements(by_chrom, chrom, pos1):
    """Element ids containing a 1-based position under BED half-open semantics."""
    hits = []
    for element in by_chrom.get(chrom, ()):
        if element["start"] >= pos1:
            break
        if element["start"] < pos1 <= element["end"]:
            hits.append(element["element_id"])
    return hits


def load_candidate_mutations(path, elements_by_chrom):
    variants = {}
    element_events = defaultdict(set)
    seen_events = set()
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample", "chrom", "pos", "ref", "alt", "vtype"}
        _need_columns(reader.fieldnames, required, "mutation table")
        for row in reader:
            chrom = normalize_chrom(row["chrom"])
            try:
                pos = int(row["pos"])
            except (TypeError, ValueError):
                continue
            element_ids = overlapping_elements(elements_by_chrom, chrom, pos)
            if not element_ids:
                continue
            ref = normalize_allele(row["ref"])
            alt = normalize_allele(row["alt"])
            key = (chrom, pos, ref, alt)
            sample = str(row["sample"] or "").strip()
            vtype = str(row["vtype"] or "").strip().upper()
            event_key = (sample, key, vtype)
            if event_key in seen_events:
                continue
            seen_events.add(event_key)

            variant = variants.setdefault(key, {
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "samples": set(),
                "vtypes": set(),
                "element_ids": set(),
            })
            variant["samples"].add(sample)
            variant["vtypes"].add(vtype)
            variant["element_ids"].update(element_ids)
            for element_id in element_ids:
                element_events[element_id].add((sample, key, vtype))
    return variants, element_events, seen_events


def _field(parts, indexes, name):
    index = indexes.get(name)
    return parts[index].strip() if index is not None and index < len(parts) else ""


def scan_cosmic(path, candidate_variants, progress_every=5_000_000):
    """Stream COSMIC once and retain only exact candidate-variant observations."""
    observations = {}
    raw_match_rows = 0
    rows_scanned = 0
    positions = {(key[0], key[1]) for key in candidate_variants}

    if not candidate_variants:
        return observations, rows_scanned, raw_match_rows

    with open_text(path) as handle:
        header_line = handle.readline()
        if not header_line:
            raise ValueError("COSMIC table is empty")
        header = header_line.rstrip("\r\n").split("\t")
        _need_columns(header, COSMIC_REQUIRED, "COSMIC table")
        indexes = {name: i for i, name in enumerate(header)}

        for line in handle:
            rows_scanned += 1
            if progress_every and rows_scanned % progress_every == 0:
                sys.stderr.write(
                    f"cosmic_compare: scanned {rows_scanned:,} COSMIC rows; "
                    f"{raw_match_rows:,} raw exact rows\n")
                sys.stderr.flush()
            parts = line.rstrip("\r\n").split("\t")
            chrom = normalize_chrom(_field(parts, indexes, "CHROMOSOME"))
            try:
                pos = int(_field(parts, indexes, "GENOME_START"))
            except (TypeError, ValueError):
                continue
            if (chrom, pos) not in positions:
                continue
            key = (
                chrom,
                pos,
                normalize_allele(_field(parts, indexes, "GENOMIC_WT_ALLELE")),
                normalize_allele(_field(parts, indexes, "GENOMIC_MUT_ALLELE")),
            )
            if key not in candidate_variants:
                continue
            raw_match_rows += 1

            cosv = _field(parts, indexes, "GENOMIC_MUTATION_ID")
            coss = _field(parts, indexes, "COSMIC_SAMPLE_ID")
            sample_name = _field(parts, indexes, "SAMPLE_NAME")
            phenotype = _field(parts, indexes, "COSMIC_PHENOTYPE_ID")
            study = _field(parts, indexes, "COSMIC_STUDY_ID")
            pmid = _field(parts, indexes, "PUBMED_PMID")
            status = _field(parts, indexes, "MUTATION_SOMATIC_STATUS")
            obs_key = (key, cosv, coss or sample_name, phenotype, study, pmid, status)
            observation = observations.setdefault(obs_key, {
                "variant_key": key,
                "cosv": cosv,
                "cosmic_sample_id": coss,
                "sample_name": sample_name,
                "phenotype_id": phenotype,
                "study_id": study,
                "pmid": pmid,
                "somatic_status": status,
                "genome_stop": _field(parts, indexes, "GENOME_STOP"),
                "gene_symbols": set(),
                "transcripts": set(),
                "mutation_nc_ids": set(),
                "legacy_ids": set(),
                "hgvsg": set(),
                "source_rows": 0,
            })
            observation["source_rows"] += 1
            for column, target in (
                ("GENE_SYMBOL", "gene_symbols"),
                ("TRANSCRIPT_ACCESSION", "transcripts"),
                ("MUTATION_NC_ID", "mutation_nc_ids"),
                ("LEGACY_MUTATION_ID", "legacy_ids"),
                ("HGVSG", "hgvsg"),
            ):
                value = _field(parts, indexes, column)
                if value:
                    observation[target].add(value)
    return observations, rows_scanned, raw_match_rows


def _joined(values):
    return ",".join(sorted(str(value) for value in values if str(value)))


def _sorted_keys(keys):
    def chrom_order(chrom):
        name = chrom[3:] if chrom.startswith("chr") else chrom
        if name.isdigit():
            return (0, int(name))
        return (1, {"X": 23, "Y": 24, "M": 25}.get(name, 1000), name)
    return sorted(keys, key=lambda key: (chrom_order(key[0]), key[1], key[2], key[3]))


def build_output_rows(cohort, elements, variants, element_events, observations):
    by_variant = defaultdict(list)
    for observation in observations.values():
        by_variant[observation["variant_key"]].append(observation)

    mutation_rows = []
    for key in _sorted_keys(variants):
        variant = variants[key]
        obs = by_variant.get(key, [])
        cosmic_samples = {o["cosmic_sample_id"] for o in obs if o["cosmic_sample_id"]}
        confirmed_samples = {
            o["cosmic_sample_id"] for o in obs
            if o["cosmic_sample_id"] and o["somatic_status"] == "Confirmed somatic variant"
        }
        studies = {o["study_id"] for o in obs if o["study_id"]}
        pmids = {o["pmid"] for o in obs if o["pmid"]}
        cosv_ids = {o["cosv"] for o in obs if o["cosv"]}
        phenotype_ids = {o["phenotype_id"] for o in obs if o["phenotype_id"]}
        statuses = {o["somatic_status"] for o in obs if o["somatic_status"]}
        cohort_samples = variant["samples"]
        self_overlap = {
            o["cosmic_sample_id"] or o["sample_name"] for o in obs
            if ((o["cosmic_sample_id"] and o["cosmic_sample_id"] in cohort_samples)
                or (o["sample_name"] and o["sample_name"] in cohort_samples))
        }
        mutation_rows.append({
            "cohort": cohort,
            "chrom": variant["chrom"],
            "pos": variant["pos"],
            "ref": shown_allele(variant["ref"]),
            "alt": shown_allele(variant["alt"]),
            "vtypes": _joined(variant["vtypes"]),
            "element_ids": _joined(variant["element_ids"]),
            "n_elements": len(variant["element_ids"]),
            "n_cohort_samples": len(cohort_samples),
            "cohort_samples": _joined(cohort_samples),
            "cosmic_match": "yes" if obs else "no",
            "n_cosmic_variants": len(cosv_ids),
            "n_cosmic_samples": len(cosmic_samples),
            "n_confirmed_somatic_samples": len(confirmed_samples),
            "n_cosmic_studies": len(studies),
            "n_pubmed_pmids": len(pmids),
            "n_exact_sample_name_overlaps": len(self_overlap),
            "independence_status": "not_assessed" if obs else "no_match",
            "cosmic_variant_ids": _joined(cosv_ids),
            "cosmic_study_ids": _joined(studies),
            "pubmed_pmids": _joined(pmids),
            "cosmic_phenotype_ids": _joined(phenotype_ids),
            "somatic_statuses": _joined(statuses),
        })

    exact_rows = []
    for observation in sorted(
            observations.values(),
            key=lambda o: (_sorted_keys([o["variant_key"]])[0], o["cosv"],
                           o["cosmic_sample_id"], o["study_id"], o["pmid"])):
        key = observation["variant_key"]
        variant = variants[key]
        sample_overlap = (
            observation["cosmic_sample_id"] in variant["samples"]
            or observation["sample_name"] in variant["samples"]
        )
        exact_rows.append({
            "cohort": cohort,
            "chrom": key[0],
            "pos": key[1],
            "genome_stop": observation["genome_stop"],
            "ref": shown_allele(key[2]),
            "alt": shown_allele(key[3]),
            "vtypes": _joined(variant["vtypes"]),
            "element_ids": _joined(variant["element_ids"]),
            "cohort_samples": _joined(variant["samples"]),
            "cosmic_variant_id": observation["cosv"],
            "cosmic_sample_id": observation["cosmic_sample_id"],
            "cosmic_sample_name": observation["sample_name"],
            "cosmic_phenotype_id": observation["phenotype_id"],
            "cosmic_study_id": observation["study_id"],
            "pubmed_pmid": observation["pmid"],
            "somatic_status": observation["somatic_status"],
            "hgvsg": _joined(observation["hgvsg"]),
            "gene_symbols": _joined(observation["gene_symbols"]),
            "transcript_accessions": _joined(observation["transcripts"]),
            "mutation_nc_ids": _joined(observation["mutation_nc_ids"]),
            "legacy_mutation_ids": _joined(observation["legacy_ids"]),
            "source_transcript_rows": observation["source_rows"],
            "exact_sample_name_overlap": "yes" if sample_overlap else "no",
            "independence_status": "possible_self_overlap" if sample_overlap else "not_assessed",
        })

    element_rows = []
    for element in sorted(
            elements.values(),
            key=lambda e: (int(e["consensus_rank"]) if str(e["consensus_rank"]).isdigit()
                           else 10**12, e["element_id"])):
        element_id = element["element_id"]
        keys = {event[1] for event in element_events.get(element_id, set())}
        matched_keys = {key for key in keys if by_variant.get(key)}
        cosmic_samples = {
            observation["cosmic_sample_id"]
            for key in matched_keys for observation in by_variant[key]
            if observation["cosmic_sample_id"]
        }
        element_rows.append({
            "cohort": cohort,
            **element,
            "n_candidate_sample_events": len(element_events.get(element_id, set())),
            "n_distinct_candidate_mutations": len(keys),
            "n_cosmic_matched_mutations": len(matched_keys),
            "n_cosmic_samples": len(cosmic_samples),
            "cosmic_match": "yes" if matched_keys else "no",
        })
    return mutation_rows, exact_rows, element_rows


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


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cosmic", required=True,
                        help="COSMIC NonCodingVariants GRCh38 TSV or TSV.GZ")
    parser.add_argument("--consensus", required=True,
                        help="consensus TSV containing is_consensus and BED coordinates")
    parser.add_argument("--mutations", required=True,
                        help="canonical 1-based mutation TSV or TSV.GZ")
    parser.add_argument("--cohort", required=True, help="cohort label for outputs")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--progress-every", type=int, default=5_000_000,
                        help="report scan progress every N COSMIC rows; 0 disables")
    args = parser.parse_args()

    elements, by_chrom = load_candidate_elements(args.consensus)
    variants, element_events, events = load_candidate_mutations(args.mutations, by_chrom)
    sys.stderr.write(
        f"cosmic_compare: {len(elements)} candidate elements; {len(events)} sample events; "
        f"{len(variants)} distinct variants\n")
    observations, rows_scanned, raw_match_rows = scan_cosmic(
        args.cosmic, variants, args.progress_every)
    mutation_rows, exact_rows, element_rows = build_output_rows(
        args.cohort, elements, variants, element_events, observations)

    prefix = Path(args.outdir) / f"{args.cohort}.cosmic"
    mutation_fields = [
        "cohort", "chrom", "pos", "ref", "alt", "vtypes", "element_ids",
        "n_elements", "n_cohort_samples", "cohort_samples", "cosmic_match",
        "n_cosmic_variants", "n_cosmic_samples", "n_confirmed_somatic_samples",
        "n_cosmic_studies", "n_pubmed_pmids", "n_exact_sample_name_overlaps",
        "independence_status", "cosmic_variant_ids", "cosmic_study_ids",
        "pubmed_pmids", "cosmic_phenotype_ids", "somatic_statuses",
    ]
    exact_fields = [
        "cohort", "chrom", "pos", "genome_stop", "ref", "alt", "vtypes",
        "element_ids", "cohort_samples", "cosmic_variant_id", "cosmic_sample_id",
        "cosmic_sample_name", "cosmic_phenotype_id", "cosmic_study_id",
        "pubmed_pmid", "somatic_status", "hgvsg", "gene_symbols",
        "transcript_accessions", "mutation_nc_ids", "legacy_mutation_ids",
        "source_transcript_rows", "exact_sample_name_overlap", "independence_status",
    ]
    element_fields = [
        "cohort", "consensus_rank", "element_id", "chrom", "start", "end",
        "acat_q", "vote_count", "n_candidate_sample_events",
        "n_distinct_candidate_mutations", "n_cosmic_matched_mutations",
        "n_cosmic_samples", "cosmic_match",
    ]
    write_tsv(f"{prefix}.mutation_summary.tsv", mutation_rows, mutation_fields)
    write_tsv(f"{prefix}.exact_matches.tsv", exact_rows, exact_fields)
    write_tsv(f"{prefix}.element_summary.tsv", element_rows, element_fields)

    matched_mutations = sum(row["cosmic_match"] == "yes" for row in mutation_rows)
    confirmed_mutations = sum(row["n_confirmed_somatic_samples"] > 0 for row in mutation_rows)
    vtype_events = Counter()
    for variant in variants.values():
        for vtype in variant["vtypes"]:
            vtype_events[vtype] += len(variant["samples"])
    summary = {
        "cohort": args.cohort,
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "cosmic_file": str(Path(args.cosmic).resolve()),
        "consensus_file": str(Path(args.consensus).resolve()),
        "mutations_file": str(Path(args.mutations).resolve()),
        "coordinate_rule": "BED candidate start < 1-based mutation position <= BED end",
        "exact_match_rule": "GRCh38 chromosome + GENOME_START + reference + alternate allele",
        "candidate_elements": len(elements),
        "candidate_sample_events": len(events),
        "candidate_distinct_mutations": len(variants),
        "candidate_sample_events_by_vtype": dict(sorted(vtype_events.items())),
        "cosmic_rows_scanned": rows_scanned,
        "cosmic_raw_exact_rows": raw_match_rows,
        "cosmic_deduplicated_observations": len(observations),
        "cosmic_matched_distinct_mutations": matched_mutations,
        "cosmic_unmatched_distinct_mutations": len(variants) - matched_mutations,
        "cosmic_confirmed_somatic_distinct_mutations": confirmed_mutations,
        "interpretation": (
            "COSMIC presence records prior observation; it is not a driver classification. "
            "Study independence remains not_assessed unless separately resolved."
        ),
    }
    summary_path = Path(f"{prefix}.summary.json")
    temporary = summary_path.with_name(summary_path.name + ".tmp")
    temporary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, summary_path)
    sys.stderr.write(json.dumps(summary, sort_keys=True) + "\n")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, csv.Error) as error:
        sys.exit(f"ERROR: {error}")
