#!/usr/bin/env python3
"""COSMIC known-cancer context of the candidate mutations (COSMICCONTEXT).

Exact GRCh38 matching of the patient-resolved candidate mutations (the EXPLAIN candidate
table: one row per sample x mutation x significant element fragment) against the user's own
COSMIC Non-Coding Variants file. The layer annotates; it never filters or re-ranks an element,
and COSMIC presence is a prior observation, not a driver label.

Ported from benchmarks/cosmic_compare.py with these changes:
  * input is the candidate table (fragment-resolved), not consensus spans + the mutation table;
  * the mutation id column may be MUTATION_NC_ID (the v104 file) or COSMIC_NCV_ID (the v104
    download-page documentation); every column is read by name;
  * lines may end in CRLF (the v104 file does);
  * the file name must be COSMIC's own, Cosmic_NonCodingVariants_v<NNN>_GRCh38.tsv[.gz]: the
    release must equal --release and a GRCh37 file stops the run;
  * every CHROMOSOME value must be a COSMIC GRCh38 name (1-22, X, Y, MT);
  * the MD5 of the file is computed while it is read and checked against --expected-md5;
  * a provenance JSON records release, file name, checksums, row counts and the rules.

Matching: same chromosome (COSMIC 5 -> chr5, MT -> chrM), GENOME_START, reference and
alternate allele, with the candidate table's '-' equal to COSMIC's empty allele. No window and
no interval-only match is ever called a match. COSMIC NCV rows are transcript-level, so
observations are collapsed by (variant, COSV, COSS or sample name, phenotype, study, PMID,
somatic status).

Same-cohort flag: an observation is possible_self_overlap when its COSMIC SAMPLE_NAME or
COSMIC_SAMPLE_ID equals a cohort sample id of that mutation exactly; every other matched
observation is not_assessed. No fuzzy, donor-prefix or crosswalk matching (the PCAWG crosswalk
stays in benchmarks/cosmic_sample_audit.py).

Outputs (prefix <outdir>/<cohort>.cosmic):
  .mutation_summary.tsv   counts per distinct candidate variant (no COSMIC identifiers)
  .element_summary.tsv    counts per significant element (no COSMIC identifiers)
  .provenance.json        release, file name, MD5s, row counts, rules
  .exact_matches.tsv      one row per collapsed COSMIC observation: COSMIC RECORDS (restricted)
  NOTICE.txt              licence notice for the restricted table
"""

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

VERSION = "1.0.0"

FILE_NAME = re.compile(r"^Cosmic_NonCodingVariants_(v\d+)_(GRCh3[78])\.tsv(\.gz)?$")
RELEASE = re.compile(r"^v\d+$")
MD5 = re.compile(r"^[0-9a-f]{32}$")
COSMIC_CHROMOSOMES = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}
ID_COLUMNS = ("MUTATION_NC_ID", "COSMIC_NCV_ID")
CONFIRMED = "Confirmed somatic variant"

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
CANDIDATE_REQUIRED = {"sample", "chrom", "pos", "ref", "alt", "vtype", "element_id"}

NOTICE = """This folder contains data from COSMIC {release} (Catalogue Of Somatic Mutations In Cancer,
https://cancer.sanger.ac.uk/cosmic), matched against the candidate mutations of this run with
the user's own licensed copy of the COSMIC Non-Coding Variants file.

Do not redistribute these files. The COSMIC licence (https://www.cosmickb.org/terms/) allows
academic users to disclose COSMIC data only to the extent necessary to enable publication of
their scientific papers ("Publication Enablement"), with an acknowledgement of COSMIC that
names the release ({release}); it forbids any other transfer, sharing or display to third
parties or on a free-to-access platform, including public websites and public repositories.
Commercial use needs a commercial licence. COSMIC identifiers (COSV, COSS, ...) are part of
COSMIC.

The count tables one folder up (cosmic/) hold no COSMIC identifiers; cite COSMIC and the
release when you report them.
"""


class HashingReader(io.RawIOBase):
    """Read-through wrapper that computes the MD5 of the bytes read from a binary file."""

    def __init__(self, handle):
        self.handle = handle
        self.md5 = hashlib.md5()

    def readable(self):
        return True

    def readinto(self, buffer):
        data = self.handle.read(len(buffer))
        n = len(data)
        buffer[:n] = data
        self.md5.update(data)
        return n

    def drain(self):
        """Hash any bytes the decompressor did not consume (trailing data)."""
        for block in iter(lambda: self.handle.read(1 << 20), b""):
            self.md5.update(block)
        return self.md5.hexdigest()


def check_file_name(path, release):
    """Release and assembly from COSMIC's file name; stop on GRCh37 or a release mismatch."""
    name = Path(path).name
    match = FILE_NAME.match(name)
    if not match:
        raise ValueError(
            f"COSMIC file name '{name}' is not COSMIC's own "
            "(Cosmic_NonCodingVariants_v<release>_GRCh38.tsv.gz). Keep the name of the file "
            "extracted from the download tar, so that the release and the assembly can be checked.")
    file_release, assembly = match.group(1), match.group(2)
    if assembly != "GRCh38":
        raise ValueError(
            f"COSMIC file '{name}' is the {assembly} file set. The pipeline analyses GRCh38 and "
            "never lifts COSMIC: download the GRCh38 Non-Coding Variants file.")
    if not RELEASE.match(release or ""):
        raise ValueError(f"--cosmic_release '{release}' is not of the form v<number> (e.g. v104)")
    if file_release != release:
        raise ValueError(
            f"COSMIC file '{name}' is release {file_release}, but --cosmic_release is {release}.")
    return file_release, assembly


def normalize_chrom(value):
    """Candidate-table chromosome (UCSC names) -> UCSC name; COSMIC names are mapped separately."""
    chrom = str(value or "").strip()
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    if chrom.upper() in {"M", "MT"}:
        return "chrM"
    return f"chr{chrom}" if chrom else ""


def cosmic_chrom(value, line_number):
    """COSMIC GRCh38 CHROMOSOME (1-22, X, Y, MT) -> UCSC name; anything else stops the run."""
    chrom = value.strip()
    if chrom not in COSMIC_CHROMOSOMES:
        raise ValueError(
            f"COSMIC line {line_number}: CHROMOSOME '{chrom}' is not a COSMIC GRCh38 chromosome "
            "name (1-22, X, Y, MT). Is this the COSMIC Non-Coding Variants file, unmodified?")
    return "chrM" if chrom == "MT" else "chr" + chrom


def normalize_allele(value):
    """Map the candidate table's '-' and COSMIC's empty allele onto the same value."""
    allele = str(value or "").strip().upper()
    return "" if allele == "-" else allele


def shown_allele(value):
    return value if value else "-"


def _need_columns(fieldnames, required, label):
    missing = required - set(fieldnames or [])
    if missing:
        raise ValueError(f"{label} lacks required columns: {sorted(missing)}")


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return open(path, "r", encoding="utf-8", newline="")


def load_candidates(path):
    """Distinct candidate variants and the sample events per element from the candidate table."""
    variants = {}
    element_events = defaultdict(set)
    events = set()
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        _need_columns(reader.fieldnames, CANDIDATE_REQUIRED, "candidate table")
        for row in reader:
            chrom = normalize_chrom(row["chrom"])
            pos = int(row["pos"])
            ref = normalize_allele(row["ref"])
            alt = normalize_allele(row["alt"])
            key = (chrom, pos, ref, alt)
            sample = str(row["sample"] or "").strip()
            vtype = str(row["vtype"] or "").strip().upper()
            element_id = str(row["element_id"] or "").strip()
            variant = variants.setdefault(key, {"samples": set(), "vtypes": set(), "element_ids": set()})
            variant["samples"].add(sample)
            variant["vtypes"].add(vtype)
            variant["element_ids"].add(element_id)
            events.add((sample, key, vtype))
            element_events[element_id].add((sample, key, vtype))
    return variants, element_events, events


def load_consensus(path):
    """Significant elements (is_consensus True) with rank, span and consensus_q, in file order."""
    elements = {}
    if not path:
        return elements
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        _need_columns(reader.fieldnames, {"element_id", "is_consensus"}, "consensus table")
        for row in reader:
            if str(row["is_consensus"]).strip() != "True":
                continue
            elements[row["element_id"]] = {key: row.get(key) or "NA"
                                           for key in ("consensus_rank", "chrom", "start", "end", "consensus_q")}
    return elements


def _field(parts, indexes, name):
    index = indexes.get(name)
    return parts[index].strip() if index is not None and index < len(parts) else ""


def scan_cosmic(path, candidate_variants, progress_every=5_000_000):
    """Stream the COSMIC file once: header and chromosome checks, MD5, exact observations."""
    observations = {}
    stats = {"rows_scanned": 0, "raw_exact_rows": 0, "crlf_lines": 0}
    positions = {(key[0], key[1]) for key in candidate_variants}
    with open(path, "rb") as raw_handle:
        hashing = HashingReader(raw_handle)
        buffered = io.BufferedReader(hashing, buffer_size=1 << 20)
        binary = gzip.GzipFile(fileobj=buffered, mode="rb") if str(path).endswith(".gz") else buffered
        handle = io.TextIOWrapper(binary, encoding="utf-8", errors="replace", newline="")
        header_line = handle.readline()
        if not header_line:
            raise ValueError("COSMIC table is empty")
        header = header_line.rstrip("\r\n").split("\t")
        _need_columns(header, COSMIC_REQUIRED, "COSMIC table")
        id_column = next((name for name in ID_COLUMNS if name in header), None)
        if id_column is None:
            raise ValueError(f"COSMIC table lacks a mutation id column ({' or '.join(ID_COLUMNS)})")
        stats["id_column"] = id_column
        indexes = {name: i for i, name in enumerate(header)}
        ci, si = indexes["CHROMOSOME"], indexes["GENOME_START"]
        known = {}
        line_number = 1
        for line in handle:
            line_number += 1
            stats["rows_scanned"] += 1
            if line.endswith("\r\n"):
                stats["crlf_lines"] += 1
            if progress_every and stats["rows_scanned"] % progress_every == 0:
                sys.stderr.write(f"cosmic_context: scanned {stats['rows_scanned']:,} COSMIC rows; "
                                 f"{stats['raw_exact_rows']:,} exact rows\n")
                sys.stderr.flush()
            parts = line.rstrip("\r\n").split("\t")
            value = parts[ci] if ci < len(parts) else ""
            chrom = known.get(value)
            if chrom is None:
                chrom = cosmic_chrom(value, line_number)
                known[value] = chrom
            try:
                pos = int(parts[si])
            except (IndexError, ValueError):
                raise ValueError(f"COSMIC line {line_number}: GENOME_START '{parts[si] if si < len(parts) else ''}' is not an integer")
            if (chrom, pos) not in positions:
                continue
            key = (chrom, pos,
                   normalize_allele(_field(parts, indexes, "GENOMIC_WT_ALLELE")),
                   normalize_allele(_field(parts, indexes, "GENOMIC_MUT_ALLELE")))
            if key not in candidate_variants:
                continue
            stats["raw_exact_rows"] += 1
            cosv = _field(parts, indexes, "GENOMIC_MUTATION_ID")
            coss = _field(parts, indexes, "COSMIC_SAMPLE_ID")
            sample_name = _field(parts, indexes, "SAMPLE_NAME")
            phenotype = _field(parts, indexes, "COSMIC_PHENOTYPE_ID")
            study = _field(parts, indexes, "COSMIC_STUDY_ID")
            pmid = _field(parts, indexes, "PUBMED_PMID")
            status = _field(parts, indexes, "MUTATION_SOMATIC_STATUS")
            obs_key = (key, cosv, coss or sample_name, phenotype, study, pmid, status)
            observation = observations.setdefault(obs_key, {
                "variant_key": key, "cosv": cosv, "cosmic_sample_id": coss, "sample_name": sample_name,
                "phenotype_id": phenotype, "study_id": study, "pmid": pmid, "somatic_status": status,
                "genome_stop": _field(parts, indexes, "GENOME_STOP"),
                "gene_symbols": set(), "transcripts": set(), "mutation_nc_ids": set(),
                "legacy_ids": set(), "hgvsg": set(), "source_rows": 0,
            })
            observation["source_rows"] += 1
            for column, target in (("GENE_SYMBOL", "gene_symbols"), ("TRANSCRIPT_ACCESSION", "transcripts"),
                                   (id_column, "mutation_nc_ids"), ("LEGACY_MUTATION_ID", "legacy_ids"),
                                   ("HGVSG", "hgvsg")):
                value = _field(parts, indexes, column)
                if value:
                    observation[target].add(value)
        # read to the end of the file so that the MD5 covers every byte
        for _ in iter(lambda: binary.read(1 << 20), b""):
            pass
        stats["md5"] = hashing.drain()
    stats["line_endings"] = ("CRLF" if stats["crlf_lines"] == stats["rows_scanned"] else
                             "LF" if stats["crlf_lines"] == 0 else "mixed")
    return observations, stats


def _joined(values):
    return ",".join(sorted(str(value) for value in values if str(value)))


def _chrom_order(chrom):
    name = chrom[3:] if chrom.startswith("chr") else chrom
    if name.isdigit():
        return (0, int(name), "")
    return (1, {"X": 23, "Y": 24, "M": 25}.get(name, 1000), name)


def _sorted_keys(keys):
    return sorted(keys, key=lambda key: (_chrom_order(key[0]), key[1], key[2], key[3]))


def _self_overlap(observation, samples):
    return bool((observation["cosmic_sample_id"] and observation["cosmic_sample_id"] in samples)
                or (observation["sample_name"] and observation["sample_name"] in samples))


def build_rows(cohort, release, consensus, variants, element_events, observations):
    by_variant = defaultdict(list)
    for observation in observations.values():
        by_variant[observation["variant_key"]].append(observation)

    mutation_rows = []
    for key in _sorted_keys(variants):
        variant = variants[key]
        obs = by_variant.get(key, [])
        samples = {o["cosmic_sample_id"] or o["sample_name"] for o in obs if o["cosmic_sample_id"] or o["sample_name"]}
        confirmed = {o["cosmic_sample_id"] or o["sample_name"] for o in obs
                     if (o["cosmic_sample_id"] or o["sample_name"]) and o["somatic_status"] == CONFIRMED}
        overlap = {o["cosmic_sample_id"] or o["sample_name"] for o in obs if _self_overlap(o, variant["samples"])}
        status = "no_match" if not obs else ("possible_self_overlap" if overlap else "not_assessed")
        mutation_rows.append({
            "cohort": cohort, "cosmic_release": release,
            "chrom": key[0], "pos": key[1], "ref": shown_allele(key[2]), "alt": shown_allele(key[3]),
            "vtypes": _joined(variant["vtypes"]), "element_ids": _joined(variant["element_ids"]),
            "n_elements": len(variant["element_ids"]), "n_cohort_samples": len(variant["samples"]),
            "cohort_samples": _joined(variant["samples"]),
            "cosmic_match": "yes" if obs else "no",
            "n_cosmic_mutation_ids": len({o["cosv"] for o in obs if o["cosv"]}),
            "n_cosmic_samples": len(samples),
            "n_confirmed_somatic_samples": len(confirmed),
            "n_cosmic_studies": len({o["study_id"] for o in obs if o["study_id"]}),
            "n_pubmed_pmids": len({o["pmid"] for o in obs if o["pmid"]}),
            "n_exact_sample_name_overlaps": len(overlap),
            "independence_status": status,
        })

    exact_rows = []
    for o in sorted(observations.values(),
                    key=lambda o: (_sorted_keys([o["variant_key"]])[0], o["cosv"], o["cosmic_sample_id"],
                                   o["sample_name"], o["study_id"], o["pmid"], o["somatic_status"])):
        key = o["variant_key"]
        variant = variants[key]
        overlap = _self_overlap(o, variant["samples"])
        exact_rows.append({
            "cohort": cohort, "cosmic_release": release,
            "chrom": key[0], "pos": key[1], "genome_stop": o["genome_stop"],
            "ref": shown_allele(key[2]), "alt": shown_allele(key[3]),
            "vtypes": _joined(variant["vtypes"]), "element_ids": _joined(variant["element_ids"]),
            "cohort_samples": _joined(variant["samples"]),
            "cosmic_mutation_id": o["cosv"], "cosmic_sample_id": o["cosmic_sample_id"],
            "cosmic_sample_name": o["sample_name"], "cosmic_phenotype_id": o["phenotype_id"],
            "cosmic_study_id": o["study_id"], "pubmed_pmid": o["pmid"], "somatic_status": o["somatic_status"],
            "hgvsg": _joined(o["hgvsg"]), "gene_symbols": _joined(o["gene_symbols"]),
            "transcript_accessions": _joined(o["transcripts"]), "mutation_nc_ids": _joined(o["mutation_nc_ids"]),
            "legacy_mutation_ids": _joined(o["legacy_ids"]), "source_transcript_rows": o["source_rows"],
            "exact_sample_name_overlap": "yes" if overlap else "no",
            "independence_status": "possible_self_overlap" if overlap else "not_assessed",
        })

    # significant elements: consensus order when given, otherwise the candidate table's ids
    element_ids = list(consensus) if consensus else sorted(element_events)
    for element_id in sorted(element_events):
        if element_id not in element_ids:
            element_ids.append(element_id)
    element_rows = []
    for element_id in element_ids:
        info = consensus.get(element_id, {})
        events = element_events.get(element_id, set())
        keys = {event[1] for event in events}
        matched = {key for key in keys if by_variant.get(key)}
        obs = [o for key in matched for o in by_variant[key]]
        element_rows.append({
            "cohort": cohort, "cosmic_release": release, "element_id": element_id,
            "consensus_rank": info.get("consensus_rank", "NA"), "chrom": info.get("chrom", "NA"),
            "start": info.get("start", "NA"), "end": info.get("end", "NA"),
            "consensus_q": info.get("consensus_q", "NA"),
            "n_candidate_sample_events": len(events),
            "n_distinct_candidate_mutations": len(keys),
            "n_cosmic_matched_mutations": len(matched),
            "n_confirmed_somatic_mutations": len({o["variant_key"] for o in obs if o["somatic_status"] == CONFIRMED}),
            "n_cosmic_samples": len({o["cosmic_sample_id"] or o["sample_name"] for o in obs
                                     if o["cosmic_sample_id"] or o["sample_name"]}),
            "n_exact_sample_name_overlaps": len({o["cosmic_sample_id"] or o["sample_name"] for o in obs
                                                 if _self_overlap(o, variants[o["variant_key"]]["samples"])}),
            "cosmic_match": "yes" if matched else "no",
        })
    return mutation_rows, exact_rows, element_rows


MUTATION_FIELDS = [
    "cohort", "cosmic_release", "chrom", "pos", "ref", "alt", "vtypes", "element_ids", "n_elements",
    "n_cohort_samples", "cohort_samples", "cosmic_match", "n_cosmic_mutation_ids", "n_cosmic_samples",
    "n_confirmed_somatic_samples", "n_cosmic_studies", "n_pubmed_pmids", "n_exact_sample_name_overlaps",
    "independence_status",
]
EXACT_FIELDS = [
    "cohort", "cosmic_release", "chrom", "pos", "genome_stop", "ref", "alt", "vtypes", "element_ids",
    "cohort_samples", "cosmic_mutation_id", "cosmic_sample_id", "cosmic_sample_name", "cosmic_phenotype_id",
    "cosmic_study_id", "pubmed_pmid", "somatic_status", "hgvsg", "gene_symbols", "transcript_accessions",
    "mutation_nc_ids", "legacy_mutation_ids", "source_transcript_rows", "exact_sample_name_overlap",
    "independence_status",
]
ELEMENT_FIELDS = [
    "cohort", "cosmic_release", "element_id", "consensus_rank", "chrom", "start", "end", "consensus_q",
    "n_candidate_sample_events", "n_distinct_candidate_mutations", "n_cosmic_matched_mutations",
    "n_confirmed_somatic_mutations", "n_cosmic_samples", "n_exact_sample_name_overlaps", "cosmic_match",
]


def write_tsv(path, rows, fieldnames):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with open(temporary, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cosmic", required=True, help="Cosmic_NonCodingVariants_v<NNN>_GRCh38.tsv.gz (user-supplied)")
    parser.add_argument("--release", required=True, help="COSMIC release of the file, e.g. v104")
    parser.add_argument("--expected-md5", default="", help="MD5 of the .tsv.gz; a mismatch stops the run")
    parser.add_argument("--candidates", required=True, help="<cohort>.candidates.tsv of EXPLAIN")
    parser.add_argument("--consensus", default="", help="<cohort>.consensus.tsv (rank, span, consensus_q)")
    parser.add_argument("--cohort", required=True)
    parser.add_argument("--outdir", default=".")
    parser.add_argument("--progress-every", type=int, default=5_000_000)
    args = parser.parse_args(argv)

    release, assembly = check_file_name(args.cosmic, args.release)
    expected = (args.expected_md5 or "").strip().lower()
    if expected and not MD5.match(expected):
        raise ValueError(f"--cosmic_ncv_md5 '{args.expected_md5}' is not a 32-character hexadecimal MD5")
    consensus = load_consensus(args.consensus)
    variants, element_events, events = load_candidates(args.candidates)
    sys.stderr.write(f"cosmic_context: {len(element_events)} elements; {len(events)} sample events; "
                     f"{len(variants)} distinct variants\n")
    observations, stats = scan_cosmic(args.cosmic, variants, args.progress_every)
    if expected and stats["md5"] != expected:
        raise ValueError(f"MD5 of {Path(args.cosmic).name} is {stats['md5']}, expected {expected} "
                         "(--cosmic_ncv_md5): the file is incomplete or not the one you meant")
    mutation_rows, exact_rows, element_rows = build_rows(
        args.cohort, release, consensus, variants, element_events, observations)

    outdir = Path(args.outdir)
    prefix = outdir / f"{args.cohort}.cosmic"
    write_tsv(f"{prefix}.mutation_summary.tsv", mutation_rows, MUTATION_FIELDS)
    write_tsv(f"{prefix}.element_summary.tsv", element_rows, ELEMENT_FIELDS)
    write_tsv(f"{prefix}.exact_matches.tsv", exact_rows, EXACT_FIELDS)
    (outdir / "NOTICE.txt").write_text(NOTICE.format(release=release))

    vtype_events = Counter(event[2] for event in events)
    provenance = {
        "layer": "cosmic_context",
        "script": f"cosmic_context.py {VERSION}",
        "cohort": args.cohort,
        "cosmic_release": release,
        "cosmic_assembly": assembly,
        "cosmic_file": Path(args.cosmic).name,
        "cosmic_file_md5": stats["md5"],
        "cosmic_file_md5_expected": expected or None,
        "md5_check": "match" if expected else "not_requested",
        "cosmic_mutation_id_column": stats["id_column"],
        "cosmic_line_endings": stats["line_endings"],
        "cosmic_rows_scanned": stats["rows_scanned"],
        "cosmic_raw_exact_rows": stats["raw_exact_rows"],
        "cosmic_collapsed_observations": len(observations),
        "candidate_elements": len(element_rows),
        "candidate_sample_events": len(events),
        "candidate_sample_events_by_vtype": dict(sorted(vtype_events.items())),
        "candidate_distinct_mutations": len(variants),
        "cosmic_matched_distinct_mutations": sum(row["cosmic_match"] == "yes" for row in mutation_rows),
        "cosmic_unmatched_distinct_mutations": sum(row["cosmic_match"] == "no" for row in mutation_rows),
        "cosmic_confirmed_somatic_distinct_mutations": sum(row["n_confirmed_somatic_samples"] > 0 for row in mutation_rows),
        "possible_self_overlap_distinct_mutations": sum(row["independence_status"] == "possible_self_overlap"
                                                        for row in mutation_rows),
        "rules": {
            "input": "EXPLAIN candidate table (sample x mutation x significant element fragment)",
            "exact_match": "GRCh38 chromosome + GENOME_START + reference + alternate allele; '-' == empty COSMIC allele",
            "collapse": "COSMIC rows collapsed by variant, COSV, COSS (or sample name), phenotype, study, PMID, somatic status",
            "self_overlap": "COSMIC SAMPLE_NAME or COSMIC_SAMPLE_ID equal to a cohort sample id (exact); otherwise not_assessed",
            "use": "annotation only: never a filter, a prior or a weight",
        },
        "outputs": {
            "counts (no COSMIC identifiers)": [f"{args.cohort}.cosmic.mutation_summary.tsv",
                                               f"{args.cohort}.cosmic.element_summary.tsv"],
            "restricted (COSMIC records, do not redistribute)": [f"{args.cohort}.cosmic.exact_matches.tsv"],
        },
        "acknowledgement": f"Contains data from COSMIC {release} (Sondka et al. 2024, doi:10.1093/nar/gkad986) "
                           "used under the COSMIC licence; cite COSMIC and the release when reporting.",
    }
    path = Path(f"{prefix}.provenance.json")
    path.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    sys.stderr.write(json.dumps({k: v for k, v in provenance.items() if isinstance(v, int)}, sort_keys=True) + "\n")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, csv.Error, EOFError, gzip.BadGzipFile) as error:
        sys.exit(f"ERROR: {error}")
