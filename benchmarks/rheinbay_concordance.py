#!/usr/bin/env python3
"""Compare Rheinbay et al. 2020 Table 5 loci with noncodingdrivers consensus results.

The publication used hg19/PCAWG functional elements whereas noncodingdrivers uses
GRCh38 Ensembl Regulatory Build elements.  This utility therefore compares physical
loci after liftOver; it does not compare element identifiers or exact q-values.
"""

import argparse
import csv
import gzip
import json
import re
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path


XML_NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
TABLE5_SHEET = "xl/worksheets/sheet5.xml"


def read_table5(workbook, tissue):
    """Read one tissue from Supplementary Table 5 without an Excel dependency."""
    with zipfile.ZipFile(workbook) as archive:
        strings = []
        root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
        for item in root.findall(XML_NS + "si"):
            strings.append("".join(node.text or "" for node in item.iter(XML_NS + "t")))

        sheet = ET.fromstring(archive.read(TABLE5_SHEET))
        matrix = []
        for row in sheet.iter(XML_NS + "row"):
            values = {}
            for cell in row.findall(XML_NS + "c"):
                column = re.match(r"[A-Z]+", cell.attrib["r"]).group(0)
                value_node = cell.find(XML_NS + "v")
                value = "" if value_node is None else (value_node.text or "")
                if cell.attrib.get("t") == "s" and value:
                    value = strings[int(value)]
                values[column] = value
            matrix.append(values)

    if not matrix:
        return []
    headers = matrix[0]
    rows = []
    for values in matrix[1:]:
        record = {header: values.get(column, "") for column, header in headers.items()}
        if record.get("tissue") == tissue:
            rows.append(record)
    return rows


def bed12_blocks(fields):
    """Return the true BED12 blocks as zero-based, half-open intervals."""
    chrom, start, end = fields[0], int(fields[1]), int(fields[2])
    if len(fields) < 12:
        return [(chrom, start, end)]
    count = int(fields[9])
    sizes = [int(value) for value in fields[10].rstrip(",").split(",") if value]
    offsets = [int(value) for value in fields[11].rstrip(",").split(",") if value]
    if count != len(sizes) or count != len(offsets):
        raise ValueError(f"invalid BED12 block columns for {fields[3]}")
    return [(chrom, start + offset, start + offset + size)
            for size, offset in zip(sizes, offsets)]


def read_elements(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    elements = {}
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            elements[fields[3]] = {
                "chrom": fields[0], "start": int(fields[1]), "end": int(fields[2]),
                "blocks": bed12_blocks(fields),
            }
    return elements


def unique_mapping(lifter, chrom, pos):
    mappings = lifter.convert_coordinate(chrom, pos) or []
    primary = [item for item in mappings if re.fullmatch(r"chr(?:[1-9]|1[0-9]|2[0-2]|X|Y)", item[0])]
    candidates = primary or mappings
    if len(candidates) != 1:
        return None
    return candidates[0]


def lift_blocks(blocks, lifter):
    """Lift BED blocks by mapping both inclusive endpoints."""
    lifted = []
    for chrom, start, end in blocks:
        left = unique_mapping(lifter, chrom, start)
        right = unique_mapping(lifter, chrom, end - 1)
        if left is None or right is None:
            return [], "unmapped_or_ambiguous"
        if left[0] != right[0] or left[2] != right[2]:
            return [], "split_mapping"
        target_start = min(int(left[1]), int(right[1]))
        target_end = max(int(left[1]), int(right[1])) + 1
        lifted.append((left[0], target_start, target_end))
    return lifted, "lifted"


def intervals_text(intervals):
    return ";".join(f"{chrom}:{start}-{end}" for chrom, start, end in intervals)


def overlaps(left, right):
    return left[0] == right[0] and left[1] < right[2] and right[1] < left[2]


def truth(value):
    return str(value).strip().lower() in {"1", "true", "t", "yes"}


def number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def read_tsv(path):
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def published_status(row):
    post_q = number(row.get("Post-filter q-value"))
    final_pass = truth(row.get("Final filter judgement"))
    if not final_pass:
        return "filtered"
    if post_q is not None and post_q < 0.1:
        return "retained_significant"
    return "retained_not_significant"


def joined(rows, field):
    values = sorted({row.get(field, "") for row in rows if row.get(field, "")})
    return ";".join(values)


OUTPUT_FIELDS = [
    "published_id", "target", "tissue", "gene", "pre_filter_q", "post_filter_q",
    "published_mutations", "published_patients", "final_filter_judgement",
    "published_status", "hg19_span", "hg19_blocks", "coordinate_status",
    "hg38_blocks", "pipeline_recovered", "pipeline_match_count",
    "pipeline_element_ids", "pipeline_intervals", "pipeline_ranks", "pipeline_consensus_q", "pipeline_acat_q",
    "pipeline_n_snv", "pipeline_candidate_rows", "pipeline_distinct_patients",
    "pipeline_distinct_variants", "spectrum_statuses", "motif_statuses",
    "nearest_genes", "nearest_gene_statuses",
]


def compare(table_rows, elements, lifter, consensus_rows, candidate_rows):
    significant = [row for row in consensus_rows if truth(row.get("is_consensus"))]
    candidate_by_element = {}
    for row in candidate_rows:
        candidate_by_element.setdefault(row.get("element_id", ""), []).append(row)

    output = []
    for published in table_rows:
        source = elements.get(published["ID"])
        lifted = []
        if source is None:
            coord_status = "not_in_retained_pcawg_catalogue"
            hg19_span = hg19_blocks = ""
        else:
            hg19_span = f"{source['chrom']}:{source['start']}-{source['end']}"
            hg19_blocks = intervals_text(source["blocks"])
            lifted, coord_status = lift_blocks(source["blocks"], lifter)

        matches = []
        for result in significant:
            interval = (result["chrom"], int(float(result["start"])), int(float(result["end"])))
            if any(overlaps(block, interval) for block in lifted):
                matches.append(result)
        matches.sort(key=lambda row: int(float(row["consensus_rank"])))
        match_ids = [row["element_id"] for row in matches]
        candidate_matches = [row for element_id in match_ids
                             for row in candidate_by_element.get(element_id, [])]
        patient_keys = {row.get("sample", "") for row in candidate_matches if row.get("sample")}
        variant_keys = {(row.get("chrom"), row.get("pos"), row.get("ref"), row.get("alt"))
                        for row in candidate_matches}

        output.append({
            "published_id": published.get("ID", ""),
            "target": published.get("target", ""),
            "tissue": published.get("tissue", ""),
            "gene": published.get("gene", ""),
            "pre_filter_q": published.get("Pre-filter q-value", ""),
            "post_filter_q": published.get("Post-filter q-value", ""),
            "published_mutations": published.get("Mutations", ""),
            "published_patients": published.get("Patients mutated", ""),
            "final_filter_judgement": published.get("Final filter judgement", ""),
            "published_status": published_status(published),
            "hg19_span": hg19_span,
            "hg19_blocks": hg19_blocks,
            "coordinate_status": coord_status,
            "hg38_blocks": intervals_text(lifted),
            "pipeline_recovered": str(bool(matches)).lower(),
            "pipeline_match_count": str(len(matches)),
            "pipeline_element_ids": ";".join(match_ids),
            "pipeline_intervals": ";".join(
                f"{row['chrom']}:{row['start']}-{row['end']}" for row in matches),
            "pipeline_ranks": ";".join(row["consensus_rank"] for row in matches),
            "pipeline_consensus_q": ";".join(row.get("consensus_q", "") for row in matches),
            "pipeline_acat_q": ";".join(row.get("acat_q", "") for row in matches),
            "pipeline_n_snv": ";".join(row.get("n_snv", "") for row in matches),
            "pipeline_candidate_rows": str(len(candidate_matches)),
            "pipeline_distinct_patients": str(len(patient_keys)),
            "pipeline_distinct_variants": str(len(variant_keys)),
            "spectrum_statuses": joined(candidate_matches, "spectrum_status"),
            "motif_statuses": joined(candidate_matches, "motif_status"),
            "nearest_genes": joined(candidate_matches, "nearest_gene"),
            "nearest_gene_statuses": joined(candidate_matches, "nearest_gene_status"),
        })
    return output


def summarize(rows):
    resolved = [row for row in rows if row["coordinate_status"] == "lifted"]
    retained = [row for row in rows if row["published_status"] == "retained_significant"]
    resolved_retained = [row for row in retained if row["coordinate_status"] == "lifted"]
    recovered = [row for row in resolved if truth(row["pipeline_recovered"])]
    recovered_retained = [row for row in resolved_retained if truth(row["pipeline_recovered"])]
    annotations = {element_id for row in recovered
                   for element_id in row["pipeline_element_ids"].split(";") if element_id}
    intervals = set()
    for row in recovered:
        for value in row["pipeline_intervals"].split(";"):
            match = re.fullmatch(r"([^:]+):(\d+)-(\d+)", value)
            if match:
                intervals.add((match.group(1), int(match.group(2)), int(match.group(3))))
    locus_total = 0
    current_chrom, current_end = None, -1
    for chrom, start, end in sorted(intervals):
        if chrom != current_chrom or start >= current_end:
            locus_total += 1
            current_chrom, current_end = chrom, end
        else:
            current_end = max(current_end, end)
    return {
        "table5_liver_rows": len(rows),
        "coordinate_resolved": len(resolved),
        "coordinate_unresolved": len(rows) - len(resolved),
        "resolved_table5_rows_recovered": len(recovered),
        "distinct_pipeline_annotations_matched": len(annotations),
        "distinct_pipeline_loci_matched": locus_total,
        "published_retained_significant": len(retained),
        "retained_significant_coordinate_resolved": len(resolved_retained),
        "retained_significant_recovered": len(recovered_retained),
        "recovered_genes": [row["gene"] for row in recovered],
        "recovered_retained_genes": [row["gene"] for row in recovered_retained],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workbook", required=True, type=Path)
    parser.add_argument("--elements", required=True, type=Path)
    parser.add_argument("--chain", required=True, type=Path)
    parser.add_argument("--consensus", required=True, type=Path)
    parser.add_argument("--candidates", type=Path,
                        help="candidate-mutation table (explanation layer); without it the patient and variant counts are 0")
    parser.add_argument("--tissue", default="Liver-HCC")
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--summary", type=Path)
    args = parser.parse_args()

    try:
        from pyliftover import LiftOver
    except ImportError as error:
        raise SystemExit("pyliftover is required for hg19-to-hg38 comparison") from error

    table_rows = read_table5(args.workbook, args.tissue)
    if not table_rows:
        raise SystemExit(f"no Table 5 rows found for tissue {args.tissue!r}")
    rows = compare(table_rows, read_elements(args.elements), LiftOver(str(args.chain)),
                   read_tsv(args.consensus), read_tsv(args.candidates) if args.candidates else [])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    summary = summarize(rows)
    if args.summary:
        args.summary.parent.mkdir(parents=True, exist_ok=True)
        args.summary.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
