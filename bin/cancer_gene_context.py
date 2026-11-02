#!/usr/bin/env python3
"""Cancer gene-list labels of the significant elements (CANCERGENECONTEXT).

Labels three gene assignments of every significant element against a cancer gene list (the
COSMIC Cancer Gene Census or any user list). Annotation only: no element is filtered, ranked or
weighted, and the nearest-gene table itself is not changed.

  * the nearest gene (any biotype) of EXPLAIN's <cohort>.nearest_gene.tsv;
  * the nearest protein-coding gene: bedtools closest -d -t first from the element midpoint to
    the protein_coding TSSs of the same TSS BED (an unnamed lncRNA is often the nearest gene);
  * the ENCODE-rE2G target genes of <cohort>.enhancer_gene_links.tsv, when that layer ran.

Gene list layouts (tab- or comma-separated, optionally gzipped, CRLF tolerated):
  * Census: GENE_SYMBOL, TIER, ROLE_IN_CANCER (and SYNONYMS, whose ENSG ids are used);
  * plain: a gene_symbol column, optional ensembl_gene_id, tier, role columns;
  * plain without a header: one gene symbol per line.
Matching: the Ensembl gene id without version when the list has one for the gene, otherwise
the exact, case-sensitive gene symbol.
"""

import argparse
import csv
import gzip
import hashlib
import io
import re
import sys
from pathlib import Path

ENSG = re.compile(r"^ENSG\d+(\.\d+)?$")
SYMBOL_COLUMNS = ("GENE_SYMBOL", "gene_symbol")
ENSEMBL_COLUMNS = ("ensembl_gene_id", "ENSEMBL_GENE_ID")
TIER_COLUMNS = ("TIER", "tier")
ROLE_COLUMNS = ("ROLE_IN_CANCER", "role")
NEAREST_REQUIRED = {"element_id", "nearest_gene", "nearest_gene_ensembl_id", "gene_type"}
LINKS_REQUIRED = {"element_id", "target_gene", "target_gene_ensembl_id", "is_self_promoter"}

FIELDS = [
    "element_id", "nearest_gene", "nearest_gene_ensembl_id", "nearest_gene_type", "nearest_in_list",
    "nearest_match_by", "nearest_list_tier", "nearest_list_role", "nearest_protein_coding_gene",
    "nearest_protein_coding_gene_ensembl_id", "nearest_protein_coding_distance", "nearest_pc_in_list",
    "nearest_pc_match_by", "nearest_pc_list_tier", "nearest_pc_list_role", "re2g_targets_in_list",
    "re2g_distal_targets_in_list", "list_name", "list_md5", "list_genes",
]


def strip_version(gene_id):
    return (gene_id or "").split(".")[0]


def read_text(path):
    data = Path(path).read_bytes()
    md5 = hashlib.md5(data).hexdigest()
    if data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)
    text = data.decode("utf-8-sig", errors="replace")
    return text, md5


def first_present(header, names):
    return next((name for name in names if name in header), None)


class GeneList:
    def __init__(self, path):
        text, self.md5 = read_text(path)
        lines = [line for line in io.StringIO(text, newline=None).read().splitlines()
                 if line.strip() and not line.startswith("#")]
        if not lines:
            raise ValueError(f"cancer gene list {path} is empty")
        delimiter = "\t" if "\t" in lines[0] else ("," if "," in lines[0] else "\t")
        rows = list(csv.reader(lines, delimiter=delimiter))
        header = [h.strip() for h in rows[0]]
        symbol_col = first_present(header, SYMBOL_COLUMNS)
        self.by_symbol = {}
        self.by_id = {}
        if symbol_col is None:
            if any(len(row) != 1 for row in rows):
                raise ValueError(f"cancer gene list {path}: no {' or '.join(SYMBOL_COLUMNS)} column, and not a one-column list")
            self.layout = "plain_no_header"
            for row in rows:
                self._add(row[0].strip(), set(), "", "")
            return
        self.layout = "census" if symbol_col == "GENE_SYMBOL" else "plain"
        index = {name: i for i, name in enumerate(header)}
        ens_col = first_present(header, ENSEMBL_COLUMNS)
        tier_col = first_present(header, TIER_COLUMNS)
        role_col = first_present(header, ROLE_COLUMNS)
        syn_col = "SYNONYMS" if "SYNONYMS" in header else None

        def get(row, col):
            return row[index[col]].strip() if col and index[col] < len(row) else ""
        for row in rows[1:]:
            ids = set()
            for col in (ens_col, syn_col):
                for token in re.split(r"[,;| ]+", get(row, col)):
                    if ENSG.match(token):
                        ids.add(strip_version(token))
            self._add(get(row, symbol_col), ids, get(row, tier_col), get(row, role_col))

    def _add(self, symbol, ids, tier, role):
        if not symbol:
            return
        entry = {"symbol": symbol, "ids": ids, "tier": tier, "role": role}
        self.by_symbol.setdefault(symbol, entry)
        for gene_id in ids:
            self.by_id.setdefault(gene_id, entry)

    def __len__(self):
        return len(self.by_symbol)

    def match(self, symbol, gene_id):
        """(entry, 'ensembl_id' | 'symbol') or (None, '')."""
        gene_id = strip_version(gene_id)
        if gene_id and gene_id in self.by_id:
            return self.by_id[gene_id], "ensembl_id"
        entry = self.by_symbol.get(symbol or "")
        if entry and not (entry["ids"] and gene_id and gene_id not in entry["ids"]):
            return entry, "symbol"
        return None, ""


def read_tsv(path, required, label):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{label} lacks required columns: {sorted(missing)}")
        return list(reader)


def read_nearest_pc(path):
    """bedtools closest -d -t first: midpoint BED6+2, protein-coding TSS BED6+2, distance."""
    nearest = {}
    with open(path) as handle:
        for line in handle:
            p = line.rstrip("\n").split("\t")
            if len(p) < 17:
                raise ValueError(f"{path}: expected 17 columns (bedtools closest -d), got {len(p)}")
            if p[11] == "." or p[16] == "-1":
                nearest[p[3]] = ("", "", "")
            else:
                nearest[p[3]] = (p[14], p[11], p[16])
    return nearest


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--nearest", required=True, help="<cohort>.nearest_gene.tsv of EXPLAIN")
    parser.add_argument("--nearest-pc", default="", help="bedtools closest output against protein_coding TSSs")
    parser.add_argument("--links", default="", help="<cohort>.enhancer_gene_links.tsv (ENCODE-rE2G layer)")
    parser.add_argument("--gene-list", required=True)
    parser.add_argument("--list-name", default="")
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)

    genes = GeneList(args.gene_list)
    list_name = args.list_name or Path(args.gene_list).name
    nearest = read_tsv(args.nearest, NEAREST_REQUIRED, "nearest-gene table")
    nearest_pc = read_nearest_pc(args.nearest_pc) if args.nearest_pc else None
    links = read_tsv(args.links, LINKS_REQUIRED, "enhancer-gene links table") if args.links else None

    def label(symbol, gene_id):
        if not symbol and not gene_id:
            return "NA", "", "", ""
        entry, by = genes.match(symbol, gene_id)
        return ("yes", by, entry["tier"], entry["role"]) if entry else ("no", "", "", "")

    out = []
    for row in nearest:
        element_id = row["element_id"]
        in_list, by, tier, role = label(row["nearest_gene"], row["nearest_gene_ensembl_id"])
        record = {
            "element_id": element_id, "nearest_gene": row["nearest_gene"],
            "nearest_gene_ensembl_id": row["nearest_gene_ensembl_id"], "nearest_gene_type": row["gene_type"],
            "nearest_in_list": in_list, "nearest_match_by": by, "nearest_list_tier": tier, "nearest_list_role": role,
            "list_name": list_name, "list_md5": genes.md5, "list_genes": len(genes),
        }
        if nearest_pc is None:
            record.update({k: "NA" for k in ("nearest_protein_coding_gene", "nearest_protein_coding_gene_ensembl_id",
                                             "nearest_protein_coding_distance", "nearest_pc_in_list")})
        else:
            symbol, gene_id, distance = nearest_pc.get(element_id, ("", "", ""))
            pc_in, pc_by, pc_tier, pc_role = label(symbol, gene_id)
            record.update({"nearest_protein_coding_gene": symbol, "nearest_protein_coding_gene_ensembl_id": gene_id,
                           "nearest_protein_coding_distance": distance, "nearest_pc_in_list": pc_in,
                           "nearest_pc_match_by": pc_by, "nearest_pc_list_tier": pc_tier, "nearest_pc_list_role": pc_role})
        if links is None:
            record["re2g_targets_in_list"] = record["re2g_distal_targets_in_list"] = "NA"
        else:
            hits, distal = set(), set()
            for link in links:
                if link["element_id"] != element_id:
                    continue
                entry, _ = genes.match(link["target_gene"], link["target_gene_ensembl_id"])
                if entry:
                    hits.add(link["target_gene"])
                    if link["is_self_promoter"] != "TRUE":
                        distal.add(link["target_gene"])
            record["re2g_targets_in_list"] = ",".join(sorted(hits))
            record["re2g_distal_targets_in_list"] = ",".join(sorted(distal))
        out.append(record)

    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n", restval="")
        writer.writeheader()
        writer.writerows(out)
    sys.stderr.write(f"cancer_gene_context: {len(out)} elements; list {list_name} ({genes.layout}, {len(genes)} genes); "
                     f"nearest in list {sum(r['nearest_in_list'] == 'yes' for r in out)}\n")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, csv.Error) as error:
        sys.exit(f"ERROR: {error}")
