"""The reference datasheet (assets/references/GRCh38.yml) against its schema and the provenance record.

Every checksum in the datasheet must also be recorded in tests/regression/PROVENANCE.md, and the
assets that ship with the pipeline must match their datasheet checksums. With
NONCODINGDRIVERS_REFERENCE_BASE set (the Karakoram data/refs folder), the staged files are
checked as well.
"""
import csv
import gzip
import hashlib
import json
import os
import re

import pytest
import yaml
from jsonschema import Draft202012Validator

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SHEET = os.path.join(ROOT, "assets", "references", "GRCh38.yml")
SCHEMA = os.path.join(ROOT, "assets", "schema_references.json")
PROVENANCE = os.path.join(ROOT, "tests", "regression", "PROVENANCE.md")
PLACEHOLDER = "${reference_base}"
# Checksums of a user-downloaded asset that has not been inspected yet (AlphaGenome AVI, to-do S7A.15)
PENDING = "pending-user-download"
PENDING_ALLOWED = {"alphagenome_avi_md5", "alphagenome_avi_tbi_md5", "alphagenome_avi_archive_md5"}
MANIFEST = os.path.join(ROOT, "assets", "references", "encode_re2g.GRCh38.tsv")


def _entry():
    data = yaml.safe_load(open(SHEET))
    assert isinstance(data, list) and len(data) == 1, "the datasheet is a one-entry list"
    return data[0]


def _md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def test_datasheet_validates_against_schema():
    schema = json.load(open(SCHEMA))
    Draft202012Validator.check_schema(schema)
    errors = sorted(Draft202012Validator(schema).iter_errors(yaml.safe_load(open(SHEET))), key=str)
    assert not errors, "\n".join(e.message for e in errors)


def test_every_property_is_a_meta_key():
    """samplesheetToList returns the whole entry as one meta map only if every field is a meta key."""
    props = json.load(open(SCHEMA))["items"]["properties"]
    assert all(v.get("meta") == [k] for k, v in props.items())


def test_paths_are_relative_to_reference_base_or_the_pipeline():
    path_keys = ["fasta", "fasta_fai", "fasta_dict", "element_catalogue", "gtf", "vep_cache", "oncodrivefml_genome",
                 "regulatory_peaks_dir", "liftover_chain_from_hg19", "liftover_source_fasta", "liftover_source_fasta_fai",
                 "blacklist_bed", "segdup_bed", "low_mappability_bed", "umap_bed", "gnomad_echtvar", "alphagenome_avi",
                 "enhancer_gene_dir", "ccre_v4_bed", "ccre_celltype_dir"]
    e = _entry()
    for k in path_keys:
        assert e[k].startswith(PLACEHOLDER + "/"), k
    for k in ("rename_map", "jaspar_motifs"):
        assert e[k].startswith("assets/") and os.path.exists(os.path.join(ROOT, e[k])), k


def test_checksums_are_recorded_in_provenance():
    recorded = set(re.findall(r"\b[0-9a-f]{32}\b", open(PROVENANCE).read()))
    e = _entry()
    missing = {k: v for k, v in e.items() if k.endswith("_md5") and v not in recorded and v != PENDING}
    assert not missing, f"datasheet checksums not in PROVENANCE.md: {missing}"


def test_pending_checksums_are_only_the_user_downloaded_avi_file():
    """A placeholder checksum is allowed for the AlphaGenome AVI download only, never for an asset the pipeline validated."""
    e = _entry()
    pending = {k for k, v in e.items() if v == PENDING}
    assert pending <= PENDING_ALLOWED, f"placeholder checksums outside the AVI keys: {sorted(pending - PENDING_ALLOWED)}"
    assert "never bundled" in e["alphagenome_avi_licence"] and "permissive-use" in e["alphagenome_avi_licence"]


def test_bundled_assets_match_their_checksums():
    e = _entry()
    assert _md5(os.path.join(ROOT, e["jaspar_motifs"])) == e["jaspar_md5"]


def test_rename_map_targets_ucsc_names():
    e = _entry()
    rows = [line.rstrip("\n").split("\t") for line in open(os.path.join(ROOT, e["rename_map"]))]
    assert all(len(r) == 2 for r in rows), "two tab-separated fields per line, no comments"
    assert all(r[1].startswith("chr") for r in rows)
    sources = [r[0] for r in rows]
    assert len(sources) == len(set(sources)), "each source name maps once"
    mapping = dict(rows)
    assert mapping["1"] == "chr1" and mapping["MT"] == "chrM" and mapping["NC_000001.11"] == "chr1"


@pytest.mark.skipif(not os.environ.get("NONCODINGDRIVERS_REFERENCE_BASE"), reason="staged references not available")
def test_staged_references_match_the_datasheet():
    base = os.environ["NONCODINGDRIVERS_REFERENCE_BASE"].rstrip("/")
    e = _entry()
    pairs = {"fasta": "fasta_md5", "fasta_fai": "fasta_fai_md5", "fasta_dict": "fasta_dict_md5",
             "element_catalogue": "element_catalogue_md5", "gtf": "gtf_md5", "liftover_chain_from_hg19": "liftover_chain_md5",
             "liftover_source_fasta": "liftover_source_fasta_md5", "liftover_source_fasta_fai": "liftover_source_fasta_fai_md5",
             "blacklist_bed": "blacklist_md5", "segdup_bed": "segdup_md5", "low_mappability_bed": "low_mappability_md5",
             "umap_bed": "umap_md5", "gnomad_echtvar": "gnomad_echtvar_md5", "ccre_v4_bed": "ccre_v4_md5"}
    # the AVI file once its checksums are recorded (placeholders until the user's download, S7A.15)
    for path_key, md5_key in (("alphagenome_avi", "alphagenome_avi_md5"),):
        if e[md5_key] != PENDING:
            pairs[path_key] = md5_key
    bad = {}
    for path_key, md5_key in pairs.items():
        path = e[path_key].replace(PLACEHOLDER, base)
        if _md5(path) != e[md5_key]:
            bad[path_key] = path
    assert not bad, f"checksum differs: {bad}"
    rows = sum(1 for _ in open(e["element_catalogue"].replace(PLACEHOLDER, base)))
    assert rows == e["element_catalogue_rows"]
    assert sum(1 for _ in open(e["fasta_fai"].replace(PLACEHOLDER, base))) == e["contig_count"]
    with gzip.open(e["ccre_v4_bed"].replace(PLACEHOLDER, base), "rt") as fh:
        assert sum(1 for _ in fh) == e["ccre_v4_rows"]
    # every staged ENCODE-rE2G and SCREEN V4 file matches its PROVENANCE.tsv row and the manifest
    manifest = {r["accession"]: r for r in csv.DictReader(open(MANIFEST), delimiter="\t")}
    for folder in ("encode_re2g", "screen_v4"):
        for row in csv.DictReader(open(os.path.join(base, folder, "PROVENANCE.tsv")), delimiter="\t"):
            if row["layer"] == "re2g":
                path = os.path.join(base, folder, row["tissue"], manifest[row["accession"]]["file_name"])
                assert manifest[row["accession"]]["md5"] == row["md5"]
            elif row["layer"] == "ccre_v4":
                path = e["ccre_v4_bed"].replace(PLACEHOLDER, base)
            else:
                path = [os.path.join(base, folder, "celltype", row["tissue"], f) for f in os.listdir(os.path.join(base, folder, "celltype", row["tissue"]))
                        if row["accession"] in f][0]
            assert _md5(path) == row["md5"], path


def test_encode_re2g_manifest():
    """Coverage table of ENCODE-rE2G per tissue folder: DNaseOnly thresholded links, unique accessions,
    liver = the 11-file D-E1 set, every other tissue its normal tissue samples."""
    rows = list(csv.DictReader(open(MANIFEST), delimiter="\t"))
    real = [r for r in rows if re.match(r"^ENCFF\w+$", r["accession"])]
    assert len({r["accession"] for r in real}) == len(real)
    assert all(re.match(r"^[0-9a-f]{32}$", r["md5"]) and int(r["bytes"]) > 0 for r in real)
    assert all(r["file_name"] == f"{r['file_name'].split('.')[0]}.{r['accession']}.bed.gz" for r in real)
    assert {r["model"] for r in real} == {"DNaseOnly"}
    selected = [r for r in real if r["selected"] == "yes"]
    liver = {r["accession"] for r in selected if r["tissue"] == "liver"}
    assert liver == {"ENCFF960GTR", "ENCFF775WBY", "ENCFF239HUO", "ENCFF167GAJ", "ENCFF184TKS", "ENCFF441ZFG",
                     "ENCFF284LMY", "ENCFF955UNH", "ENCFF043RAT", "ENCFF429OPN", "ENCFF603ILE"}
    assert all(r["biosample_classification"] == "tissue" for r in selected if r["tissue"] != "liver")
    assert len(selected) == 407
