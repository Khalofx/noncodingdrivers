"""bin/postfilter_flags.py and bin/write_driver_evidence.py (D33): flags, q_postfilter and the schema."""
import csv
import json
import os
import re
import subprocess
import sys

import pandas as pd
from jsonschema import Draft202012Validator

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, "tests", "data", "driver_evidence")
SCHEMA = os.path.join(ROOT, "assets", "schema_driver_evidence.json")
SOURCES = os.path.join(ROOT, "assets", "evidence_sources.tsv")


def _run(tmp, mode="pcawg", scope="nominal", extra=None, out="test.driver_evidence.tsv"):
    pf = tmp / "test.postfilter.tsv"
    subprocess.run([sys.executable, os.path.join(ROOT, "bin", "postfilter_flags.py"), "--consensus", os.path.join(D, "test.consensus.tsv"),
                    "--element-flags", os.path.join(D, "test.element_flags.tsv"), "--mode", mode, "--output", str(pf)], check=True)
    ev = tmp / out
    subprocess.run([sys.executable, os.path.join(ROOT, "bin", "write_driver_evidence.py"), "--cohort", "test",
                    "--consensus", os.path.join(D, "test.consensus.tsv"), "--element-flags", os.path.join(D, "test.element_flags.tsv"),
                    "--postfilter", str(pf), "--sources", SOURCES, "--coordinates", "lifted", "--pipeline-version", "1.0.0dev",
                    "--scope", scope, "--output", str(ev)] + (extra or []), check=True)
    return pd.read_csv(pf, sep="\t", dtype=str).set_index("element_id"), ev


def test_postfilter_flags_and_q_postfilter(tmp_path):
    pf, _ = _run(tmp_path)
    assert pf.loc["promoter_A", "postfilter_status"] == "pass"
    assert pf.loc["enhancer_B", "mappability_status"] == "fail"          # 6 of 7 mutations in low mappability
    assert pf.loc["enhancer_C", "min_support_status"] == "fail"          # 2 mutations, 2 patients
    assert pf.loc["CTCF_D", "max_sample_fraction_status"] == "fail"      # one patient
    assert pf.loc["promoter_E", "mappability_status"] == "not_applied"   # no mutations
    assert pf.loc["promoter_E", "postfilter_status"] == "fail"           # min support fails at 0
    # failing elements get p = 1 before Benjamini-Hochberg; consensus_q is untouched
    assert float(pf.loc["enhancer_B", "q_postfilter"]) == 1.0
    assert float(pf.loc["promoter_A", "q_postfilter"]) < 1e-7


def test_flag_only_has_no_q_postfilter(tmp_path):
    pf, _ = _run(tmp_path, mode="flag_only")
    assert "q_postfilter" not in pf.columns


def test_evidence_table_validates_and_is_scoped(tmp_path):
    _, ev = _run(tmp_path)
    rows = list(csv.DictReader(open(ev), delimiter="\t"))
    errors = list(Draft202012Validator(json.load(open(SCHEMA))).iter_errors(rows))
    assert not errors, errors[0].message
    ids = {r["element_id"] for r in rows}
    assert ids == {"promoter_A", "enhancer_B", "enhancer_C", "promoter_E"}   # nominal scope: min p < 0.01 or consensus_q < 0.25
    seg = next(r for r in rows if r["element_id"] == "enhancer_B" and r["evidence_item"] == "frac_in_segdup")
    assert (seg["value"], seg["status"]) == ("1", "flagged")
    cf = next(r for r in rows if r["element_id"] == "enhancer_B" and r["evidence_item"] == "callable_fraction")
    assert cf["native_or_lifted"] == "native"
    assert all(r["native_or_lifted"] == "lifted" for r in rows if r["evidence_item"] != "callable_fraction")
    assert not any(re.search(r"tier", r["evidence_item"]) or re.search(r"tier|driver", r["status"]) for r in rows)   # no tier, no "driver" label
    _, ev_all = _run(tmp_path, scope="all")
    assert {r["element_id"] for r in csv.DictReader(open(ev_all), delimiter="\t")} == set(pd.read_csv(os.path.join(D, "test.consensus.tsv"), sep="\t").element_id)


def test_source_versions_match_module_pins():
    src = {r["evidence_item"]: r for r in csv.DictReader(open(SOURCES), delimiter="\t")}
    text = lambda p: open(os.path.join(ROOT, p)).read()   # noqa: E731
    assert "noncodingdrivers-activedriverwgs:" + src["p_activedriverwgs"]["source_version"] in text("modules/local/activedriverwgs/main.nf")
    assert src["p_oncodrivefml"]["source_version"] in text("modules/local/oncodrivefml/main.nf")
    assert "noncodingdrivers-oncodriveclustl:" + src["p_oncodriveclustl"]["source_version"] in text("modules/local/oncodriveclustl/main.nf")
    assert "bedtools=2.31.1" in text("modules/nf-core/bedtools/intersect/environment.yml")
    assert "echtvar=0.2.4" in text("modules/nf-core/echtvar/anno/environment.yml") and "echtvar 0.2.4" in src["frac_germline_pop_common"]["source_version"]
    assert "empiricalbrownsmethod=1.34.0" in text("modules/local/combinepvalues/environment.yml").lower()


def test_default_consensus_callers_leave_the_table_unchanged(tmp_path):
    # --consensus-callers naming every caller of the table (the default) writes no role rows: same bytes
    _, ev = _run(tmp_path)
    _, ev_named = _run(tmp_path, extra=["--consensus-callers", "activedriverwgs,oncodrivefml"], out="named.tsv")
    assert open(ev, "rb").read() == open(ev_named, "rb").read()
    assert not any(r["evidence_item"].startswith("caller_role_") for r in csv.DictReader(open(ev), delimiter="\t"))


def test_caller_reported_alongside_gets_role_rows_and_its_own_items(tmp_path):
    _, ev = _run(tmp_path, extra=["--consensus-callers", "activedriverwgs"], out="alongside.tsv")
    rows = list(csv.DictReader(open(ev), delimiter="\t"))
    errors = list(Draft202012Validator(json.load(open(SCHEMA))).iter_errors(rows))
    assert not errors, errors[0].message
    roles = {r["evidence_item"]: (r["element_id"], r["value"], r["unit"], r["status"]) for r in rows if r["evidence_item"].startswith("caller_role_")}
    assert roles == {"caller_role_activedriverwgs": (".", "1", "status", "combined"),
                     "caller_role_oncodrivefml": (".", "0", "status", "reported_alongside")}
    # the alongside caller keeps its p and q rows
    assert any(r["element_id"] == "promoter_A" and r["evidence_item"] == "q_oncodrivefml" for r in rows)


def test_oncodriveclustl_alongside_items_validate(tmp_path):
    # OncodriveCLUSTL reported alongside: its p/q rows, its cluster counts and the role rows; schema-valid
    ev = tmp_path / "clustl.driver_evidence.tsv"
    subprocess.run([sys.executable, os.path.join(ROOT, "bin", "write_driver_evidence.py"), "--cohort", "test",
                    "--consensus", os.path.join(D, "test.consensus_clustl.tsv"), "--sources", SOURCES, "--coordinates", "native",
                    "--pipeline-version", "1.0.0dev", "--consensus-callers", "activedriverwgs,oncodrivefml",
                    "--oncodriveclustl", os.path.join(D, "test.oncodriveclustl.tsv"), "--output", str(ev)], check=True)
    rows = list(csv.DictReader(open(ev), delimiter="\t"))
    errors = list(Draft202012Validator(json.load(open(SCHEMA))).iter_errors(rows))
    assert not errors, errors[0].message
    get = lambda e, i: next(r for r in rows if r["element_id"] == e and r["evidence_item"] == i)   # noqa: E731
    assert (get("promoter_A", "p_oncodriveclustl")["value"], get("promoter_A", "q_oncodriveclustl")["status"]) == ("0.0004", "significant")
    assert (get("promoter_A", "oncodriveclustl_n_clusters")["value"], get("promoter_A", "oncodriveclustl_clustered_mutations")["value"]) == ("2", "8")
    assert (get("enhancer_C", "oncodriveclustl_n_clusters")["value"], get("enhancer_C", "oncodriveclustl_n_clusters")["status"]) == ("0", "observed")
    assert (get("promoter_E", "oncodriveclustl_n_clusters")["value"], get("promoter_E", "oncodriveclustl_n_clusters")["status"]) == ("NA", "not_tested")
    assert get(".", "caller_role_oncodriveclustl")["status"] == "reported_alongside"
    assert get(".", "caller_role_activedriverwgs")["status"] == "combined"
    # the caller items keep the fixed order: ActiveDriverWGS, OncodriveFML, then OncodriveCLUSTL
    items = [r["evidence_item"] for r in rows if r["element_id"] == "promoter_A"]
    assert items.index("q_oncodrivefml") < items.index("p_oncodriveclustl") < items.index("oncodriveclustl_n_clusters")


LAYER_PREFIXES = ("avi_", "opencravat_", "cosmic_", "cancer_gene_", "re2g_", "ccre_v4_", "null_fpr_")
LAYER_ARGS = ["--avi-summary", os.path.join(D, "test.avi.element_summary.tsv"), "--opencravat", os.path.join(D, "test.opencravat.tsv"),
              "--cosmic-elements", os.path.join(D, "test.cosmic.element_summary.tsv"),
              "--cancer-genes", os.path.join(D, "test.cancer_gene_context.tsv"),
              "--re2g-summary", os.path.join(D, "test.enhancer_gene_summary.tsv"), "--ccre-summary", os.path.join(D, "test.ccre_v4_summary.tsv"),
              "--null-fpr", os.path.join(D, "test.null_fpr.tsv")]


def test_layer_summary_items_validate(tmp_path):
    # every layer on (synthetic tables): one fixed set of summary items per significant element and layer
    _, ev = _run(tmp_path, extra=LAYER_ARGS, out="layers.tsv")
    rows = list(csv.DictReader(open(ev), delimiter="\t"))
    errors = list(Draft202012Validator(json.load(open(SCHEMA))).iter_errors(rows))
    assert not errors, errors[0].message
    get = lambda e, i: next(r for r in rows if r["element_id"] == e and r["evidence_item"] == i)   # noqa: E731
    layer = [r for r in rows if r["evidence_item"].startswith(LAYER_PREFIXES) and r["element_id"] != "."]
    # 16 items per significant element of the nominal scope (CTCF_D is outside it, its AVI row is not read)
    assert {r["element_id"] for r in layer} == {"promoter_A", "enhancer_B", "enhancer_C", "promoter_E"}
    assert len(layer) == 4 * 16
    assert (get("promoter_A", "avi_max_phred")["value"], get("promoter_A", "avi_max_phred")["unit"]) == ("31.2", "score")
    assert (get("promoter_A", "avi_n_scored")["value"], get("promoter_A", "avi_n_unscored")["value"]) == ("2", "1")
    assert (get("enhancer_B", "avi_max_phred")["value"], get("enhancer_B", "avi_max_phred")["status"]) == ("NA", "not_applied")
    assert get("enhancer_C", "avi_n_scored")["status"] == "not_applied"            # no row in the AVI summary
    # OpenCRAVAT: candidate variants with an annotation per annotator; a variant shared by two elements counts in both
    assert [get("promoter_A", f"opencravat_n_{a}")["value"] for a in ("ccre_screen", "repeat", "clinvar")] == ["2", "1", "1"]
    assert [get("enhancer_B", f"opencravat_n_{a}")["value"] for a in ("ccre_screen", "repeat", "clinvar")] == ["1", "1", "1"]
    assert get("promoter_E", "opencravat_n_clinvar")["value"] == "0"
    # COSMIC: counts only, release in source_version
    assert [get("promoter_A", i)["value"] for i in ("cosmic_n_matched_mutations", "cosmic_n_samples", "cosmic_n_sample_name_overlaps")] == ["2", "150", "1"]
    assert get("promoter_A", "cosmic_n_samples")["source_version"] == "COSMIC v999"
    assert all(re.fullmatch(r"NA|[0-9]+", r["value"]) for r in rows if r["evidence_item"].startswith("cosmic_"))
    assert not any(re.search(r"COS[VMS][0-9]", "\t".join(r.values())) for r in rows)
    # cancer gene list labels
    assert (get("promoter_A", "cancer_gene_nearest_in_list")["value"], get("promoter_A", "cancer_gene_nearest_in_list")["status"]) == ("1", "in_list")
    assert get("enhancer_C", "cancer_gene_nearest_in_list")["status"] == "not_applied"   # no nearest gene
    assert get("promoter_E", "cancer_gene_re2g_targets_in_list")["value"] == "2"
    assert get("promoter_E", "cancer_gene_nearest_pc_in_list")["source_version"].startswith("synthetic list (md5 0123456789abcdef")
    # ENCODE-rE2G over two tissues: union of distal targets, top by score, nearest status over the tissues
    assert (get("enhancer_B", "re2g_n_distal_targets")["value"], get("enhancer_B", "re2g_top_distal_target")["value"]) == ("3", "GENEZ")
    assert get("enhancer_B", "re2g_top_distal_target")["status"] == "nearest_in_targets"
    assert (get("promoter_A", "re2g_top_distal_target")["value"], get("promoter_A", "re2g_top_distal_target")["status"]) == ("NA", "no_link")
    assert get("enhancer_C", "re2g_top_distal_target")["status"] == "no_nearest_gene"
    assert get("promoter_E", "re2g_n_distal_targets")["status"] == "not_applied"
    assert get("enhancer_B", "re2g_n_distal_targets")["source_version"].endswith("tissues liver,heart")
    # SCREEN V4
    assert (get("promoter_E", "ccre_v4_dominant_class")["value"], get("promoter_E", "ccre_v4_covered_fraction")["value"]) == ("CA-H3K4me3", "1")
    assert get("enhancer_B", "ccre_v4_dominant_class")["value"] == "none"
    # element-level summaries of GRCh38 resources are native; mutation-based ones carry the cohort's provenance
    assert get("promoter_A", "ccre_v4_covered_fraction")["native_or_lifted"] == "native"
    assert get("promoter_A", "avi_n_scored")["native_or_lifted"] == "lifted"
    # OncodriveCLUSTL null false-positive rate (cohort level)
    assert (get(".", "null_fpr_oncodriveclustl")["value"], get(".", "null_fpr_oncodriveclustl")["threshold"]) == ("0.046875", "0.05")
    # layer items follow the element's other items, in the fixed layer order
    items = [r["evidence_item"] for r in rows if r["element_id"] == "promoter_A"]
    assert items.index("postfilter_status") < items.index("avi_max_phred") < items.index("opencravat_n_ccre_screen") \
        < items.index("cosmic_n_matched_mutations") < items.index("cancer_gene_nearest_in_list") < items.index("re2g_n_distal_targets") \
        < items.index("ccre_v4_dominant_class")


def test_layers_add_rows_only(tmp_path):
    # with the layer rows removed, the layers-on table is the layers-off table, byte for byte; scope all
    # adds no layer row for an element that is not significant
    for scope in ("nominal", "all"):
        _, off = _run(tmp_path, scope=scope, out=f"off.{scope}.tsv")
        _, on = _run(tmp_path, scope=scope, extra=LAYER_ARGS, out=f"on.{scope}.tsv")
        kept = [line for line in open(on) if not line.split("\t")[2].startswith(LAYER_PREFIXES)]
        assert "".join(kept) == open(off).read()
        sig = {r["element_id"] for r in csv.DictReader(open(on), delimiter="\t") if r["evidence_item"].startswith("ccre_v4_")}
        assert sig == {"promoter_A", "enhancer_B", "enhancer_C", "promoter_E"}


def test_label_values_only_for_label_items(tmp_path):
    # the schema keeps numeric values for every other item
    _, ev = _run(tmp_path, extra=LAYER_ARGS, out="labels.tsv")
    rows = list(csv.DictReader(open(ev), delimiter="\t"))
    bad = dict(next(r for r in rows if r["evidence_item"] == "ccre_v4_covered_fraction"), value="dELS")
    assert list(Draft202012Validator(json.load(open(SCHEMA))).iter_errors([bad]))
    src = {r["evidence_item"]: r for r in csv.DictReader(open(SOURCES), delimiter="\t")}
    assert {i for i, r in src.items() if r["unit"] == "label"} == {"re2g_top_distal_target", "ccre_v4_dominant_class"}
