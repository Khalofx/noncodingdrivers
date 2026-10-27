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


def _run(tmp, mode="pcawg", scope="nominal"):
    pf = tmp / "test.postfilter.tsv"
    subprocess.run([sys.executable, os.path.join(ROOT, "bin", "postfilter_flags.py"), "--consensus", os.path.join(D, "test.consensus.tsv"),
                    "--element-flags", os.path.join(D, "test.element_flags.tsv"), "--mode", mode, "--output", str(pf)], check=True)
    ev = tmp / "test.driver_evidence.tsv"
    subprocess.run([sys.executable, os.path.join(ROOT, "bin", "write_driver_evidence.py"), "--cohort", "test",
                    "--consensus", os.path.join(D, "test.consensus.tsv"), "--element-flags", os.path.join(D, "test.element_flags.tsv"),
                    "--postfilter", str(pf), "--sources", SOURCES, "--coordinates", "lifted", "--pipeline-version", "1.0.0dev",
                    "--scope", scope, "--output", str(ev)], check=True)
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
    assert "bedtools=2.31.1" in text("modules/nf-core/bedtools/intersect/environment.yml")
    assert "echtvar=0.2.4" in text("modules/nf-core/echtvar/anno/environment.yml") and "echtvar 0.2.4" in src["frac_germline_pop_common"]["source_version"]
    assert "empiricalbrownsmethod=1.34.0" in text("modules/local/combinepvalues/environment.yml").lower()
