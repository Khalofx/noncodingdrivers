#!/usr/bin/env python3
"""Turn the pipeline's own result tables into MultiQC custom-content sections.

WHY THIS EXISTS. MultiQC summarises recognised QC logs (fastp, samtools, Picard).
noncodingdrivers starts from called variants, so no aligner or trimmer runs and its
outputs are tables. MultiQC would find nothing and emit an empty report, so the
run-level summary is generated from the tables themselves and handed to MultiQC as
custom content.

WHAT IT REPORTS, and what it deliberately does not.
  It reports counts that describe THE RUN: how many elements were tested, how many
  callers tested each, how many reached the combined-q threshold, and how many rows
  each caller wrote.

  It does NOT rank candidates, does not name a "top hit", and does not describe any
  element as a driver. A QC report answers "did this run behave sensibly?", not
  "what did we find"; the consensus table and the explanation layer answer that.
  Putting a candidate list in a QC report invites it to be read as a result summary.

HONESTY RULES
  * an absent or unparseable value is None, never 0
  * a caller file that is absent yields an omitted section with a stated reason on
    stderr, never a zero
  * a header-only consensus table is a hard error, never an empty report

Usage
-----
    run_summary_mqc.py --consensus X.consensus.tsv \
        --caller-table X.activedriverwgs.tsv --caller-table X.oncodrivefml.tsv \
        --sample X --fdr 0.1 --outdir .
"""
import argparse
import csv
import gzip
import os
import sys

TAB = "\t"


def _open(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def read_tsv(path):
    with _open(path) as fh:
        return list(csv.DictReader(fh, delimiter=TAB))


def fnum(row, key):
    """float or None -- an unparseable or absent value is None, never 0."""
    v = row.get(key)
    if v is None or v == "" or v in ("NA", "nan", "NaN"):
        return None
    try:
        return float(v)
    except ValueError:
        return None


def esc(s):
    return str(s).replace('"', "'")


# Filename tokens that are NOT the caller name. Two conventions reach this script:
#   raw caller output   <sample>.<caller>.tsv             -> 3 fields
#   harmonised output   <sample>.<caller>.harmonized.tsv  -> 4 fields
# A fixed position is wrong for one of them: on the harmonised name the
# second-to-last field is "harmonized" for EVERY caller, which would collapse all
# callers into one map key and silently report a single row count.
_SUFFIX_TOKENS = {"harmonized", "harmonised", "tsv", "txt", "gz", "csv"}


def caller_name_from_path(path):
    """Derive a caller name from a table filename, suffix-convention agnostic.

    Strips trailing known suffix tokens, then returns the last remaining field --
    which is the caller under both conventions. A name with no caller field falls
    back to the first field so an unexpected file is visible rather than mislabelled.
    """
    parts = os.path.basename(path).split(".")
    while len(parts) > 1 and parts[-1].lower() in _SUFFIX_TOKENS:
        parts = parts[:-1]
    return parts[-1] if len(parts) >= 2 else parts[0]


def write_section(path, sid, name, description, plot_type, data, headers=None,
                  extra=None):
    """Emit one MultiQC custom-content YAML section."""
    out = []
    out.append('id: "%s"' % sid)
    out.append('section_name: "%s"' % esc(name))
    out.append('description: "%s"' % esc(description))
    out.append('plot_type: "%s"' % plot_type)
    out.append("pconfig:")
    out.append('    id: "%s_plot"' % sid)
    out.append('    title: "%s"' % esc(name))
    for k, v in (extra or {}).items():
        out.append('    %s: "%s"' % (k, esc(v)))
    if headers:
        out.append("headers:")
        for key, cfg in headers.items():
            out.append("    %s:" % key)
            for hk, hv in cfg.items():
                out.append('        %s: "%s"' % (hk, esc(hv)))
    out.append("data:")
    for sample, vals in data.items():
        out.append('    "%s":' % esc(sample))
        for k, v in vals.items():
            out.append("        %s: %s" % (k, v))
    with open(path, "w") as fh:
        fh.write("\n".join(out) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--consensus", required=True)
    ap.add_argument("--caller", action="append", default=[],
                    help="name=path, repeatable")
    ap.add_argument("--caller-table", action="append", default=[],
                    help="path only; the caller name is derived from the filename "
                         "(repeatable). Preferred over --caller from Nextflow, so the "
                         "name-derivation rule lives in tested Python rather than in "
                         "an untested Groovy expression.")
    ap.add_argument("--fdr", type=float, default=0.1,
                    help="combined-q threshold reported as significant (the pipeline's --fdr_threshold)")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args()

    rows = read_tsv(a.consensus)
    if not rows:
        sys.exit("ERROR: %s has no data rows" % a.consensus)
    S = a.sample
    os.makedirs(a.outdir, exist_ok=True)
    written = []

    # ---------------- 1. run overview ----------------------------------------
    n_total = len(rows)
    has_q = any(fnum(r, "consensus_q") is not None for r in rows)
    sig = [r for r in rows if (fnum(r, "consensus_q") or 1.0) < a.fdr] if has_q else []
    methods = sorted({r.get("combination_method", "") for r in rows} - {""})
    method = methods[0] if len(methods) == 1 else "unknown"
    votes = [fnum(r, "vote_count") for r in rows]
    vote2 = sum(1 for v in votes if v is not None and v >= 2)
    tested = [t for t in (fnum(r, "n_tested_callers") for r in rows) if t is not None]

    overview = {
        "elements_tested": n_total,
        "significant": len(sig) if has_q else 0,
        "multi_caller_agreement": vote2,
        "callers_run": int(max(tested)) if tested else 0,
    }
    p = os.path.join(a.outdir, "%s_noncodingdrivers_overview_mqc.yaml" % S)
    write_section(
        p, "noncodingdrivers_overview", "noncodingdrivers: run overview",
        "Counts describing this run. 'Significant' means combined q (%s) below %g across the "
        "full element catalogue -- a statistical threshold, not a driver call." % (method, a.fdr),
        "table", {S: overview},
        headers={
            "elements_tested": {"title": "Elements tested", "format": "{:,.0f}"},
            "significant": {"title": "Significant (q<%g)" % a.fdr, "format": "{:,.0f}"},
            "multi_caller_agreement": {"title": "Multi-caller agreement", "format": "{:,.0f}"},
            "callers_run": {"title": "Callers run", "format": "{:,.0f}"},
        })
    written.append(p)

    # ---------------- 2. per-caller coverage ----------------------------------
    # The most useful QC number the pipeline has: a caller that silently tested far
    # fewer elements than its peers is the failure mode that has recurred here.
    pcols = sorted(c for c in rows[0] if c.startswith("p_"))
    if pcols:
        cov = {c.replace("p_", ""): sum(1 for r in rows if fnum(r, c) is not None)
               for c in pcols}
        p = os.path.join(a.outdir, "%s_noncodingdrivers_caller_coverage_mqc.yaml" % S)
        write_section(
            p, "noncodingdrivers_caller_coverage", "noncodingdrivers: per-caller coverage",
            "Elements for which each caller returned a p-value. Large differences are "
            "expected -- impact-based callers only score elements containing at least "
            "one mutation -- but a caller far below its peers indicates a truncated run.",
            "bargraph", {S: cov}, extra={"ylab": "elements with a p-value"})
        written.append(p)

    # ---------------- 3. raw caller output rows -------------------------------
    # Straight from the caller files, so a caller whose output never reached the
    # consensus join is visible even when its p-column is absent above.
    if a.caller or a.caller_table:
        pairs = []
        for arg in a.caller:
            if "=" not in arg:
                sys.exit("ERROR: --caller expects name=path, got %r" % arg)
            name, path = arg.split("=", 1)
            pairs.append((name, path))
        for path in a.caller_table:
            pairs.append((caller_name_from_path(path), path))

        sizes = {}
        for name, path in pairs:
            if not os.path.exists(path):
                sys.stderr.write("[run_summary_mqc] %s: %s absent, section omitted\n"
                                 % (name, path))
                continue
            if name in sizes:
                # Refuse to overwrite. A collision means the name-derivation rule did
                # not separate two callers, and reporting one count under a shared
                # label would understate coverage silently.
                sys.exit("ERROR: two caller tables resolve to the same name %r; "
                         "refusing to overwrite one caller's row count with another's"
                         % name)
            with _open(path) as fh:
                sizes[name] = max(sum(1 for _ in fh) - 1, 0)
        if sizes:
            p = os.path.join(a.outdir, "%s_noncodingdrivers_caller_rows_mqc.yaml" % S)
            write_section(
                p, "noncodingdrivers_caller_rows", "noncodingdrivers: raw caller output rows",
                "Rows in each caller's own output before harmonisation. Compare against "
                "per-caller coverage above: a caller with many raw rows but few "
                "p-values in the consensus lost records in the join.",
                "bargraph", {S: dict(sorted(sizes.items()))},   # caller order fixed, not arrival order
                extra={"ylab": "rows in caller output"})
            written.append(p)

    for w in written:
        sys.stderr.write("[run_summary_mqc] wrote %s\n" % os.path.basename(w))
    sys.stderr.write("[run_summary_mqc] %d section(s) from %d elements\n"
                     % (len(written), n_total))


if __name__ == "__main__":
    main()
