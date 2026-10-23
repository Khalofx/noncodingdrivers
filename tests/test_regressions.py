"""Regression checks of the kept bin/ scripts.

Checks of scripts that were replaced by packaged tools (the in-house p-value combination, motif
scoring, coding-sequence exclusion) left with those scripts. The p-value combination is covered
by the COMBINEPVALUES nf-test.
"""
import os

BIN = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "bin")


def test_activedriverwgs_runner_bisects_instead_of_dropping_shards():
    """Discarding a whole shard on any error silently drops a third of the catalogue from a
    cohort run, and untested elements can never be flagged. The runner must bisect a failing
    shard and record the element ids it could not test."""
    body = open(os.path.join(BIN, "run_activedriverwgs.R")).read()
    assert "bisect" in body.lower(), "shard failures must bisect, not drop"
    assert "untested" in body.lower(), "untestable element ids must be recorded"
