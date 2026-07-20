"""Unified trace-event taxonomy (spec §3.8; Plan-07 reconciliation)."""
from dataclasses import fields as dc_fields

import pytest

import quiverlab
from quiverlab.trace.events import (
    Dispatch, ReductionStep, AmbiguityEvent, ResolutionTerm,
    DifferentialEvent, LiftStep, RankStep,
)


def _names(cls):
    return {f.name for f in dc_fields(cls)}


def test_reexports_are_the_shipped_classes_not_copies():
    import quiverlab.groebner.events as ge
    import quiverlab.resolutions_cs.trace as ct
    assert Dispatch is ge.Dispatch
    assert ReductionStep is ge.ReductionStep
    assert AmbiguityEvent is ct.AmbiguityEvent
    assert ResolutionTerm is ct.ResolutionTerm
    assert DifferentialEvent is ct.DifferentialEvent
    assert LiftStep is ct.LiftStep


def test_rankstep_is_new_and_shaped():
    assert _names(RankStep) == {
        "degree", "side", "nrows", "ncols", "rank", "field", "matrix", "elided", "note"}
    r = RankStep(degree=1, side="cochain", nrows=2, ncols=2, rank=1, field="CC",
                 matrix=[["0", "0"], ["2", "0"]])
    assert r.rank == 1 and r.elided is False and r.note == "" and r.matrix[1][0] == "2"


@pytest.mark.verbose_default
def test_verbose_flag_defaults_true():
    assert quiverlab.verbose is True
