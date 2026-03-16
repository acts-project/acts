# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import numpy as np
import pytest

import acts

bh = pytest.importorskip("boost_histogram")


# ---------------------------------------------------------------------------
# Axis construction
# ---------------------------------------------------------------------------


def test_axis_regular():
    ax = acts.Axis(10, 0.0, 1.0, "x")
    assert ax.size == 10
    assert ax.label == "x"
    edges = ax.edges
    assert len(edges) == 11
    assert edges[0] == pytest.approx(0.0)
    assert edges[-1] == pytest.approx(1.0)


def test_axis_regular_no_label():
    ax = acts.Axis(5, -1.0, 1.0)
    assert ax.size == 5
    assert ax.label == ""


def test_axis_variable():
    edges_in = [0.0, 1.0, 3.0, 6.0]
    ax = acts.Axis(edges_in, "eta")
    assert ax.size == 3
    assert ax.label == "eta"
    edges_out = ax.edges
    assert len(edges_out) == 4
    for a, b in zip(edges_out, edges_in):
        assert a == pytest.approx(b)


# ---------------------------------------------------------------------------
# Histogram construction and fill
# ---------------------------------------------------------------------------


def test_histogram1_construction():
    ax = acts.Axis(10, 0.0, 1.0, "x")
    h = acts.Histogram1("h1", "title 1", [ax])
    assert h.name == "h1"
    assert h.title == "title 1"
    assert h.rank == 1


def test_histogram2_construction():
    ax = acts.Axis(10, 0.0, 1.0, "x")
    ay = acts.Axis(5, 0.0, 5.0, "y")
    h = acts.Histogram2("h2", "title 2", [ax, ay])
    assert h.rank == 2


def test_histogram3_construction():
    ax = acts.Axis(4, 0.0, 1.0, "x")
    ay = acts.Axis(4, 0.0, 1.0, "y")
    az = acts.Axis(4, 0.0, 1.0, "z")
    h = acts.Histogram3("h3", "title 3", [ax, ay, az])
    assert h.rank == 3


def test_histogram1_fill_and_values():
    ax = acts.Axis(4, 0.0, 4.0, "x")
    h = acts.Histogram1("h", "t", [ax])
    h.fill([0.5])
    h.fill([0.5])
    h.fill([1.5])
    vals = h.histogram().values()
    assert vals.shape == (4,)
    assert vals[0] == pytest.approx(2.0)
    assert vals[1] == pytest.approx(1.0)
    assert vals[2] == pytest.approx(0.0)
    assert vals[3] == pytest.approx(0.0)


def test_histogram2_fill_and_values():
    ax = acts.Axis(2, 0.0, 2.0, "x")
    ay = acts.Axis(2, 0.0, 2.0, "y")
    h = acts.Histogram2("h2", "t", [ax, ay])
    h.fill([0.5, 0.5])
    h.fill([0.5, 1.5])
    h.fill([1.5, 0.5])
    vals = h.histogram().values()
    assert vals.shape == (2, 2)
    assert vals[0, 0] == pytest.approx(1.0)
    assert vals[0, 1] == pytest.approx(1.0)
    assert vals[1, 0] == pytest.approx(1.0)
    assert vals[1, 1] == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# BoostHistogram axis access
# ---------------------------------------------------------------------------


def test_boost_histogram_axis():
    ax = acts.Axis(10, 0.0, 1.0, "myaxis")
    h = acts.Histogram1("h", "t", [ax])
    bh_raw = h.histogram()
    assert bh_raw.rank == 1
    ax0 = bh_raw.axis(0)
    assert ax0.size == 10
    assert ax0.label == "myaxis"


# ---------------------------------------------------------------------------
# Conversion to boost_histogram via _to_boost_histogram_() protocol
# ---------------------------------------------------------------------------


def test_histogram1_to_bh():
    ax = acts.Axis(4, 0.0, 4.0, "x")
    h = acts.Histogram1("h", "t", [ax])
    h.fill([0.5])
    h.fill([1.5])
    bh_h = bh.Histogram(h)
    assert bh_h.ndim == 1
    assert bh_h.axes[0].extent == 6  # 4 bins + 2 flow
    assert bh_h.sum() == pytest.approx(2.0)


def test_histogram2_to_bh():
    ax = acts.Axis(3, 0.0, 3.0, "x")
    ay = acts.Axis(2, 0.0, 2.0, "y")
    h = acts.Histogram2("h2", "t", [ax, ay])
    h.fill([0.5, 0.5])
    h.fill([0.5, 1.5])
    bh_h = bh.Histogram(h)
    assert bh_h.ndim == 2
    assert bh_h.sum() == pytest.approx(2.0)


def test_histogram1_to_bh_axes_metadata():
    ax = acts.Axis(5, 0.0, 5.0, "mylabel")
    h = acts.Histogram1("h", "t", [ax])
    bh_h = bh.Histogram(h)
    assert bh_h.axes[0].metadata == "mylabel"


# ---------------------------------------------------------------------------
# ProfileHistogram1
# ---------------------------------------------------------------------------


def test_profile_histogram1_fill_and_means():
    ax = acts.Axis(10, 0.0, 10.0, "x")
    ph = acts.ProfileHistogram1("p", "t", [ax], "y value")
    assert ph.name == "p"
    assert ph.sampleAxisTitle == "y value"
    assert ph.rank == 1
    # Fill bin 2 (x in [2,3)) with sample values 10, 20, 30 -> mean=20
    ph.fill([2.5], 10.0)
    ph.fill([2.5], 20.0)
    ph.fill([2.5], 30.0)
    means = ph.histogram().means()
    assert means.shape == (10,)
    assert means[2] == pytest.approx(20.0)


def test_profile_histogram1_variances():
    ax = acts.Axis(10, 0.0, 10.0, "x")
    ph = acts.ProfileHistogram1("p", "t", [ax], "y")
    ph.fill([0.5], 1.0)
    ph.fill([0.5], 3.0)
    # sample variance of [1, 3] = ((1-2)^2 + (3-2)^2) / (2-1) = 2
    variances = ph.histogram().variances()
    assert variances.shape == (10,)
    assert variances[0] == pytest.approx(2.0)


def test_profile_histogram1_to_bh():
    ax = acts.Axis(5, 0.0, 5.0, "x")
    ph = acts.ProfileHistogram1("p", "t", [ax], "y")
    ph.fill([0.5], 10.0)
    ph.fill([0.5], 20.0)
    bh_h = bh.Histogram(ph)
    assert bh_h.ndim == 1
    view = bh_h.view()
    assert view["value"][0] == pytest.approx(15.0)


# ---------------------------------------------------------------------------
# Efficiency1 and Efficiency2
# ---------------------------------------------------------------------------


def test_efficiency1_fill():
    ax = acts.Axis(4, 0.0, 4.0, "x")
    eff = acts.Efficiency1("eff", "t", [ax])
    assert eff.rank == 1
    eff.fill([0.5], True)
    eff.fill([0.5], True)
    eff.fill([0.5], False)
    eff.fill([1.5], False)
    accepted = eff.accepted().values()
    total = eff.total().values()
    assert accepted.shape == (4,)
    assert total.shape == (4,)
    assert accepted[0] == pytest.approx(2.0)
    assert total[0] == pytest.approx(3.0)
    assert accepted[1] == pytest.approx(0.0)
    assert total[1] == pytest.approx(1.0)


def test_efficiency2_fill():
    ax = acts.Axis(2, 0.0, 2.0, "x")
    ay = acts.Axis(2, 0.0, 2.0, "y")
    eff = acts.Efficiency2("eff2", "t", [ax, ay])
    assert eff.rank == 2
    eff.fill([0.5, 0.5], True)
    eff.fill([0.5, 0.5], False)
    eff.fill([1.5, 1.5], True)
    total = eff.total().values()
    accepted = eff.accepted().values()
    assert total.shape == (2, 2)
    assert total[0, 0] == pytest.approx(2.0)
    assert total[1, 1] == pytest.approx(1.0)
    assert accepted[0, 0] == pytest.approx(1.0)
    assert accepted[1, 1] == pytest.approx(1.0)
