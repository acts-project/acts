# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pytest

import acts

bh = pytest.importorskip("boost_histogram")
mplhep = pytest.importorskip("mplhep")


def test_histogram1_to_boost_histogram():
    bh_h = bh.Histogram(acts._demo_histogram1())
    assert bh_h.ndim == 1
    assert bh_h.sum() > 0
    assert bh_h.axes[0].metadata == "x [a.u.]"


def test_profile1_to_boost_histogram():
    bh_h = bh.Histogram(acts._demo_profile1())
    assert bh_h.ndim == 1
    assert bh_h.view()["count"].sum() > 0


def test_efficiency1_to_boost_histogram():
    h = acts._demo_efficiency1()
    bh_acc = bh.Histogram(h.accepted)
    bh_tot = bh.Histogram(h.total)
    assert bh_acc.sum() > 0
    assert bh_tot.sum() >= bh_acc.sum()


def test_plot_histogram():
    fig, ax = plt.subplots()
    result = acts._demo_histogram1().plot(ax=ax)
    assert result is not None
    assert ax.get_xlabel() == "x [a.u.]"
    assert ax.get_title() == "Demo Histogram"
    plt.close(fig)


def test_plot_histogram_no_ax():
    result = acts._demo_histogram1().plot()
    assert result is not None
    plt.close("all")


def test_plot_profile():
    fig, ax = plt.subplots()
    result = acts._demo_profile1().plot(ax=ax)
    assert result is not None
    assert ax.get_xlabel() == "x [a.u.]"
    assert ax.get_ylabel() == "y [a.u.]"
    assert ax.get_title() == "Demo Profile"
    plt.close(fig)


def test_plot_efficiency():
    fig, ax = plt.subplots()
    acts._demo_efficiency1().plot(ax=ax)
    assert ax.get_title() == "Demo Efficiency"
    plt.close(fig)
