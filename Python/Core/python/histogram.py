# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.


def _patch_histogram_types(m):
    """Patch ACTS histogram types with the boost-histogram protocol.

    Adds _to_boost_histogram_() to BoostHistogram, BoostProfileHistogram,
    Histogram1/2/3, ProfileHistogram1, and Efficiency1/2 so that
    bh.Histogram(acts_obj) works.

    boost_histogram is imported lazily on first conversion; an ImportError
    is printed and re-raised if the package is not installed.
    """

    def _import_bh():
        try:
            import boost_histogram as bh

            return bh
        except ImportError as e:
            print(
                "boost_histogram is not installed; cannot convert ACTS histogram to bh.Histogram."
            )
            raise

    def _bh_axes(bh, self):
        return [
            bh.axis.Variable(list(self.axis(i).edges), metadata=self.axis(i).label)
            for i in range(self.rank)
        ]

    # BoostHistogram -> bh.Histogram with Double storage
    def _boost_hist_to_bh(self):
        bh = _import_bh()
        h = bh.Histogram(*_bh_axes(bh, self), storage=bh.storage.Double())
        h.view(flow=False)[:] = self.values()
        return h

    m.BoostHistogram._to_boost_histogram_ = _boost_hist_to_bh

    # BoostProfileHistogram -> bh.Histogram with Mean storage
    def _profile_to_bh(self):
        bh = _import_bh()
        h = bh.Histogram(*_bh_axes(bh, self), storage=bh.storage.Mean())
        view = h.view()
        view["value"] = self.means()
        view["variance"] = self.variances()
        return h

    m.BoostProfileHistogram._to_boost_histogram_ = _profile_to_bh

    # Histogram<Dim>: delegate to .histogram() which is a BoostHistogram
    for cls in (m.Histogram1, m.Histogram2, m.Histogram3):
        cls._to_boost_histogram_ = lambda self: self.histogram()._to_boost_histogram_()

    # ProfileHistogram<1>: delegate to .histogram() which is a BoostProfileHistogram
    m.ProfileHistogram1._to_boost_histogram_ = (
        lambda self: self.histogram()._to_boost_histogram_()
    )
