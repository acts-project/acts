// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/Histogram.hpp"

#include <array>
#include <vector>

using namespace Acts;
using namespace Acts::Experimental;

BOOST_AUTO_TEST_SUITE(HistogramSuite)

BOOST_AUTO_TEST_CASE(Histogram1D_UniformBinning) {
  auto axis = AxisVariant(BoostRegularAxis(10, 0.0, 10.0, "x"));
  Histogram1 hist("test", "Test Histogram", {axis});

  BOOST_CHECK_EQUAL(hist.name(), "test");
  BOOST_CHECK_EQUAL(hist.title(), "Test Histogram");

  hist.fill({5.0});
  hist.fill({5.0});

  const auto& bh = hist.histogram();
  BOOST_CHECK_EQUAL(bh.axis(0).size(), 10);

  // Verify count in bin containing 5.0
  auto binIndex = bh.axis(0).index(5.0);
  double binContent = bh.at(binIndex).value();
  BOOST_CHECK_CLOSE(binContent, 2.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(Histogram1D_VariableBinning) {
  // Create histogram with variable binning
  std::vector<double> edges = {0.0, 1.0, 3.0, 7.0, 10.0};
  auto axis = AxisVariant(BoostVariableAxis(edges, "x"));
  Histogram1 hist("test_var", "Test Variable Binning", {axis});

  // Fill value that falls in bin [1, 3)
  hist.fill({2.0});

  // Access boost histogram
  const auto& bh = hist.histogram();
  // Check number of regular bins (not including underflow/overflow)
  BOOST_CHECK_EQUAL(bh.axis(0).size(), 4);

  // Verify the value is in the correct bin
  auto binIndex = bh.axis(0).index(2.0);
  BOOST_CHECK_EQUAL(binIndex, 1);
  double binContent = bh.at(binIndex).value();
  BOOST_CHECK_CLOSE(binContent, 1.0, 1e-10);

  // Verify other bins are empty
  for (int i = 0; i < bh.axis(0).size(); ++i) {
    if (i != binIndex) {
      double content = bh.at(i).value();
      BOOST_CHECK_EQUAL(content, 0.0);
    }
  }
}

BOOST_AUTO_TEST_CASE(Histogram2D_FillAndAccess) {
  auto xAxis = AxisVariant(BoostRegularAxis(10, 0.0, 10.0, "x"));
  auto yAxis = AxisVariant(BoostRegularAxis(10, -5.0, 5.0, "y"));
  Histogram2 hist("test_2d", "Test 2D Histogram", {xAxis, yAxis});

  BOOST_CHECK_EQUAL(hist.name(), "test_2d");
  BOOST_CHECK_EQUAL(hist.title(), "Test 2D Histogram");

  hist.fill({5.0, 2.0});

  const auto& bh = hist.histogram();
  auto xIdx = bh.axis(0).index(5.0);
  auto yIdx = bh.axis(1).index(2.0);
  double binContent = bh.at(xIdx, yIdx).value();
  BOOST_CHECK_CLOSE(binContent, 1.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(Histogram2D_VariableBinning) {
  // Create 2D histogram with variable binning on both axes
  std::vector<double> xEdges = {0.0, 1.0, 3.0, 5.0};
  std::vector<double> yEdges = {-2.0, -1.0, 0.0, 1.0, 2.0};
  auto xAxis = AxisVariant(BoostVariableAxis(xEdges, "eta"));
  auto yAxis = AxisVariant(BoostVariableAxis(yEdges, "res"));
  Histogram2 hist("res_vs_eta", "Residual vs Eta", {xAxis, yAxis});

  // Fill multiple entries
  hist.fill({2.0, 0.5});
  hist.fill({2.0, 0.5});
  hist.fill({0.5, -1.5});

  // Access boost histogram
  const auto& bh = hist.histogram();
  BOOST_CHECK_EQUAL(bh.axis(0).size(), 3);
  BOOST_CHECK_EQUAL(bh.axis(1).size(), 4);

  // Verify first filled bin (2.0, 0.5) - filled twice
  auto xIdx1 = bh.axis(0).index(2.0);
  auto yIdx1 = bh.axis(1).index(0.5);
  double binContent1 = bh.at(xIdx1, yIdx1).value();
  BOOST_CHECK_CLOSE(binContent1, 2.0, 1e-10);

  // Verify second filled bin (0.5, -1.5) - filled once
  auto xIdx2 = bh.axis(0).index(0.5);
  auto yIdx2 = bh.axis(1).index(-1.5);
  double binContent2 = bh.at(xIdx2, yIdx2).value();
  BOOST_CHECK_CLOSE(binContent2, 1.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(Histogram1D_UnderflowOverflow) {
  // Create histogram to test underflow/overflow handling
  auto axis = AxisVariant(BoostRegularAxis(10, 0.0, 10.0, "x"));
  Histogram1 hist("test_flow", "Test Flow", {axis});

  // Fill values in range, underflow, and overflow
  hist.fill({5.0});   // in range
  hist.fill({-1.0});  // underflow
  hist.fill({15.0});  // overflow

  const auto& bh = hist.histogram();

  // boost::histogram has underflow/overflow bins by default
  // Regular bins: 0..9, underflow: -1, overflow: 10
  auto inRangeIdx = bh.axis(0).index(5.0);
  double binContent = bh.at(inRangeIdx).value();
  BOOST_CHECK_CLOSE(binContent, 1.0, 1e-10);

  // Note: accessing underflow/overflow requires special handling
  // which is implementation detail - converters will handle this
}

BOOST_AUTO_TEST_CASE(Histogram1D_EmptyHistogram) {
  auto axis = AxisVariant(BoostRegularAxis(10, -5.0, 5.0, "x"));
  Histogram1 hist("empty", "Empty Histogram", {axis});

  const auto& bh = hist.histogram();
  for (int i = 0; i < bh.axis(0).size(); ++i) {
    double content = bh.at(i).value();
    BOOST_CHECK_EQUAL(content, 0.0);
  }
}

BOOST_AUTO_TEST_CASE(Histogram_SetAndGetBinContent) {
  auto axis = AxisVariant(BoostRegularAxis(4, 0.0, 4.0, "x"));
  Histogram1 hist("set_get", "Set/Get", {axis});

  hist.setBinContent({2}, 17.5);
  BOOST_CHECK_CLOSE(hist.binContent({2}), 17.5, 1e-10);

  // Setting must overwrite, not accumulate
  hist.setBinContent({2}, 3.0);
  BOOST_CHECK_CLOSE(hist.binContent({2}), 3.0, 1e-10);

  // A fill and an explicit set must be visible through the same accessor
  hist.fill({0.5});
  BOOST_CHECK_CLOSE(hist.binContent({0}), 1.0, 1e-10);

  // Untouched bins stay empty
  BOOST_CHECK_EQUAL(hist.binContent({1}), 0.0);
  BOOST_CHECK_EQUAL(hist.binContent({3}), 0.0);
}

BOOST_AUTO_TEST_CASE(Histogram2D_SetAndGetBinContent) {
  auto xAxis = AxisVariant(BoostRegularAxis(3, 0.0, 3.0, "x"));
  auto yAxis = AxisVariant(BoostRegularAxis(2, 0.0, 2.0, "y"));
  Histogram2 hist("set_get_2d", "Set/Get 2D", {xAxis, yAxis});

  hist.setBinContent({2, 1}, 7.0);
  BOOST_CHECK_CLOSE(hist.binContent({2, 1}), 7.0, 1e-10);
  BOOST_CHECK_EQUAL(hist.binContent({0, 0}), 0.0);
}

// Regression test: projectionX/Y used to build the projected axis but never
// copy the bin contents, so they returned an empty histogram.
BOOST_AUTO_TEST_CASE(Histogram2D_ProjectionX_CopiesContents) {
  // Asymmetric binning so an axis mix-up cannot pass unnoticed
  std::vector<double> xEdges = {0.0, 1.0, 3.0, 5.0};
  std::vector<double> yEdges = {-2.0, -1.0, 0.0, 1.0, 2.0};
  auto xAxis = AxisVariant(BoostVariableAxis(xEdges, "eta"));
  auto yAxis = AxisVariant(BoostVariableAxis(yEdges, "res"));
  Histogram2 hist("res_vs_eta", "Residual vs Eta", {xAxis, yAxis});

  // Known pattern: content[xBin][yBin]
  const std::array<std::array<double, 4>, 3> pattern = {
      {{1, 2, 0, 3}, {0, 4, 5, 0}, {6, 0, 0, 7}}};
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 4; ++j) {
      hist.setBinContent({i, j}, pattern[i][j]);
    }
  }

  const Histogram1 projX = projectionX(hist);
  BOOST_CHECK_EQUAL(projX.histogram().axis(0).size(), 3);

  double total = 0;
  for (int i = 0; i < 3; ++i) {
    // Projection onto X sums over the Y bins of each X column
    double expected = 0;
    for (int j = 0; j < 4; ++j) {
      expected += pattern[i][j];
    }
    BOOST_CHECK_CLOSE(projX.binContent({i}), expected, 1e-10);
    total += projX.binContent({i});
  }

  // The bug produced an all-zero histogram, so guard the integral explicitly
  BOOST_CHECK_CLOSE(total, 28.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(Histogram2D_ProjectionY_CopiesContents) {
  std::vector<double> xEdges = {0.0, 1.0, 3.0, 5.0};
  std::vector<double> yEdges = {-2.0, -1.0, 0.0, 1.0, 2.0};
  auto xAxis = AxisVariant(BoostVariableAxis(xEdges, "eta"));
  auto yAxis = AxisVariant(BoostVariableAxis(yEdges, "res"));
  Histogram2 hist("res_vs_eta", "Residual vs Eta", {xAxis, yAxis});

  const std::array<std::array<double, 4>, 3> pattern = {
      {{1, 2, 0, 3}, {0, 4, 5, 0}, {6, 0, 0, 7}}};
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 4; ++j) {
      hist.setBinContent({i, j}, pattern[i][j]);
    }
  }

  const Histogram1 projY = projectionY(hist);
  BOOST_CHECK_EQUAL(projY.histogram().axis(0).size(), 4);

  double total = 0;
  for (int j = 0; j < 4; ++j) {
    // Projection onto Y sums over the X bins of each Y row
    double expected = 0;
    for (int i = 0; i < 3; ++i) {
      expected += pattern[i][j];
    }
    BOOST_CHECK_CLOSE(projY.binContent({j}), expected, 1e-10);
    total += projY.binContent({j});
  }

  BOOST_CHECK_CLOSE(total, 28.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(Histogram2D_Projection_PreservesAxis) {
  std::vector<double> xEdges = {0.0, 1.0, 3.0, 5.0};
  std::vector<double> yEdges = {-2.0, -1.0, 0.0, 1.0, 2.0};
  auto xAxis = AxisVariant(BoostVariableAxis(xEdges, "eta"));
  auto yAxis = AxisVariant(BoostVariableAxis(yEdges, "res"));
  Histogram2 hist("res_vs_eta", "Residual vs Eta", {xAxis, yAxis});

  const Histogram1 projX = projectionX(hist);
  BOOST_CHECK_EQUAL(projX.name(), "res_vs_eta_projX");
  BOOST_CHECK_EQUAL(projX.title(), "Residual vs Eta projection X");
  BOOST_CHECK_EQUAL(projX.histogram().axis(0).metadata(), "eta");
  BOOST_CHECK(extractBinEdges(projX.histogram().axis(0)) == xEdges);

  const Histogram1 projY = projectionY(hist);
  BOOST_CHECK_EQUAL(projY.name(), "res_vs_eta_projY");
  BOOST_CHECK_EQUAL(projY.title(), "Residual vs Eta projection Y");
  BOOST_CHECK_EQUAL(projY.histogram().axis(0).metadata(), "res");
  BOOST_CHECK(extractBinEdges(projY.histogram().axis(0)) == yEdges);
}

BOOST_AUTO_TEST_CASE(Histogram2D_Projection_IncludesFlowBins) {
  // boost::histogram::algorithm::project sums over the flow bins of the
  // reduced axis, unlike ROOT's TH2::ProjectionX/Y. Pin that behaviour down so
  // a future change to the projection helpers cannot alter it silently.
  auto xAxis = AxisVariant(BoostRegularAxis(2, 0.0, 2.0, "x"));
  auto yAxis = AxisVariant(BoostRegularAxis(2, 0.0, 2.0, "y"));
  Histogram2 hist("flow", "Flow", {xAxis, yAxis});

  hist.fill({0.5, 0.5});   // both in range
  hist.fill({0.5, 99.0});  // Y overflow, X bin 0
  hist.fill({-5.0, 0.5});  // X underflow, Y bin 0

  const Histogram1 projX = projectionX(hist);
  // X bin 0 picks up the in-range entry *and* the Y-overflow entry
  BOOST_CHECK_CLOSE(projX.binContent({0}), 2.0, 1e-10);
  BOOST_CHECK_EQUAL(projX.binContent({1}), 0.0);

  const Histogram1 projY = projectionY(hist);
  // Y bin 0 picks up the in-range entry *and* the X-underflow entry
  BOOST_CHECK_CLOSE(projY.binContent({0}), 2.0, 1e-10);
  BOOST_CHECK_EQUAL(projY.binContent({1}), 0.0);
}

BOOST_AUTO_TEST_CASE(SliceLastAxis_2D) {
  // Asymmetric binning so swapping the axes cannot pass
  std::vector<double> xEdges = {0.0, 1.0, 3.0, 5.0};
  std::vector<double> yEdges = {-2.0, -1.0, 0.0, 1.0, 2.0};
  auto xAxis = AxisVariant(BoostVariableAxis(xEdges, "eta"));
  auto yAxis = AxisVariant(BoostVariableAxis(yEdges, "res"));
  Histogram2 hist("res_vs_eta", "Residual vs Eta", {xAxis, yAxis});

  const std::array<std::array<double, 4>, 3> pattern = {
      {{1, 2, 0, 3}, {0, 4, 5, 0}, {6, 0, 0, 7}}};
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 4; ++j) {
      hist.setBinContent({i, j}, pattern[i][j]);
    }
  }

  for (int i = 0; i < 3; ++i) {
    const Histogram1 slice = sliceLastAxis(hist, i);

    // The slice spans the last axis and inherits its binning and metadata
    BOOST_CHECK_EQUAL(slice.histogram().axis(0).size(), 4);
    BOOST_CHECK_EQUAL(slice.histogram().axis(0).metadata(), "res");
    BOOST_CHECK(extractBinEdges(slice.histogram().axis(0)) == yEdges);

    for (int j = 0; j < 4; ++j) {
      BOOST_CHECK_CLOSE(slice.binContent({j}), pattern[i][j], 1e-10);
    }
  }
}

BOOST_AUTO_TEST_CASE(SliceLastAxis_2D_IsNotAProjection) {
  // A slice must see one column only, not sum over the sliced axis (which a
  // naive reduce()+project() would do).
  auto xAxis = AxisVariant(BoostRegularAxis(3, 0.0, 3.0, "x"));
  auto yAxis = AxisVariant(BoostRegularAxis(2, 0.0, 2.0, "y"));
  Histogram2 hist("noproj", "No projection", {xAxis, yAxis});

  hist.setBinContent({0, 0}, 1.0);
  hist.setBinContent({1, 0}, 10.0);
  hist.setBinContent({2, 1}, 100.0);
  // Content outside the x range must not leak in either
  hist.fill({-5.0, 0.5});
  hist.fill({99.0, 0.5});

  const Histogram1 slice0 = sliceLastAxis(hist, 0);
  BOOST_CHECK_CLOSE(slice0.binContent({0}), 1.0, 1e-10);
  BOOST_CHECK_EQUAL(slice0.binContent({1}), 0.0);

  const Histogram1 slice2 = sliceLastAxis(hist, 2);
  BOOST_CHECK_EQUAL(slice2.binContent({0}), 0.0);
  BOOST_CHECK_CLOSE(slice2.binContent({1}), 100.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(SliceLastAxis_3D) {
  auto xAxis = AxisVariant(BoostRegularAxis(2, 0.0, 2.0, "eta"));
  auto yAxis = AxisVariant(BoostRegularAxis(3, 0.0, 3.0, "pt"));
  auto zAxis = AxisVariant(BoostRegularAxis(4, -2.0, 2.0, "res"));
  Histogram3 hist("res_vs_eta_pt", "Residual", {xAxis, yAxis, zAxis});

  // Distinct value per (i, j, k) so any index mix-up shows up
  const auto encode = [](int i, int j, int k) {
    return 100.0 * i + 10.0 * j + k + 1;
  };
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 4; ++k) {
        hist.setBinContent({i, j, k}, encode(i, j, k));
      }
    }
  }

  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 3; ++j) {
      const Histogram1 slice = sliceLastAxis(hist, i, j);
      BOOST_CHECK_EQUAL(slice.histogram().axis(0).size(), 4);
      BOOST_CHECK_EQUAL(slice.histogram().axis(0).metadata(), "res");
      for (int k = 0; k < 4; ++k) {
        BOOST_CHECK_CLOSE(slice.binContent({k}), encode(i, j, k), 1e-10);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(Histogram1D_SetBinWithError) {
  auto axis = AxisVariant(BoostRegularAxis(5, 0.0, 5.0, "eta"));
  Histogram1 hist("mean", "Mean", {axis});

  // Untouched bins read back as zero content and zero error
  BOOST_CHECK_EQUAL(hist.binContent({0}), 0.0);
  BOOST_CHECK_EQUAL(hist.binError({0}), 0.0);

  hist.setBin({2}, -1.25, 0.5);
  BOOST_CHECK_CLOSE(hist.binContent({2}), -1.25, 1e-10);
  BOOST_CHECK_CLOSE(hist.binError({2}), 0.5, 1e-10);

  // Setting overwrites rather than accumulates
  hist.setBin({2}, 3.0, 0.25);
  BOOST_CHECK_CLOSE(hist.binContent({2}), 3.0, 1e-10);
  BOOST_CHECK_CLOSE(hist.binError({2}), 0.25, 1e-10);

  // Neighbours stay untouched
  BOOST_CHECK_EQUAL(hist.binContent({1}), 0.0);
  BOOST_CHECK_EQUAL(hist.binError({3}), 0.0);

  // The axis is carried through for converters
  BOOST_CHECK_EQUAL(hist.histogram().axis(0).size(), 5);
  BOOST_CHECK_EQUAL(hist.histogram().axis(0).metadata(), "eta");
}

BOOST_AUTO_TEST_CASE(Histogram2D_SetBinWithError) {
  std::vector<double> xEdges = {0.0, 1.0, 3.0};
  auto xAxis = AxisVariant(BoostVariableAxis(xEdges, "eta"));
  auto yAxis = AxisVariant(BoostRegularAxis(3, 0.0, 3.0, "pt"));
  Histogram2 hist("width", "Width", {xAxis, yAxis});

  hist.setBin({1, 2}, 0.75, 0.1);
  BOOST_CHECK_CLOSE(hist.binContent({1, 2}), 0.75, 1e-10);
  BOOST_CHECK_CLOSE(hist.binError({1, 2}), 0.1, 1e-10);
  BOOST_CHECK_EQUAL(hist.binContent({0, 0}), 0.0);

  BOOST_CHECK(extractBinEdges(hist.histogram().axis(0)) == xEdges);
  BOOST_CHECK_EQUAL(hist.histogram().axis(1).metadata(), "pt");
}

BOOST_AUTO_TEST_CASE(Histogram_ZeroErrorIsAllowed) {
  auto axis = AxisVariant(BoostRegularAxis(2, 0.0, 2.0, "x"));
  Histogram1 hist("zero", "Zero", {axis});

  hist.setBin({0}, 5.0, 0.0);
  BOOST_CHECK_CLOSE(hist.binContent({0}), 5.0, 1e-10);
  BOOST_CHECK_EQUAL(hist.binError({0}), 0.0);
}

BOOST_AUTO_TEST_CASE(Histogram_Fill_GivesSqrtNError) {
  // A plain fill() should give the usual counting-statistics error, matching
  // what ROOT's TH1::Fill (with Sumw2 enabled) would report.
  auto axis = AxisVariant(BoostRegularAxis(4, 0.0, 4.0, "x"));
  Histogram1 hist("counts", "Counts", {axis});

  for (int i = 0; i < 9; ++i) {
    hist.fill({0.5});
  }

  BOOST_CHECK_CLOSE(hist.binContent({0}), 9.0, 1e-10);
  BOOST_CHECK_CLOSE(hist.binError({0}), 3.0, 1e-10);
  BOOST_CHECK_EQUAL(hist.binContent({1}), 0.0);
  BOOST_CHECK_EQUAL(hist.binError({1}), 0.0);
}

BOOST_AUTO_TEST_CASE(SliceLastAxis_PropagatesErrors) {
  auto xAxis = AxisVariant(BoostRegularAxis(2, 0.0, 2.0, "eta"));
  auto yAxis = AxisVariant(BoostRegularAxis(3, 0.0, 3.0, "res"));
  Histogram2 hist("res_vs_eta", "Residual vs Eta", {xAxis, yAxis});

  hist.setBin({0, 0}, 1.0, 0.1);
  hist.setBin({0, 1}, 2.0, 0.2);
  hist.setBin({0, 2}, 3.0, 0.3);

  const Histogram1 slice = sliceLastAxis(hist, 0);
  BOOST_CHECK_CLOSE(slice.binContent({0}), 1.0, 1e-10);
  BOOST_CHECK_CLOSE(slice.binError({0}), 0.1, 1e-10);
  BOOST_CHECK_CLOSE(slice.binContent({1}), 2.0, 1e-10);
  BOOST_CHECK_CLOSE(slice.binError({1}), 0.2, 1e-10);
  BOOST_CHECK_CLOSE(slice.binContent({2}), 3.0, 1e-10);
  BOOST_CHECK_CLOSE(slice.binError({2}), 0.3, 1e-10);
}

BOOST_AUTO_TEST_SUITE_END()
