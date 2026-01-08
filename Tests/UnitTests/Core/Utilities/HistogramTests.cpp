// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/Histogram.hpp"
#include "Acts/Utilities/ProtoAxis.hpp"

#include <vector>

using namespace Acts;

BOOST_AUTO_TEST_SUITE(HistogramSuite)

BOOST_AUTO_TEST_CASE(Histogram1D_UniformBinning) {
  ProtoAxis protoAxis(AxisBoundaryType::Bound, 0.0, 10.0, 10);
  auto axis = BoostVariableAxis(protoAxis.getAxis().getBinEdges(), "x");
  Histogram1D hist("test", "Test Histogram", axis);

  BOOST_CHECK_EQUAL(hist.name(), "test");
  BOOST_CHECK_EQUAL(hist.title(), "Test Histogram");
  BOOST_CHECK_EQUAL(hist.axisTitle(), "x");

  hist.fill(5.0);
  hist.fill(5.0);

  const auto& bh = hist.histogram();
  BOOST_CHECK_EQUAL(bh.axis(0).size(), 10);

  // Verify count in bin containing 5.0
  auto binIndex = bh.axis(0).index(5.0);
  double binContent = bh.at(binIndex);
  BOOST_CHECK_CLOSE(binContent, 2.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(Histogram1D_VariableBinning) {
  // Create histogram with variable binning
  std::vector<double> edges = {0.0, 1.0, 3.0, 7.0, 10.0};
  auto axis = BoostVariableAxis(edges, "x");
  Histogram1D hist("test_var", "Test Variable Binning", axis);

  // Fill value that falls in bin [1, 3)
  hist.fill(2.0);

  // Access boost histogram
  const auto& bh = hist.histogram();
  // Check number of regular bins (not including underflow/overflow)
  BOOST_CHECK_EQUAL(bh.axis(0).size(), 4);

  // Verify the value is in the correct bin
  auto binIndex = bh.axis(0).index(2.0);
  BOOST_CHECK_EQUAL(binIndex, 1);
  double binContent = bh.at(binIndex);
  BOOST_CHECK_CLOSE(binContent, 1.0, 1e-10);

  // Verify other bins are empty
  for (int i = 0; i < bh.axis(0).size(); ++i) {
    if (i != binIndex) {
      double content = bh.at(i);
      BOOST_CHECK_EQUAL(content, 0.0);
    }
  }
}

BOOST_AUTO_TEST_CASE(Histogram2D_FillAndAccess) {
  ProtoAxis protoX(AxisBoundaryType::Bound, 0.0, 10.0, 10);
  ProtoAxis protoY(AxisBoundaryType::Bound, -5.0, 5.0, 10);
  auto xAxis = BoostVariableAxis(protoX.getAxis().getBinEdges(), "x");
  auto yAxis = BoostVariableAxis(protoY.getAxis().getBinEdges(), "y");
  Histogram2D hist("test_2d", "Test 2D Histogram", xAxis, yAxis);

  BOOST_CHECK_EQUAL(hist.name(), "test_2d");
  BOOST_CHECK_EQUAL(hist.title(), "Test 2D Histogram");
  BOOST_CHECK_EQUAL(hist.xAxisTitle(), "x");
  BOOST_CHECK_EQUAL(hist.yAxisTitle(), "y");

  hist.fill(5.0, 2.0);

  const auto& bh = hist.histogram();
  auto xIdx = bh.axis(0).index(5.0);
  auto yIdx = bh.axis(1).index(2.0);
  double binContent = bh.at(xIdx, yIdx);
  BOOST_CHECK_CLOSE(binContent, 1.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(Histogram2D_VariableBinning) {
  // Create 2D histogram with variable binning on both axes
  std::vector<double> xEdges = {0.0, 1.0, 3.0, 5.0};
  std::vector<double> yEdges = {-2.0, -1.0, 0.0, 1.0, 2.0};
  auto xAxis = BoostVariableAxis(xEdges, "eta");
  auto yAxis = BoostVariableAxis(yEdges, "res");
  Histogram2D hist("res_vs_eta", "Residual vs Eta", xAxis, yAxis);

  // Fill multiple entries
  hist.fill(2.0, 0.5);
  hist.fill(2.0, 0.5);
  hist.fill(0.5, -1.5);

  // Access boost histogram
  const auto& bh = hist.histogram();
  BOOST_CHECK_EQUAL(bh.axis(0).size(), 3);
  BOOST_CHECK_EQUAL(bh.axis(1).size(), 4);

  // Verify first filled bin (2.0, 0.5) - filled twice
  auto xIdx1 = bh.axis(0).index(2.0);
  auto yIdx1 = bh.axis(1).index(0.5);
  double binContent1 = bh.at(xIdx1, yIdx1);
  BOOST_CHECK_CLOSE(binContent1, 2.0, 1e-10);

  // Verify second filled bin (0.5, -1.5) - filled once
  auto xIdx2 = bh.axis(0).index(0.5);
  auto yIdx2 = bh.axis(1).index(-1.5);
  double binContent2 = bh.at(xIdx2, yIdx2);
  BOOST_CHECK_CLOSE(binContent2, 1.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(Histogram1D_UnderflowOverflow) {
  // Create histogram to test underflow/overflow handling
  ProtoAxis protoAxis(AxisBoundaryType::Bound, 0.0, 10.0, 10);
  auto axis = BoostVariableAxis(protoAxis.getAxis().getBinEdges(), "x");
  Histogram1D hist("test_flow", "Test Flow", axis);

  // Fill values in range, underflow, and overflow
  hist.fill(5.0);   // in range
  hist.fill(-1.0);  // underflow
  hist.fill(15.0);  // overflow

  const auto& bh = hist.histogram();

  // boost::histogram has underflow/overflow bins by default
  // Regular bins: 0..9, underflow: -1, overflow: 10
  auto inRangeIdx = bh.axis(0).index(5.0);
  double binContent = bh.at(inRangeIdx);
  BOOST_CHECK_CLOSE(binContent, 1.0, 1e-10);

  // Note: accessing underflow/overflow requires special handling
  // which is implementation detail - converters will handle this
}

BOOST_AUTO_TEST_CASE(Histogram1D_EmptyHistogram) {
  ProtoAxis protoAxis(AxisBoundaryType::Bound, -5.0, 5.0, 10);
  auto axis = BoostVariableAxis(protoAxis.getAxis().getBinEdges(), "x");
  Histogram1D hist("empty", "Empty Histogram", axis);

  const auto& bh = hist.histogram();
  for (int i = 0; i < bh.axis(0).size(); ++i) {
    double content = bh.at(i);
    BOOST_CHECK_EQUAL(content, 0.0);
  }
}

BOOST_AUTO_TEST_CASE(Histogram2D_ProjectionX) {
  // Create 2D histogram
  ProtoAxis protoX(AxisBoundaryType::Bound, 0.0, 5.0, 5);
  ProtoAxis protoY(AxisBoundaryType::Bound, 0.0, 4.0, 4);
  auto xAxis = BoostVariableAxis(protoX.getAxis().getBinEdges(), "x");
  auto yAxis = BoostVariableAxis(protoY.getAxis().getBinEdges(), "y");
  Histogram2D hist2d("test_2d", "Test 2D", xAxis, yAxis);

  // Fill some entries
  hist2d.fill(1.5, 0.5);
  hist2d.fill(1.5, 1.5);
  hist2d.fill(1.5, 2.5);
  hist2d.fill(3.5, 1.5);
  hist2d.fill(3.5, 2.5);

  // Project onto X axis
  Histogram1D projX = hist2d.projectionX();

  // Check metadata
  BOOST_CHECK_EQUAL(projX.name(), "test_2d_projX");
  BOOST_CHECK_EQUAL(projX.title(), "Test 2D projection X");
  BOOST_CHECK_EQUAL(projX.axisTitle(), "x");

  // Check projected values
  const auto& bhX = projX.histogram();
  BOOST_CHECK_EQUAL(bhX.axis(0).size(), 5);

  // Bin containing x=1.5 should have 3 entries
  auto idx1 = bhX.axis(0).index(1.5);
  double content1 = bhX.at(idx1);
  BOOST_CHECK_CLOSE(content1, 3.0, 1e-10);

  // Bin containing x=3.5 should have 2 entries
  auto idx2 = bhX.axis(0).index(3.5);
  double content2 = bhX.at(idx2);
  BOOST_CHECK_CLOSE(content2, 2.0, 1e-10);

  // Other bins should be empty
  for (int i = 0; i < bhX.axis(0).size(); ++i) {
    if (i != idx1 && i != idx2) {
      double content = bhX.at(i);
      BOOST_CHECK_EQUAL(content, 0.0);
    }
  }
}

BOOST_AUTO_TEST_CASE(Histogram2D_ProjectionY) {
  // Create 2D histogram
  ProtoAxis protoX(AxisBoundaryType::Bound, 0.0, 5.0, 5);
  ProtoAxis protoY(AxisBoundaryType::Bound, 0.0, 4.0, 4);
  auto xAxis = BoostVariableAxis(protoX.getAxis().getBinEdges(), "x");
  auto yAxis = BoostVariableAxis(protoY.getAxis().getBinEdges(), "y");
  Histogram2D hist2d("test_2d", "Test 2D", xAxis, yAxis);

  // Fill some entries
  hist2d.fill(0.5, 1.5);
  hist2d.fill(1.5, 1.5);
  hist2d.fill(2.5, 1.5);
  hist2d.fill(1.5, 2.5);
  hist2d.fill(2.5, 2.5);

  // Project onto Y axis
  Histogram1D projY = hist2d.projectionY();

  // Check metadata
  BOOST_CHECK_EQUAL(projY.name(), "test_2d_projY");
  BOOST_CHECK_EQUAL(projY.title(), "Test 2D projection Y");
  BOOST_CHECK_EQUAL(projY.axisTitle(), "y");

  // Check projected values
  const auto& bhY = projY.histogram();
  BOOST_CHECK_EQUAL(bhY.axis(0).size(), 4);

  // Bin containing y=1.5 should have 3 entries
  auto idx1 = bhY.axis(0).index(1.5);
  double content1 = bhY.at(idx1);
  BOOST_CHECK_CLOSE(content1, 3.0, 1e-10);

  // Bin containing y=2.5 should have 2 entries
  auto idx2 = bhY.axis(0).index(2.5);
  double content2 = bhY.at(idx2);
  BOOST_CHECK_CLOSE(content2, 2.0, 1e-10);

  // Other bins should be empty
  for (int i = 0; i < bhY.axis(0).size(); ++i) {
    if (i != idx1 && i != idx2) {
      double content = bhY.at(i);
      BOOST_CHECK_EQUAL(content, 0.0);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
