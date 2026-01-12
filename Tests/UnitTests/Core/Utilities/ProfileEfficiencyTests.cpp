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

using namespace Acts;
using namespace Acts::Experimental;

BOOST_AUTO_TEST_SUITE(ProfileEfficiencySuite)

BOOST_AUTO_TEST_CASE(ProfileHistogram_BasicFill) {
  ProtoAxis protoAxis(AxisBoundaryType::Bound, 0.0, 10.0, 10);
  auto xAxis = BoostVariableAxis(protoAxis.getAxis().getBinEdges(), "x");
  ProfileHistogram1D profile("test_prof", "Test Profile", {xAxis}, "y value");

  BOOST_CHECK_EQUAL(profile.name(), "test_prof");
  BOOST_CHECK_EQUAL(profile.title(), "Test Profile");
  BOOST_CHECK_EQUAL(profile.sampleAxisTitle(), "y value");

  // Fill x=2.5 with y=10.0, y=20.0, y=30.0
  profile.fill({2.5}, 10.0);
  profile.fill({2.5}, 20.0);
  profile.fill({2.5}, 30.0);

  const auto& bh = profile.histogram();
  auto xIdx = bh.axis(0).index(2.5);

  // Mean should be 20.0
  double mean = bh.at(xIdx).value();
  BOOST_CHECK_CLOSE(mean, 20.0, 1e-6);

  // Count should be 3
  double count = bh.at(xIdx).count();
  BOOST_CHECK_CLOSE(count, 3.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(ProfileHistogram_MultipleBins) {
  ProtoAxis protoAxis(AxisBoundaryType::Bound, -2.5, 2.5, 5);
  auto xAxis = BoostVariableAxis(protoAxis.getAxis().getBinEdges(), "eta");
  ProfileHistogram1D profile("res_vs_eta", "Residual vs Eta", {xAxis},
                             "residual");

  // Fill different eta bins with different mean values
  profile.fill({-2.0}, 1.0);  // bin 0
  profile.fill({-2.0}, 3.0);
  profile.fill({0.0}, 5.0);  // bin 2 (middle)
  profile.fill({0.0}, 7.0);
  profile.fill({2.0}, 9.0);  // bin 4
  profile.fill({2.0}, 11.0);

  const auto& bh = profile.histogram();

  // Check mean in first bin
  auto idx0 = bh.axis(0).index(-2.0);
  BOOST_CHECK_CLOSE(bh.at(idx0).value(), 2.0, 1e-6);

  // Check mean in middle bin
  auto idx2 = bh.axis(0).index(0.0);
  BOOST_CHECK_CLOSE(bh.at(idx2).value(), 6.0, 1e-6);

  // Check mean in last bin
  auto idx4 = bh.axis(0).index(2.0);
  BOOST_CHECK_CLOSE(bh.at(idx4).value(), 10.0, 1e-6);
}

BOOST_AUTO_TEST_CASE(Efficiency1D_BasicFill) {
  ProtoAxis protoAxis(AxisBoundaryType::Bound, -3.0, 3.0, 10);
  auto axis =
      AxisVariant(BoostVariableAxis(protoAxis.getAxis().getBinEdges(), "eta"));
  Efficiency1D eff("eff_vs_eta", "Efficiency vs Eta", {axis});

  BOOST_CHECK_EQUAL(eff.name(), "eff_vs_eta");
  BOOST_CHECK_EQUAL(eff.title(), "Efficiency vs Eta");

  // Fill eta=0.5: 7 accepted, 3 failed
  for (int i = 0; i < 7; ++i) {
    eff.fill({0.5}, true);
  }
  for (int i = 0; i < 3; ++i) {
    eff.fill({0.5}, false);
  }

  // Get efficiency for bin containing 0.5
  auto binIdx = eff.acceptedHistogram().axis(0).index(0.5);
  double accepted = static_cast<double>(eff.acceptedHistogram().at(binIdx));
  double total = static_cast<double>(eff.totalHistogram().at(binIdx));

  BOOST_CHECK_CLOSE(accepted, 7.0, 1e-10);
  BOOST_CHECK_CLOSE(total, 10.0, 1e-10);
  BOOST_CHECK_CLOSE(accepted / total, 0.7, 1e-6);
}

BOOST_AUTO_TEST_CASE(Efficiency1D_MultipleBins) {
  ProtoAxis protoAxis(AxisBoundaryType::Bound, 0.0, 5.0, 5);
  auto axis =
      AxisVariant(BoostVariableAxis(protoAxis.getAxis().getBinEdges(), "pt"));
  Efficiency1D eff("eff_vs_pt", "Efficiency vs pT", {axis});

  // Bin 0: 50% efficiency
  eff.fill({0.5}, true);
  eff.fill({0.5}, false);

  // Bin 2: 80% efficiency
  for (int i = 0; i < 8; ++i) {
    eff.fill({2.5}, true);
  }
  for (int i = 0; i < 2; ++i) {
    eff.fill({2.5}, false);
  }

  // Bin 4: 100% efficiency
  for (int i = 0; i < 5; ++i) {
    eff.fill({4.5}, true);
  }

  const auto& accepted = eff.acceptedHistogram();
  const auto& total = eff.totalHistogram();

  auto idx0 = accepted.axis(0).index(0.5);
  BOOST_CHECK_CLOSE(static_cast<double>(accepted.at(idx0)) /
                        static_cast<double>(total.at(idx0)),
                    0.5, 1e-6);

  auto idx2 = accepted.axis(0).index(2.5);
  BOOST_CHECK_CLOSE(static_cast<double>(accepted.at(idx2)) /
                        static_cast<double>(total.at(idx2)),
                    0.8, 1e-6);

  auto idx4 = accepted.axis(0).index(4.5);
  BOOST_CHECK_CLOSE(static_cast<double>(accepted.at(idx4)) /
                        static_cast<double>(total.at(idx4)),
                    1.0, 1e-6);
}

BOOST_AUTO_TEST_CASE(Efficiency2D_BasicFill) {
  ProtoAxis protoX(AxisBoundaryType::Bound, -2.5, 2.5, 5);
  ProtoAxis protoY(AxisBoundaryType::Bound, 0.0, 5.0, 5);
  auto xAxis =
      AxisVariant(BoostVariableAxis(protoX.getAxis().getBinEdges(), "eta"));
  auto yAxis =
      AxisVariant(BoostVariableAxis(protoY.getAxis().getBinEdges(), "pt"));
  Efficiency2D eff("eff_vs_eta_pt", "Efficiency vs Eta and pT", {xAxis, yAxis});

  BOOST_CHECK_EQUAL(eff.name(), "eff_vs_eta_pt");
  BOOST_CHECK_EQUAL(eff.title(), "Efficiency vs Eta and pT");

  // Fill (0.0, 2.5): 3 accepted, 1 failed
  eff.fill({0.0, 2.5}, true);
  eff.fill({0.0, 2.5}, true);
  eff.fill({0.0, 2.5}, true);
  eff.fill({0.0, 2.5}, false);

  const auto& accepted = eff.acceptedHistogram();
  const auto& total = eff.totalHistogram();

  auto xIdx = accepted.axis(0).index(0.0);
  auto yIdx = accepted.axis(1).index(2.5);

  double acceptedCount = static_cast<double>(accepted.at(xIdx, yIdx));
  double totalCount = static_cast<double>(total.at(xIdx, yIdx));

  BOOST_CHECK_CLOSE(acceptedCount, 3.0, 1e-10);
  BOOST_CHECK_CLOSE(totalCount, 4.0, 1e-10);
  BOOST_CHECK_CLOSE(acceptedCount / totalCount, 0.75, 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()
