// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "ActsExamples/Utilities/Helpers.hpp"
#include "ActsExamples/Validation/BoostHistogramWrappers.hpp"

#include <random>

using namespace ActsExamples;

BOOST_AUTO_TEST_SUITE(BoostProfileEfficiencySuite)

BOOST_AUTO_TEST_CASE(BoostProfileHistogram_BasicFill) {
  auto xBinning = PlotHelpers::Binning::Uniform("x", 10, 0.0, 10.0);
  BoostProfileHistogram profile("test_prof", "Test Profile", xBinning,
                                "y value");

  BOOST_CHECK_EQUAL(profile.name(), "test_prof");
  BOOST_CHECK_EQUAL(profile.title(), "Test Profile");
  BOOST_CHECK_EQUAL(profile.xAxisTitle(), "x");
  BOOST_CHECK_EQUAL(profile.yAxisTitle(), "y value");

  // Fill x=2.5 with y=10.0, y=20.0, y=30.0
  profile.fill(2.5, 10.0);
  profile.fill(2.5, 20.0);
  profile.fill(2.5, 30.0);

  const auto& bh = profile.histogram();
  auto xIdx = bh.axis(0).index(2.5);

  // Mean should be 20.0
  double mean = bh.at(xIdx).value();
  BOOST_CHECK_CLOSE(mean, 20.0, 1e-6);

  // Count should be 3
  double count = bh.at(xIdx).count();
  BOOST_CHECK_CLOSE(count, 3.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(BoostProfileHistogram_MultipleBins) {
  auto xBinning = PlotHelpers::Binning::Uniform("eta", 5, -2.5, 2.5);
  BoostProfileHistogram profile("res_vs_eta", "Residual vs Eta", xBinning,
                                "residual");

  // Fill different eta bins with different mean values
  profile.fill(-2.0, 1.0);  // bin 0
  profile.fill(-2.0, 3.0);
  profile.fill(0.0, 5.0);  // bin 2 (middle)
  profile.fill(0.0, 7.0);
  profile.fill(2.0, 9.0);  // bin 4
  profile.fill(2.0, 11.0);

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

BOOST_AUTO_TEST_CASE(BoostEfficiency1D_BasicFill) {
  auto binning = PlotHelpers::Binning::Uniform("eta", 10, -3.0, 3.0);
  BoostEfficiency1D eff("eff_vs_eta", "Efficiency vs Eta", binning);

  BOOST_CHECK_EQUAL(eff.name(), "eff_vs_eta");
  BOOST_CHECK_EQUAL(eff.title(), "Efficiency vs Eta");
  BOOST_CHECK_EQUAL(eff.axisTitle(), "eta");

  // Fill eta=0.5: 7 passed, 3 failed
  for (int i = 0; i < 7; ++i) {
    eff.fill(0.5, true);
  }
  for (int i = 0; i < 3; ++i) {
    eff.fill(0.5, false);
  }

  // Get efficiency for bin containing 0.5
  auto binIdx = eff.passedHistogram().axis(0).index(0.5);
  double passed = static_cast<double>(eff.passedHistogram().at(binIdx));
  double total = static_cast<double>(eff.totalHistogram().at(binIdx));

  BOOST_CHECK_CLOSE(passed, 7.0, 1e-10);
  BOOST_CHECK_CLOSE(total, 10.0, 1e-10);
  BOOST_CHECK_CLOSE(passed / total, 0.7, 1e-6);
}

BOOST_AUTO_TEST_CASE(BoostEfficiency1D_MultipleBins) {
  auto binning = PlotHelpers::Binning::Uniform("pt", 5, 0.0, 5.0);
  BoostEfficiency1D eff("eff_vs_pt", "Efficiency vs pT", binning);

  // Bin 0: 50% efficiency
  eff.fill(0.5, true);
  eff.fill(0.5, false);

  // Bin 2: 80% efficiency
  for (int i = 0; i < 8; ++i) {
    eff.fill(2.5, true);
  }
  for (int i = 0; i < 2; ++i) {
    eff.fill(2.5, false);
  }

  // Bin 4: 100% efficiency
  for (int i = 0; i < 5; ++i) {
    eff.fill(4.5, true);
  }

  const auto& passed = eff.passedHistogram();
  const auto& total = eff.totalHistogram();

  auto idx0 = passed.axis(0).index(0.5);
  BOOST_CHECK_CLOSE(static_cast<double>(passed.at(idx0)) /
                        static_cast<double>(total.at(idx0)),
                    0.5, 1e-6);

  auto idx2 = passed.axis(0).index(2.5);
  BOOST_CHECK_CLOSE(static_cast<double>(passed.at(idx2)) /
                        static_cast<double>(total.at(idx2)),
                    0.8, 1e-6);

  auto idx4 = passed.axis(0).index(4.5);
  BOOST_CHECK_CLOSE(static_cast<double>(passed.at(idx4)) /
                        static_cast<double>(total.at(idx4)),
                    1.0, 1e-6);
}

BOOST_AUTO_TEST_CASE(BoostEfficiency2D_BasicFill) {
  auto xBinning = PlotHelpers::Binning::Uniform("eta", 5, -2.5, 2.5);
  auto yBinning = PlotHelpers::Binning::Uniform("pt", 5, 0.0, 5.0);
  BoostEfficiency2D eff("eff_vs_eta_pt", "Efficiency vs Eta and pT", xBinning,
                        yBinning);

  BOOST_CHECK_EQUAL(eff.name(), "eff_vs_eta_pt");
  BOOST_CHECK_EQUAL(eff.title(), "Efficiency vs Eta and pT");
  BOOST_CHECK_EQUAL(eff.xAxisTitle(), "eta");
  BOOST_CHECK_EQUAL(eff.yAxisTitle(), "pt");

  // Fill (0.0, 2.5): 3 passed, 1 failed
  eff.fill(0.0, 2.5, true);
  eff.fill(0.0, 2.5, true);
  eff.fill(0.0, 2.5, true);
  eff.fill(0.0, 2.5, false);

  const auto& passed = eff.passedHistogram();
  const auto& total = eff.totalHistogram();

  auto xIdx = passed.axis(0).index(0.0);
  auto yIdx = passed.axis(1).index(2.5);

  double passedCount = static_cast<double>(passed.at(xIdx, yIdx));
  double totalCount = static_cast<double>(total.at(xIdx, yIdx));

  BOOST_CHECK_CLOSE(passedCount, 3.0, 1e-10);
  BOOST_CHECK_CLOSE(totalCount, 4.0, 1e-10);
  BOOST_CHECK_CLOSE(passedCount / totalCount, 0.75, 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()
