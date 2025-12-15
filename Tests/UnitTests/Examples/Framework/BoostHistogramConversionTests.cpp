// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "ActsExamples/Io/Root/BoostHistogramToRootConverter.hpp"
#include "ActsExamples/Utilities/Helpers.hpp"
#include "ActsExamples/Validation/BoostHistogramWrappers.hpp"

#include <TH1F.h>
#include <TH2F.h>

#include <cmath>
#include <random>
#include <string>
#include <vector>

using namespace ActsExamples;

BOOST_AUTO_TEST_SUITE(BoostHistogramConversionSuite)

BOOST_AUTO_TEST_CASE(Conversion_Boost1D_to_ROOT_UniformBinning) {
  auto binning = PlotHelpers::Binning::Uniform("x [cm]", 10, 0.0, 10.0);
  BoostHistogram1D boostHist("test_hist", "Test Histogram", binning);

  // Fill with 100 random values
  std::mt19937 rng(42);
  std::normal_distribution<double> dist(5.0, 1.0);
  for (int i = 0; i < 100; ++i) {
    boostHist.fill(dist(rng));
  }

  TH1F* rootHist = BoostHistogramToRoot::toTH1F(boostHist);

  // Verify metadata
  BOOST_CHECK_EQUAL(std::string(rootHist->GetName()), "test_hist");
  BOOST_CHECK_EQUAL(std::string(rootHist->GetTitle()), "Test Histogram");
  BOOST_CHECK_EQUAL(std::string(rootHist->GetXaxis()->GetTitle()), "x [cm]");

  // Verify binning
  BOOST_CHECK_EQUAL(rootHist->GetNbinsX(), 10);

  // Verify bin contents match
  const auto& bh = boostHist.histogram();
  for (int i = 1; i <= rootHist->GetNbinsX(); ++i) {
    BOOST_CHECK_CLOSE(rootHist->GetBinContent(i), static_cast<double>(bh.at(i - 1)), 1e-10);
  }

  delete rootHist;
}

BOOST_AUTO_TEST_CASE(Conversion_Boost1D_to_ROOT_VariableBinning) {
  std::vector<double> edges = {0.0, 0.5, 1.0, 2.0, 5.0, 10.0};
  auto binning = PlotHelpers::Binning::Variable("eta", edges);
  BoostHistogram1D boostHist("res_eta", "Residual vs Eta", binning);

  // Fill bins with different counts
  boostHist.fill(0.3);
  boostHist.fill(0.3);
  boostHist.fill(1.5);
  boostHist.fill(3.0);

  TH1F* rootHist = BoostHistogramToRoot::toTH1F(boostHist);

  BOOST_CHECK_EQUAL(rootHist->GetNbinsX(), 5);

  // Verify bin contents match
  const auto& bh = boostHist.histogram();
  for (int i = 1; i <= rootHist->GetNbinsX(); ++i) {
    BOOST_CHECK_CLOSE(rootHist->GetBinContent(i), static_cast<double>(bh.at(i - 1)), 1e-10);
  }

  delete rootHist;
}

BOOST_AUTO_TEST_CASE(Conversion_Boost1D_to_ROOT_LogarithmicBinning) {
  auto binning = PlotHelpers::Binning::Logarithmic("pT [GeV]", 10, 0.1, 100.0);
  BoostHistogram1D boostHist("pt_log", "pT Logarithmic", binning);

  // Fill with 100 log-distributed values
  std::mt19937 rng(123);
  std::uniform_real_distribution<double> logDist(std::log(0.1),
                                                  std::log(100.0));
  for (int i = 0; i < 100; ++i) {
    boostHist.fill(std::exp(logDist(rng)));
  }

  TH1F* rootHist = BoostHistogramToRoot::toTH1F(boostHist);

  BOOST_CHECK_EQUAL(rootHist->GetNbinsX(), 10);

  // Verify bin contents match
  const auto& bh = boostHist.histogram();
  for (int i = 1; i <= rootHist->GetNbinsX(); ++i) {
    BOOST_CHECK_CLOSE(rootHist->GetBinContent(i), static_cast<double>(bh.at(i - 1)), 1e-10);
  }

  delete rootHist;
}

BOOST_AUTO_TEST_CASE(Conversion_Boost2D_to_ROOT) {
  auto xBinning = PlotHelpers::Binning::Uniform("eta", 10, -2.5, 2.5);
  auto yBinning = PlotHelpers::Binning::Uniform("residual", 10, -5.0, 5.0);
  BoostHistogram2D boostHist("res_vs_eta", "Residual vs Eta", xBinning,
                             yBinning);

  // Fill with 100 2D Gaussian values
  std::mt19937 rng(456);
  std::normal_distribution<double> etaDist(0.0, 1.0);
  std::normal_distribution<double> resDist(0.0, 2.0);
  for (int i = 0; i < 100; ++i) {
    boostHist.fill(etaDist(rng), resDist(rng));
  }

  TH2F* rootHist = BoostHistogramToRoot::toTH2F(boostHist);

  // Verify metadata
  BOOST_CHECK_EQUAL(std::string(rootHist->GetName()), "res_vs_eta");
  BOOST_CHECK_EQUAL(std::string(rootHist->GetXaxis()->GetTitle()), "eta");
  BOOST_CHECK_EQUAL(std::string(rootHist->GetYaxis()->GetTitle()), "residual");

  // Verify binning and contents
  BOOST_CHECK_EQUAL(rootHist->GetNbinsX(), 10);
  BOOST_CHECK_EQUAL(rootHist->GetNbinsY(), 10);

  const auto& bh = boostHist.histogram();
  for (int i = 1; i <= rootHist->GetNbinsX(); ++i) {
    for (int j = 1; j <= rootHist->GetNbinsY(); ++j) {
      BOOST_CHECK_CLOSE(rootHist->GetBinContent(i, j), static_cast<double>(bh.at(i - 1, j - 1)),
                        1e-10);
    }
  }

  delete rootHist;
}

BOOST_AUTO_TEST_CASE(Conversion_Boost2D_to_ROOT_VariableBinning) {
  std::vector<double> etaEdges = {-2.5, -1.5, -0.5, 0.5, 1.5, 2.5};
  std::vector<double> ptEdges = {0.5, 1.0, 2.0, 5.0, 10.0};

  auto xBinning = PlotHelpers::Binning::Variable("eta", etaEdges);
  auto yBinning = PlotHelpers::Binning::Variable("pT [GeV]", ptEdges);
  BoostHistogram2D boostHist("eff_vs_eta_pt", "Efficiency vs Eta and pT",
                             xBinning, yBinning);

  boostHist.fill(0.0, 3.0);
  boostHist.fill(-1.0, 7.0);
  boostHist.fill(2.0, 1.5);

  TH2F* rootHist = BoostHistogramToRoot::toTH2F(boostHist);

  BOOST_CHECK_EQUAL(rootHist->GetNbinsX(), 5);
  BOOST_CHECK_EQUAL(rootHist->GetNbinsY(), 4);

  // Verify bin contents match
  const auto& bh = boostHist.histogram();
  for (int i = 1; i <= rootHist->GetNbinsX(); ++i) {
    for (int j = 1; j <= rootHist->GetNbinsY(); ++j) {
      BOOST_CHECK_CLOSE(rootHist->GetBinContent(i, j), static_cast<double>(bh.at(i - 1, j - 1)),
                        1e-10);
    }
  }

  delete rootHist;
}

BOOST_AUTO_TEST_CASE(Conversion_EmptyHistogram) {
  auto binning = PlotHelpers::Binning::Uniform("x", 10, -10.0, 10.0);
  BoostHistogram1D boostHist("empty", "Empty Histogram", binning);

  TH1F* rootHist = BoostHistogramToRoot::toTH1F(boostHist);

  BOOST_CHECK_EQUAL(rootHist->GetNbinsX(), 10);
  for (int i = 1; i <= rootHist->GetNbinsX(); ++i) {
    BOOST_CHECK_EQUAL(rootHist->GetBinContent(i), 0.0);
  }

  delete rootHist;
}

BOOST_AUTO_TEST_CASE(Conversion_MultipleFills) {
  auto binning = PlotHelpers::Binning::Uniform("x", 10, 0.0, 10.0);
  BoostHistogram1D boostHist("multifill", "Multiple Fills", binning);

  boostHist.fill(5.0);
  boostHist.fill(5.0);
  boostHist.fill(5.0);

  TH1F* rootHist = BoostHistogramToRoot::toTH1F(boostHist);

  auto binIdx = boostHist.histogram().axis(0).index(5.0);
  BOOST_CHECK_CLOSE(rootHist->GetBinContent(binIdx + 1), 3.0, 1e-10);

  delete rootHist;
}

BOOST_AUTO_TEST_CASE(Conversion_MetadataPreservation) {
  auto binning = PlotHelpers::Binning::Uniform("Distance [mm]", 10, 0.0, 1000.0);
  BoostHistogram1D boostHist("track_distance", "Track Distance", binning);

  boostHist.fill(500.0);

  TH1F* rootHist = BoostHistogramToRoot::toTH1F(boostHist);

  BOOST_CHECK_EQUAL(std::string(rootHist->GetName()), "track_distance");
  BOOST_CHECK_EQUAL(std::string(rootHist->GetTitle()), "Track Distance");
  BOOST_CHECK_EQUAL(std::string(rootHist->GetXaxis()->GetTitle()),
                    "Distance [mm]");

  delete rootHist;
}

BOOST_AUTO_TEST_SUITE_END()
