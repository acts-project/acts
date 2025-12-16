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

#include <TEfficiency.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>

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

BOOST_AUTO_TEST_CASE(Conversion_BoostProfile_to_TProfile) {
  auto xBinning = PlotHelpers::Binning::Uniform("eta", 10, -2.5, 2.5);
  BoostProfileHistogram profile("res_mean_vs_eta", "Mean Residual vs Eta",
                                xBinning, "residual [mm]");

  // Fill with known values to check mean calculation
  profile.fill(-2.0, 1.0);
  profile.fill(-2.0, 3.0);  // Mean = 2.0
  profile.fill(0.0, 5.0);
  profile.fill(0.0, 7.0);  // Mean = 6.0
  profile.fill(2.0, 9.0);
  profile.fill(2.0, 11.0);  // Mean = 10.0

  TProfile* rootProfile = BoostHistogramToRoot::toTProfile(profile);

  // Verify metadata
  BOOST_CHECK_EQUAL(std::string(rootProfile->GetName()), "res_mean_vs_eta");
  BOOST_CHECK_EQUAL(std::string(rootProfile->GetTitle()),
                    "Mean Residual vs Eta");
  BOOST_CHECK_EQUAL(std::string(rootProfile->GetXaxis()->GetTitle()), "eta");
  BOOST_CHECK_EQUAL(std::string(rootProfile->GetYaxis()->GetTitle()),
                    "residual [mm]");

  // Verify binning
  BOOST_CHECK_EQUAL(rootProfile->GetNbinsX(), 10);

  // Verify mean values in filled bins
  const auto& bh = profile.histogram();
  auto idx0 = bh.axis(0).index(-2.0);
  auto idx2 = bh.axis(0).index(0.0);
  auto idx4 = bh.axis(0).index(2.0);

  BOOST_CHECK_CLOSE(rootProfile->GetBinContent(idx0 + 1),
                    bh.at(idx0).value(), 1e-6);
  BOOST_CHECK_CLOSE(rootProfile->GetBinContent(idx2 + 1),
                    bh.at(idx2).value(), 1e-6);
  BOOST_CHECK_CLOSE(rootProfile->GetBinContent(idx4 + 1),
                    bh.at(idx4).value(), 1e-6);

  delete rootProfile;
}

BOOST_AUTO_TEST_CASE(Conversion_BoostEfficiency1D_to_TEfficiency) {
  auto binning = PlotHelpers::Binning::Uniform("eta", 10, -3.0, 3.0);
  BoostEfficiency1D eff("eff_vs_eta", "Efficiency vs Eta", binning);

  // Fill with known pass/fail patterns
  // Bin at eta=0.5: 7 passed, 3 failed (70% efficiency)
  for (int i = 0; i < 7; ++i) {
    eff.fill(0.5, true);
  }
  for (int i = 0; i < 3; ++i) {
    eff.fill(0.5, false);
  }

  // Bin at eta=-1.5: 5 passed, 5 failed (50% efficiency)
  for (int i = 0; i < 5; ++i) {
    eff.fill(-1.5, true);
    eff.fill(-1.5, false);
  }

  TEfficiency* rootEff = BoostHistogramToRoot::toTEfficiency(eff);

  // Verify metadata
  BOOST_CHECK_EQUAL(std::string(rootEff->GetName()), "eff_vs_eta");
  BOOST_CHECK_EQUAL(std::string(rootEff->GetTitle()), "Efficiency vs Eta");

  // Verify binning
  BOOST_CHECK_EQUAL(rootEff->GetTotalHistogram()->GetNbinsX(), 10);

  // Verify efficiency values
  const auto& passed = eff.passedHistogram();
  const auto& total = eff.totalHistogram();

  auto idx1 = passed.axis(0).index(0.5);
  auto idx2 = passed.axis(0).index(-1.5);

  double eff1 = static_cast<double>(passed.at(idx1)) /
                static_cast<double>(total.at(idx1));
  double eff2 = static_cast<double>(passed.at(idx2)) /
                static_cast<double>(total.at(idx2));

  BOOST_CHECK_CLOSE(rootEff->GetEfficiency(idx1 + 1), eff1, 1e-6);
  BOOST_CHECK_CLOSE(rootEff->GetEfficiency(idx2 + 1), eff2, 1e-6);

  delete rootEff;
}

BOOST_AUTO_TEST_CASE(Conversion_BoostEfficiency2D_to_TEfficiency) {
  auto xBinning = PlotHelpers::Binning::Uniform("eta", 5, -2.5, 2.5);
  auto yBinning = PlotHelpers::Binning::Uniform("pt", 5, 0.0, 5.0);
  BoostEfficiency2D eff("eff_vs_eta_pt", "Efficiency vs Eta and pT", xBinning,
                        yBinning);

  // Fill bin (0.0, 2.5): 3 passed, 1 failed (75% efficiency)
  eff.fill(0.0, 2.5, true);
  eff.fill(0.0, 2.5, true);
  eff.fill(0.0, 2.5, true);
  eff.fill(0.0, 2.5, false);

  // Fill bin (-1.5, 1.5): 2 passed, 2 failed (50% efficiency)
  eff.fill(-1.5, 1.5, true);
  eff.fill(-1.5, 1.5, true);
  eff.fill(-1.5, 1.5, false);
  eff.fill(-1.5, 1.5, false);

  TEfficiency* rootEff = BoostHistogramToRoot::toTEfficiency(eff);

  // Verify metadata
  BOOST_CHECK_EQUAL(std::string(rootEff->GetName()), "eff_vs_eta_pt");
  BOOST_CHECK_EQUAL(std::string(rootEff->GetTitle()),
                    "Efficiency vs Eta and pT");

  // Verify 2D binning
  BOOST_CHECK_EQUAL(rootEff->GetTotalHistogram()->GetNbinsX(), 5);
  BOOST_CHECK_EQUAL(rootEff->GetTotalHistogram()->GetNbinsY(), 5);

  // Verify efficiency values
  const auto& passed = eff.passedHistogram();
  const auto& total = eff.totalHistogram();

  auto xIdx1 = passed.axis(0).index(0.0);
  auto yIdx1 = passed.axis(1).index(2.5);
  auto xIdx2 = passed.axis(0).index(-1.5);
  auto yIdx2 = passed.axis(1).index(1.5);

  double eff1 = static_cast<double>(passed.at(xIdx1, yIdx1)) /
                static_cast<double>(total.at(xIdx1, yIdx1));
  double eff2 = static_cast<double>(passed.at(xIdx2, yIdx2)) /
                static_cast<double>(total.at(xIdx2, yIdx2));

  int globalBin1 =
      rootEff->GetGlobalBin(static_cast<int>(xIdx1 + 1), static_cast<int>(yIdx1 + 1));
  int globalBin2 =
      rootEff->GetGlobalBin(static_cast<int>(xIdx2 + 1), static_cast<int>(yIdx2 + 1));

  BOOST_CHECK_CLOSE(rootEff->GetEfficiency(globalBin1), eff1, 1e-6);
  BOOST_CHECK_CLOSE(rootEff->GetEfficiency(globalBin2), eff2, 1e-6);

  delete rootEff;
}

BOOST_AUTO_TEST_SUITE_END()
