// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/ProtoAxis.hpp"
#include "ActsPlugins/Root/HistogramConverter.hpp"

#include <cmath>
#include <random>
#include <string>
#include <vector>

#include <TEfficiency.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>

using namespace Acts;
using namespace Acts::Experimental;
using namespace ActsPlugins;

BOOST_AUTO_TEST_SUITE(HistogramConversionSuite)

namespace {
/// Helper function to test 1D histogram conversion
void testHist1D(const Histogram1& boostHist) {
  TH1F* rootHist = toRoot(boostHist);

  // Verify metadata
  BOOST_CHECK_EQUAL(std::string(rootHist->GetName()), boostHist.name());
  BOOST_CHECK_EQUAL(std::string(rootHist->GetTitle()), boostHist.title());

  // Verify number of bins
  const auto& bh = boostHist.histogram();
  BOOST_CHECK_EQUAL(rootHist->GetNbinsX(), static_cast<int>(bh.axis(0).size()));

  // Verify bin contents match
  for (int i = 1; i <= rootHist->GetNbinsX(); ++i) {
    BOOST_CHECK_CLOSE(rootHist->GetBinContent(i),
                      static_cast<double>(bh.at(i - 1)), 1e-10);
  }

  delete rootHist;
}

/// Helper function to test 2D histogram conversion
void testHist2D(const Histogram2& boostHist) {
  TH2F* rootHist = toRoot(boostHist);

  // Verify metadata
  BOOST_CHECK_EQUAL(std::string(rootHist->GetName()), boostHist.name());
  BOOST_CHECK_EQUAL(std::string(rootHist->GetTitle()), boostHist.title());

  // Verify number of bins
  const auto& bh = boostHist.histogram();
  BOOST_CHECK_EQUAL(rootHist->GetNbinsX(), static_cast<int>(bh.axis(0).size()));
  BOOST_CHECK_EQUAL(rootHist->GetNbinsY(), static_cast<int>(bh.axis(1).size()));

  // Verify bin contents match
  for (int i = 1; i <= rootHist->GetNbinsX(); ++i) {
    for (int j = 1; j <= rootHist->GetNbinsY(); ++j) {
      BOOST_CHECK_CLOSE(rootHist->GetBinContent(i, j),
                        static_cast<double>(bh.at(i - 1, j - 1)), 1e-10);
    }
  }

  delete rootHist;
}
}  // namespace

BOOST_AUTO_TEST_CASE(Conversion_EmptyHistogram) {
  ProtoAxis protoAxis(AxisBoundaryType::Bound, -10.0, 10.0, 10);
  auto axis =
      AxisVariant(BoostVariableAxis(protoAxis.getAxis().getBinEdges(), "x"));
  Histogram1 boostHist("empty", "Empty Histogram", {axis});

  TH1F* rootHist = toRoot(boostHist);

  BOOST_CHECK_EQUAL(rootHist->GetNbinsX(), 10);
  for (int i = 1; i <= rootHist->GetNbinsX(); ++i) {
    BOOST_CHECK_EQUAL(rootHist->GetBinContent(i), 0.0);
  }

  delete rootHist;
}

BOOST_AUTO_TEST_CASE(Conversion_Boost1D_to_ROOT_UniformBinning) {
  ProtoAxis protoAxis(AxisBoundaryType::Bound, 0.0, 10.0, 10);
  auto axis = AxisVariant(
      BoostVariableAxis(protoAxis.getAxis().getBinEdges(), "x [cm]"));
  Histogram1 boostHist("test_hist", "Test Histogram", {axis});

  // Fill with 100 random values
  std::mt19937 rng(42);
  std::normal_distribution<double> dist(5.0, 1.0);
  for (int i = 0; i < 100; ++i) {
    boostHist.fill({dist(rng)});
  }

  testHist1D(boostHist);
}

BOOST_AUTO_TEST_CASE(Conversion_Boost1D_to_ROOT_VariableBinning) {
  std::vector<double> edges = {0.0, 0.5, 1.0, 2.0, 5.0, 10.0};
  auto axis = BoostVariableAxis(edges, "eta");
  Histogram1 boostHist("res_eta", "Residual vs Eta", {axis});

  // Fill bins with different counts
  boostHist.fill({0.3});
  boostHist.fill({0.3});
  boostHist.fill({1.5});
  boostHist.fill({3.0});

  testHist1D(boostHist);
}

BOOST_AUTO_TEST_CASE(Conversion_Boost2D_to_ROOT) {
  ProtoAxis protoX(AxisBoundaryType::Bound, -2.5, 2.5, 10);
  ProtoAxis protoY(AxisBoundaryType::Bound, -5.0, 5.0, 10);
  auto xAxis = BoostVariableAxis(protoX.getAxis().getBinEdges(), "eta");
  auto yAxis = BoostVariableAxis(protoY.getAxis().getBinEdges(), "residual");
  Histogram2 boostHist("res_vs_eta", "Residual vs Eta", {xAxis, yAxis});

  // Fill with 100 2D Gaussian values
  std::mt19937 rng(456);
  std::normal_distribution<double> etaDist(0.0, 1.0);
  std::normal_distribution<double> resDist(0.0, 2.0);
  for (int i = 0; i < 100; ++i) {
    boostHist.fill({etaDist(rng), resDist(rng)});
  }

  testHist2D(boostHist);
}

BOOST_AUTO_TEST_CASE(Conversion_Boost2D_to_ROOT_VariableBinning) {
  std::vector<double> etaEdges = {-2.5, -1.5, -0.5, 0.5, 1.5, 2.5};
  std::vector<double> ptEdges = {0.5, 1.0, 2.0, 5.0, 10.0};

  auto xAxis = BoostVariableAxis(etaEdges, "eta");
  auto yAxis = BoostVariableAxis(ptEdges, "pT [GeV]");
  Histogram2 boostHist("eff_vs_eta_pt", "Efficiency vs Eta and pT",
                       {xAxis, yAxis});

  boostHist.fill({0.0, 3.0});
  boostHist.fill({-1.0, 7.0});
  boostHist.fill({2.0, 1.5});

  testHist2D(boostHist);
}

BOOST_AUTO_TEST_CASE(Conversion_BoostProfile_to_TProfile) {
  ProtoAxis protoAxis(AxisBoundaryType::Bound, -2.5, 2.5, 10);
  auto xAxis = BoostVariableAxis(protoAxis.getAxis().getBinEdges(), "eta");
  ProfileHistogram1 profile("res_mean_vs_eta", "Mean Residual vs Eta", {xAxis},
                            "residual [mm]");

  std::map<double, double> expectedMeans;

  // Fill with known values to check mean calculation
  profile.fill({-2.0}, 1.0);
  profile.fill({-2.0}, 3.0);
  expectedMeans[-2.0] = 2.0;

  profile.fill({0.0}, 5.0);
  profile.fill({0.0}, 7.0);
  expectedMeans[0.0] = 6.0;

  profile.fill({2.0}, 9.0);
  profile.fill({2.0}, 11.0);
  profile.fill({2.0}, 13.0);
  expectedMeans[2.0] = 11.0;

  TProfile* rootProfile = toRoot(profile);

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

  for (auto x : {-2.0, 0.0, 2.0}) {
    auto binIdx = bh.axis(0).index(x);
    BOOST_CHECK_CLOSE(rootProfile->GetBinContent(binIdx + 1),
                      bh.at(binIdx).value(), 1e-6);
    BOOST_CHECK_CLOSE(rootProfile->GetBinContent(binIdx + 1),
                      expectedMeans.at(x), 1e-6);
  }

  delete rootProfile;
}

BOOST_AUTO_TEST_CASE(Conversion_BoostProfile_to_TProfile_WithErrors) {
  ProtoAxis protoAxis(AxisBoundaryType::Bound, -2.5, 2.5, 5);
  auto xAxis = BoostVariableAxis(protoAxis.getAxis().getBinEdges(), "eta");

  // Create ACTS ProfileHistogram
  ProfileHistogram1 actsProfile("profile_test", "Profile Test", {xAxis},
                                "value [units]");

  // Create ROOT TProfile with identical binning for comparison
  const auto& bh = actsProfile.histogram();
  const auto& axis = bh.axis(0);
  std::vector<double> edges(axis.size() + 1);
  for (int i = 0; i < axis.size(); ++i) {
    edges[i] = axis.bin(i).lower();
  }
  edges.back() = axis.bin(axis.size() - 1).upper();
  TProfile rootProfile("root_profile", "ROOT Profile", 5, edges.data());
  rootProfile.Sumw2();

  // Bin 1: mean=12, n=3
  for (auto y : {10.0, 12.0, 14.0}) {
    rootProfile.Fill(-2.0, y);
    actsProfile.fill({-2.0}, y);
  }

  // Bin 2: mean=6, n=2
  for (auto y : {5.0, 7.0}) {
    rootProfile.Fill(0.0, y);
    actsProfile.fill({0.0}, y);
  }

  // Bin 3: mean=23, n=4
  for (auto y : {20.0, 22.0, 24.0, 26.0}) {
    rootProfile.Fill(1.5, y);
    actsProfile.fill({1.5}, y);
  }

  // Convert ACTS profile to ROOT
  TProfile* convertedProfile = toRoot(actsProfile);

  // Compare binning
  BOOST_CHECK_EQUAL(convertedProfile->GetNbinsX(), rootProfile.GetNbinsX());

  // Compare each filled bin
  for (int i = 0; i < axis.size(); ++i) {
    int rootBin = i + 1;
    const auto& acc = bh.at(i);

    // Check count
    BOOST_CHECK_EQUAL(convertedProfile->GetBinEntries(rootBin),
                      rootProfile.GetBinEntries(rootBin));
    BOOST_CHECK_CLOSE(convertedProfile->GetBinEntries(rootBin), acc.count(),
                      1e-10);

    // Check mean value
    BOOST_CHECK_CLOSE(convertedProfile->GetBinContent(rootBin),
                      rootProfile.GetBinContent(rootBin), 1e-6);
    BOOST_CHECK_CLOSE(convertedProfile->GetBinContent(rootBin), acc.value(),
                      1e-6);

    // Check errors
    double rootError = rootProfile.GetBinError(rootBin);
    double convertedError = convertedProfile->GetBinError(rootBin);
    BOOST_CHECK_CLOSE(convertedError, rootError, 1e-6);
  }

  delete convertedProfile;
}

BOOST_AUTO_TEST_CASE(Conversion_Efficiency1D_to_TEfficiency) {
  ProtoAxis protoAxis(AxisBoundaryType::Bound, -3.0, 3.0, 10);
  auto axis = BoostVariableAxis(protoAxis.getAxis().getBinEdges(), "eta");
  Efficiency1 eff("eff_vs_eta", "Efficiency vs Eta", {axis});

  // Fill with known pass/fail patterns
  // Bin at eta=0.5: 3 accepted, 1 failed (75% efficiency)
  for (auto v : {true, true, true, false}) {
    eff.fill({0.5}, v);
  }

  // Bin at eta=-1.5: 3 accepted, 3 failed (50% efficiency)
  for (auto v : {true, true, true, false, false, false}) {
    eff.fill({-1.5}, v);
  }

  TEfficiency* rootEff = toRoot(eff);

  // Verify metadata
  BOOST_CHECK_EQUAL(std::string(rootEff->GetName()), "eff_vs_eta");
  BOOST_CHECK_EQUAL(std::string(rootEff->GetTitle()), "Efficiency vs Eta");

  // Verify binning
  BOOST_CHECK_EQUAL(rootEff->GetTotalHistogram()->GetNbinsX(), 10);

  // Verify efficiency values
  const auto& accepted = eff.acceptedHistogram();
  const auto& total = eff.totalHistogram();

  for (auto x : {0.5, -1.5}) {
    auto binIdx = accepted.axis(0).index(x);
    double expectedEff = static_cast<double>(accepted.at(binIdx)) /
                         static_cast<double>(total.at(binIdx));
    BOOST_CHECK_CLOSE(rootEff->GetEfficiency(binIdx + 1), expectedEff, 1e-6);
  }

  delete rootEff;
}

BOOST_AUTO_TEST_CASE(Conversion_Efficiency2D_to_TEfficiency) {
  ProtoAxis protoX(AxisBoundaryType::Bound, -2.5, 2.5, 5);
  ProtoAxis protoY(AxisBoundaryType::Bound, 0.0, 5.0, 5);
  auto xAxis = BoostVariableAxis(protoX.getAxis().getBinEdges(), "eta");
  auto yAxis = BoostVariableAxis(protoY.getAxis().getBinEdges(), "pt");
  Efficiency2 eff("eff_vs_eta_pt", "Efficiency vs Eta and pT", {xAxis, yAxis});

  // Fill bin (0.0, 2.5): 3 accepted, 1 failed (75% efficiency)
  for (auto v : {true, true, true, false}) {
    eff.fill({0.0, 2.5}, v);
  }

  // Fill bin (-1.5, 1.5): 3 accepted, 3 failed (50% efficiency)
  for (auto v : {true, true, true, false, false, false}) {
    eff.fill({-1.5, 1.5}, v);
  }

  TEfficiency* rootEff = toRoot(eff);

  // Verify metadata
  BOOST_CHECK_EQUAL(std::string(rootEff->GetName()), "eff_vs_eta_pt");
  BOOST_CHECK_EQUAL(std::string(rootEff->GetTitle()),
                    "Efficiency vs Eta and pT");

  // Verify 2D binning
  BOOST_CHECK_EQUAL(rootEff->GetTotalHistogram()->GetNbinsX(), 5);
  BOOST_CHECK_EQUAL(rootEff->GetTotalHistogram()->GetNbinsY(), 5);

  // Verify efficiency values
  const auto& accepted = eff.acceptedHistogram();
  const auto& total = eff.totalHistogram();

  using P = std::pair<double, double>;
  for (auto coords : {P{0.0, 2.5}, P{-1.5, 1.5}}) {
    auto xIdx = accepted.axis(0).index(coords.first);
    auto yIdx = accepted.axis(1).index(coords.second);
    double expectedEff = static_cast<double>(accepted.at(xIdx, yIdx)) /
                         static_cast<double>(total.at(xIdx, yIdx));
    int globalBin = rootEff->GetGlobalBin(static_cast<int>(xIdx + 1),
                                          static_cast<int>(yIdx + 1));
    BOOST_CHECK_CLOSE(rootEff->GetEfficiency(globalBin), expectedEff, 1e-6);
  }

  delete rootEff;
}

BOOST_AUTO_TEST_SUITE_END()
