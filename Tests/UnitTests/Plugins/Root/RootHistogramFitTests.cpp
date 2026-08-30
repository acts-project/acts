// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/Histogram.hpp"
#include "ActsPlugins/Root/HistogramConverter.hpp"
#include "ActsPlugins/Root/RootHistogramFit.hpp"

#include <cstdint>
#include <random>
#include <string>
#include <vector>

#include <TH1.h>
#include <TH2.h>

using namespace Acts::Experimental;
using ActsPlugins::RootHistogramFit;
using ActsPlugins::toRoot;

namespace {

Histogram1 makeHistogram(const std::string& name, int nBins, double xMin,
                         double xMax) {
  auto axis = AxisVariant(BoostRegularAxis(nBins, xMin, xMax, "x"));
  return Histogram1(name, name, {axis});
}

/// Sample `count` entries from N(mean, sigma) into a fresh histogram
Histogram1 sampled(const std::string& name, std::size_t count, double mean,
                   double sigma, int nBins, double xMin, double xMax,
                   std::uint32_t seed) {
  Histogram1 hist = makeHistogram(name, nBins, xMin, xMax);

  std::mt19937 generator(seed);
  std::normal_distribution<double> distribution(mean, sigma);
  for (std::size_t i = 0; i < count; ++i) {
    hist.fill({distribution(generator)});
  }

  return hist;
}

}  // namespace

BOOST_AUTO_TEST_SUITE(RootHistogramFitSuite)

BOOST_AUTO_TEST_CASE(SingleFit_RecoversKnownParameters) {
  const RootHistogramFit fitter;
  const Histogram1 hist =
      sampled("single", 100000, 0.0, 1.0, 100, -8.0, 8.0, 101);

  const auto result = fitter(hist);
  BOOST_REQUIRE(result.has_value());

  const auto& [mean, sigma, meanError, sigmaError] = *result;
  BOOST_CHECK_SMALL(mean, 0.05);
  BOOST_CHECK_CLOSE(sigma, 1.0, 5.0);
  BOOST_CHECK_GT(meanError, 0.0);
  BOOST_CHECK_GT(sigmaError, 0.0);
}

BOOST_AUTO_TEST_CASE(RestrictedRange_Works) {
  const RootHistogramFit fitter;
  const Histogram1 hist =
      sampled("restricted", 30000, 0.4, 1.2, 80, -8.0, 8.0, 202);

  const auto result = fitter(hist, RootHistogramFit::Range{-2.6, 3.4});
  BOOST_REQUIRE(result.has_value());

  const auto& [mean, sigma, meanError, sigmaError] = *result;
  BOOST_CHECK_CLOSE(mean, 0.4, 10.0);
  BOOST_CHECK_CLOSE(sigma, 1.2, 10.0);
}

BOOST_AUTO_TEST_CASE(DegenerateInputs_Fail) {
  const RootHistogramFit fitter;
  auto axis = AxisVariant(BoostRegularAxis(20, -5.0, 5.0, "x"));

  Histogram1 empty("empty", "empty", {axis});
  BOOST_CHECK(!fitter(empty).has_value());

  Histogram1 spike("spike", "spike", {axis});
  spike.setBinContent({10}, 500.0);
  BOOST_CHECK(!fitter(spike).has_value());
}

BOOST_AUTO_TEST_CASE(Histogram1D_ConvertsWithErrors) {
  std::vector<double> edges = {0.0, 1.0, 3.0, 7.0};
  auto axis = AxisVariant(BoostVariableAxis(edges, "eta"));
  Histogram1 hist("resmean_d0_vs_eta", "Mean", {axis});

  hist.setBin({0}, 1.5, 0.25);
  hist.setBin({1}, -2.5, 0.5);
  // Bin 2 deliberately left untouched

  const auto rootHist = toRoot(hist);
  BOOST_REQUIRE(rootHist != nullptr);
  BOOST_CHECK_EQUAL(rootHist->GetName(), "resmean_d0_vs_eta");
  BOOST_CHECK_EQUAL(rootHist->GetTitle(), "Mean");
  BOOST_CHECK_EQUAL(rootHist->GetNbinsX(), 3);
  BOOST_CHECK_EQUAL(std::string(rootHist->GetXaxis()->GetTitle()), "eta");

  BOOST_CHECK_CLOSE(rootHist->GetBinContent(1), 1.5, 1e-4);
  BOOST_CHECK_CLOSE(rootHist->GetBinError(1), 0.25, 1e-4);
  BOOST_CHECK_CLOSE(rootHist->GetBinContent(2), -2.5, 1e-4);
  BOOST_CHECK_CLOSE(rootHist->GetBinError(2), 0.5, 1e-4);
  BOOST_CHECK_EQUAL(rootHist->GetBinContent(3), 0.0);
  BOOST_CHECK_EQUAL(rootHist->GetBinError(3), 0.0);

  // Variable binning must be carried over
  BOOST_CHECK_CLOSE(rootHist->GetXaxis()->GetBinLowEdge(1), 0.0, 1e-6);
  BOOST_CHECK_CLOSE(rootHist->GetXaxis()->GetBinUpEdge(3), 7.0, 1e-6);
}

BOOST_AUTO_TEST_CASE(Histogram2D_ConvertsWithErrors) {
  auto xAxis = AxisVariant(BoostRegularAxis(2, 0.0, 2.0, "eta"));
  auto yAxis = AxisVariant(BoostRegularAxis(3, 0.0, 3.0, "pt"));
  Histogram2 hist("reswidth_d0_vs_eta_pt", "Width", {xAxis, yAxis});

  hist.setBin({1, 2}, 0.75, 0.1);

  const auto rootHist = toRoot(hist);
  BOOST_REQUIRE(rootHist != nullptr);
  BOOST_CHECK_EQUAL(rootHist->GetNbinsX(), 2);
  BOOST_CHECK_EQUAL(rootHist->GetNbinsY(), 3);
  BOOST_CHECK_EQUAL(std::string(rootHist->GetXaxis()->GetTitle()), "eta");
  BOOST_CHECK_EQUAL(std::string(rootHist->GetYaxis()->GetTitle()), "pt");

  BOOST_CHECK_CLOSE(rootHist->GetBinContent(2, 3), 0.75, 1e-4);
  BOOST_CHECK_CLOSE(rootHist->GetBinError(2, 3), 0.1, 1e-4);
  BOOST_CHECK_EQUAL(rootHist->GetBinContent(1, 1), 0.0);
}

BOOST_AUTO_TEST_SUITE_END()
