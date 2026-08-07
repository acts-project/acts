// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Exercises the backend-agnostic fit interface (HistogramFitFunction,
// iterativeFit, extractMeanWidthProfiles) with ActsPlugins::RootHistogramFit
// as a concrete backend -- the only fit backend available in this checkout.
// A second backend (e.g. a Python callable, or a future C++ implementation)
// can be plugged in through the same HistogramFitFunction signature without
// any of this code changing.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/Histogram.hpp"
#include "ActsExamples/Validation/IterativeFit.hpp"
#include "ActsPlugins/Root/RootHistogramFit.hpp"

#include <cmath>
#include <optional>
#include <random>
#include <string>

using namespace Acts;
using namespace Acts::Experimental;
using namespace ActsExamples;

namespace {

const ActsPlugins::RootHistogramFit rootFitter;
const HistogramFitFunction rootFn = [](const Histogram1& hist,
                                       std::optional<HistogramFitRange> range) {
  return rootFitter.fit(hist, range);
};

/// Sample `count` entries from N(mean, sigma) into a fresh histogram
Histogram1 sampled(const std::string& name, std::size_t count, double mean,
                   double sigma, int nBins, double xMin, double xMax,
                   std::uint32_t seed) {
  auto axis = AxisVariant(BoostRegularAxis(nBins, xMin, xMax, "x"));
  Histogram1 hist(name, name, {axis});

  std::mt19937 generator(seed);
  std::normal_distribution<double> distribution(mean, sigma);
  for (std::size_t i = 0; i < count; ++i) {
    hist.fill({distribution(generator)});
  }

  return hist;
}

}  // namespace

BOOST_AUTO_TEST_SUITE(HistogramFitInterfaceSuite)

BOOST_AUTO_TEST_CASE(IterativeFit_RecoversKnownParameters) {
  const Histogram1 hist =
      sampled("iterative", 20000, 0.3, 0.7, 60, -5.0, 5.0, 102);

  const auto result = iterativeFit(rootFn, hist, /*sigmaRange=*/3.0,
                                   /*iterations=*/3);
  BOOST_REQUIRE(result.has_value());

  const auto& [mean, sigma, meanError, sigmaError] = *result;
  BOOST_CHECK_CLOSE(mean, 0.3, 5.0);
  BOOST_CHECK_CLOSE(sigma, 0.7, 5.0);
  BOOST_CHECK_GT(meanError, 0.0);
  BOOST_CHECK_GT(sigmaError, 0.0);
}

// ActsPlugins::RootHistogramFit is otherwise uncalled from generic Examples
// code, so this is the test proving that the generic profile extraction
// really works with a HistogramFitFunction plugged in from outside Examples,
// not just with a fit function Examples defines itself.
BOOST_AUTO_TEST_CASE(ExtractMeanWidthProfiles_WorksWithRootFitter) {
  const int nEtaBins = 3;
  auto etaAxis = AxisVariant(BoostRegularAxis(nEtaBins, 0.0, 3.0, "eta"));
  auto resAxis = AxisVariant(BoostRegularAxis(80, -10.0, 10.0, "res"));
  Histogram2 hist("res_vs_eta", "Residual vs Eta", {etaAxis, resAxis});

  std::mt19937 generator(4242);
  for (int i = 0; i < nEtaBins; ++i) {
    std::normal_distribution<double> distribution(0.0, 0.5 + 0.3 * i);
    for (int n = 0; n < 20000; ++n) {
      hist.fill({i + 0.5, distribution(generator)});
    }
  }

  const auto profiles = extractMeanWidthProfiles(rootFn, hist, "mean", "width");

  BOOST_CHECK_EQUAL(profiles.fitFailureFraction, 0.0);
  for (int i = 0; i < nEtaBins; ++i) {
    BOOST_CHECK_CLOSE(profiles.width.value({i}), 0.5 + 0.3 * i, 5.0);
  }
}

BOOST_AUTO_TEST_SUITE_END()
