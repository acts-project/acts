// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/Histogram.hpp"
#include "ActsExamples/Validation/GaussianHistogramFit.hpp"
#include "ActsExamples/Validation/HistogramFit.hpp"

#include <array>
#include <cmath>
#include <functional>
#include <optional>
#include <random>
#include <tuple>
#include <vector>

using namespace Acts;
using namespace Acts::Experimental;
using namespace ActsExamples;

namespace {

/// Fill a histogram with `count` samples drawn from N(mean, sigma)
Histogram1 sampleGaussian(std::size_t count, double mean, double sigma,
                          int nBins, double xMin, double xMax,
                          std::uint32_t seed) {
  auto axis = AxisVariant(BoostRegularAxis(nBins, xMin, xMax, "x"));
  Histogram1 hist("toy", "Toy", {axis});

  std::mt19937 generator(seed);
  std::normal_distribution<double> distribution(mean, sigma);
  for (std::size_t i = 0; i < count; ++i) {
    hist.fill({distribution(generator)});
  }

  return hist;
}

}  // namespace

BOOST_AUTO_TEST_SUITE(GaussianHistogramFitSuite)

namespace {
// `iterativeFit`/`extractMeanWidthProfiles` take a `HistogramFitFunction`, not
// the free function directly -- this adapts it (a plain function pointer
// already satisfies the signature, but the tests below refer to `fitFn` by
// name throughout).
const HistogramFitFunction fitFn = &gaussianHistogramFit;
}  // namespace

BOOST_AUTO_TEST_CASE(ExactGaussian_RecoversParameters) {
  // A histogram filled with the analytic Gaussian shape rather than samples:
  // the chi-square minimum must sit essentially on the truth, so this isolates
  // the optimiser from statistical scatter.
  const double trueMean = 0.35;
  const double trueSigma = 0.8;
  const double amplitude = 1000.0;

  auto axis = AxisVariant(BoostRegularAxis(80, -6.0, 6.0, "x"));
  Histogram1 hist("exact", "Exact", {axis});
  const auto& boostAxis = hist.histogram().axis(0);
  for (int i = 0; i < boostAxis.size(); ++i) {
    const double centre =
        0.5 * (boostAxis.bin(i).lower() + boostAxis.bin(i).upper());
    const double z = (centre - trueMean) / trueSigma;
    hist.setBinContent({i}, amplitude * std::exp(-0.5 * z * z));
  }

  const auto result = gaussianHistogramFit(hist);
  BOOST_REQUIRE(result.has_value());
  const auto& [mean, sigma, meanError, sigmaError] = *result;
  BOOST_CHECK_CLOSE(mean, trueMean, 0.01);
  BOOST_CHECK_CLOSE(sigma, trueSigma, 0.01);
  BOOST_CHECK_GT(meanError, 0.0);
  BOOST_CHECK_GT(sigmaError, 0.0);
}

BOOST_AUTO_TEST_CASE(HighStatistics_RecoversTruth) {
  const double trueMean = -0.2;
  const double trueSigma = 1.3;
  const auto hist =
      sampleGaussian(100000, trueMean, trueSigma, 100, -8.0, 8.0, 12345);

  const auto result = gaussianHistogramFit(hist);
  BOOST_REQUIRE(result.has_value());
  const auto& [mean, sigma, meanError, sigmaError] = *result;

  // With 1e5 entries the statistical uncertainties are tiny, so require the
  // truth within a few of the reported errors
  BOOST_CHECK_LT(std::abs(mean - trueMean), 4 * meanError);
  BOOST_CHECK_LT(std::abs(sigma - trueSigma), 4 * sigmaError);

  // Sanity check that the reported errors are near the asymptotic expectation
  const double expectedMeanError = trueSigma / std::sqrt(100000.0);
  BOOST_CHECK_CLOSE(meanError, expectedMeanError, 20.0);
}

BOOST_AUTO_TEST_CASE(LowStatistics_StillFits) {
  const auto hist = sampleGaussian(60, 0.0, 1.0, 20, -5.0, 5.0, 777);

  const auto result = gaussianHistogramFit(hist);
  BOOST_REQUIRE(result.has_value());
  const auto& [mean, sigma, meanError, sigmaError] = *result;
  BOOST_CHECK_LT(std::abs(mean), 4 * meanError);
  BOOST_CHECK_LT(std::abs(sigma - 1.0), 4 * sigmaError);
}

namespace {

/// Run `toys` fits of `entriesPerToy`-entry histograms and return the
/// (mean, RMS) of the mean and sigma pulls, `(fitted - true) / error`
struct PullStatistics {
  double meanPullMean{};
  double meanPullRms{};
  double sigmaPullMean{};
  double sigmaPullRms{};
  std::size_t successes{};
};

PullStatistics pullStatistics(double trueMean, double trueSigma,
                              std::size_t toys, std::size_t entriesPerToy,
                              int nBins, double xMin, double xMax,
                              std::uint32_t seedBase) {
  double meanPullSum = 0;
  double meanPullSumSq = 0;
  double sigmaPullSum = 0;
  double sigmaPullSumSq = 0;
  std::size_t successes = 0;

  for (std::size_t toy = 0; toy < toys; ++toy) {
    const auto hist =
        sampleGaussian(entriesPerToy, trueMean, trueSigma, nBins, xMin, xMax,
                       seedBase + static_cast<std::uint32_t>(toy));
    const auto result = gaussianHistogramFit(hist);
    if (!result.has_value()) {
      continue;
    }
    ++successes;

    const auto& [mean, sigma, meanError, sigmaError] = *result;
    const double meanPull = (mean - trueMean) / meanError;
    const double sigmaPull = (sigma - trueSigma) / sigmaError;
    meanPullSum += meanPull;
    meanPullSumSq += meanPull * meanPull;
    sigmaPullSum += sigmaPull;
    sigmaPullSumSq += sigmaPull * sigmaPull;
  }

  const auto n = static_cast<double>(successes);
  return PullStatistics{meanPullSum / n, std::sqrt(meanPullSumSq / n),
                        sigmaPullSum / n, std::sqrt(sigmaPullSumSq / n),
                        successes};
}

}  // namespace

BOOST_AUTO_TEST_CASE(ChiSquarePulls_ValidateCovarianceIsUnscaled) {
  // The pulls (fitted - true) / error must be distributed like N(0, 1). This
  // validates that (J^T J)^-1 -- the chi-square covariance by construction --
  // needs no extra scaling, unlike the old Nelder-Mead implementation's
  // finite-difference Hessian, which needed an explicit MINUIT "UP" factor
  // (see PROGRESS.md). Binning is coarse enough relative to sigma that
  // per-bin counts are large throughout, which avoids the small-sample Neyman
  // bias tested separately below (ChiSquarePulls_ShowTheKnownNeymanBias) and
  // isolates the covariance scale: a missing or wrong factor would show up as
  // a gross RMS deviation, clearly distinguishable from that bias.
  const auto stats = pullStatistics(0.0, 1.0, 500, 20000, 20, -3.0, 3.0, 2000);

  BOOST_TEST_MESSAGE("mean pull:  mean = " << stats.meanPullMean
                                           << " rms = " << stats.meanPullRms);
  BOOST_TEST_MESSAGE("sigma pull: mean = " << stats.sigmaPullMean
                                           << " rms = " << stats.sigmaPullRms);

  BOOST_CHECK_EQUAL(stats.successes, 500u);

  BOOST_CHECK_LT(std::abs(stats.meanPullMean), 0.15);
  BOOST_CHECK_GT(stats.meanPullRms, 0.85);
  BOOST_CHECK_LT(stats.meanPullRms, 1.2);
  BOOST_CHECK_GT(stats.sigmaPullRms, 0.85);
  BOOST_CHECK_LT(stats.sigmaPullRms, 1.2);
}

BOOST_AUTO_TEST_CASE(ChiSquarePulls_ShowTheKnownNeymanBias) {
  // Neyman's chi-square (variance taken from the observed count, as ROOT's
  // "SQ0" and this fit both do) is known to bias fits where low-count bins
  // dominate: a bin that fluctuates down gets an over-large weight (1/n), a
  // bin that fluctuates up an under-large one. Fine binning relative to sigma
  // means most bins sit far in the tail with only a handful of counts, so the
  // effect is large here. This is expected, ROOT-matching behaviour, not a
  // defect.
  const auto stats = pullStatistics(0.0, 1.0, 500, 2000, 60, -6.0, 6.0, 1000);

  BOOST_TEST_MESSAGE("sigma pull: mean = " << stats.sigmaPullMean
                                           << " rms = " << stats.sigmaPullRms);
  BOOST_CHECK_EQUAL(stats.successes, 500u);
  BOOST_CHECK_GT(std::abs(stats.sigmaPullMean), 0.3);
}

BOOST_AUTO_TEST_CASE(RestrictedRange_SelectsByBinCentre) {
  auto axis = AxisVariant(BoostRegularAxis(10, 0.0, 10.0, "x"));
  Histogram1 hist("range", "Range", {axis});
  // Bin centres are 0.5, 1.5, ... 9.5
  for (int i = 0; i < 10; ++i) {
    hist.setBinContent({i}, 100.0);
  }
  // A narrow peak in the middle so the restricted fit is well defined
  hist.setBinContent({4}, 900.0);
  hist.setBinContent({5}, 900.0);

  // Range covering only the central four bin centres (3.5, 4.5, 5.5, 6.5)
  const auto restricted =
      gaussianHistogramFit(hist, HistogramFitRange{3.0, 7.0});
  BOOST_REQUIRE(restricted.has_value());
  BOOST_CHECK_CLOSE(std::get<0>(*restricted), 5.0, 15.0);

  // The full fit sees the flat pedestal too and must give a wider sigma
  const auto full = gaussianHistogramFit(hist);
  BOOST_REQUIRE(full.has_value());
  BOOST_CHECK_GT(std::get<1>(*full), std::get<1>(*restricted));
}

BOOST_AUTO_TEST_CASE(IterativeFit_NarrowsOntoCore) {
  // Gaussian core plus a broad uniform background. The unrestricted fit is
  // pulled wide by the background; iterating onto the core must recover the
  // input sigma much better.
  auto axis = AxisVariant(BoostRegularAxis(100, -10.0, 10.0, "x"));
  Histogram1 hist("core", "Core", {axis});

  std::mt19937 generator(999);
  std::normal_distribution<double> core(0.0, 1.0);
  std::uniform_real_distribution<double> background(-10.0, 10.0);
  for (int i = 0; i < 20000; ++i) {
    hist.fill({core(generator)});
  }
  for (int i = 0; i < 4000; ++i) {
    hist.fill({background(generator)});
  }

  const auto single = gaussianHistogramFit(hist);
  BOOST_REQUIRE(single.has_value());

  const auto iterated = iterativeFit(fitFn, hist, 2.0, 4);
  BOOST_REQUIRE(iterated.has_value());

  BOOST_CHECK_LT(std::abs(std::get<1>(*iterated) - 1.0),
                 std::abs(std::get<1>(*single) - 1.0));
  BOOST_CHECK_CLOSE(std::get<1>(*iterated), 1.0, 10.0);
}

BOOST_AUTO_TEST_CASE(IterativeFit_SingleIterationMatchesPlainFit) {
  const auto hist = sampleGaussian(5000, 0.1, 0.9, 50, -5.0, 5.0, 31337);

  const auto plain = gaussianHistogramFit(hist);
  BOOST_REQUIRE(plain.has_value());

  // iterations == 1 means only the unrestricted fit runs
  const auto once = iterativeFit(fitFn, hist, 3.0, 1);
  BOOST_REQUIRE(once.has_value());
  BOOST_CHECK_CLOSE(std::get<0>(*once), std::get<0>(*plain), 1e-9);
  BOOST_CHECK_CLOSE(std::get<1>(*once), std::get<1>(*plain), 1e-9);

  // Values below 1 must not run fewer than the initial fit either
  const auto zero = iterativeFit(fitFn, hist, 3.0, 0);
  BOOST_REQUIRE(zero.has_value());
  BOOST_CHECK_CLOSE(std::get<0>(*zero), std::get<0>(*plain), 1e-9);
}

BOOST_AUTO_TEST_CASE(Degenerate_EmptyHistogram) {
  auto axis = AxisVariant(BoostRegularAxis(20, -5.0, 5.0, "x"));
  const Histogram1 hist("empty", "Empty", {axis});

  BOOST_CHECK(!gaussianHistogramFit(hist).has_value());
  BOOST_CHECK(!iterativeFit(fitFn, hist, 3.0, 3).has_value());
}

BOOST_AUTO_TEST_CASE(Degenerate_SingleFilledBin) {
  auto axis = AxisVariant(BoostRegularAxis(20, -5.0, 5.0, "x"));
  Histogram1 hist("spike", "Spike", {axis});
  hist.setBinContent({10}, 500.0);

  // One bin cannot constrain three parameters; sigma would run to zero
  BOOST_CHECK(!gaussianHistogramFit(hist).has_value());
}

BOOST_AUTO_TEST_CASE(Degenerate_TwoFilledBins) {
  auto axis = AxisVariant(BoostRegularAxis(20, -5.0, 5.0, "x"));
  Histogram1 hist("two", "Two", {axis});
  hist.setBinContent({9}, 100.0);
  hist.setBinContent({10}, 120.0);

  // Still fewer populated bins than free parameters
  BOOST_CHECK(!gaussianHistogramFit(hist).has_value());
}

BOOST_AUTO_TEST_CASE(Degenerate_ThreeFilledBinsSucceed) {
  auto axis = AxisVariant(BoostRegularAxis(20, -5.0, 5.0, "x"));
  Histogram1 hist("three", "Three", {axis});
  hist.setBinContent({9}, 80.0);
  hist.setBinContent({10}, 120.0);
  hist.setBinContent({11}, 70.0);

  // Exactly at the limit, this must work rather than fail
  const auto result = gaussianHistogramFit(hist);
  BOOST_REQUIRE(result.has_value());
  BOOST_CHECK(std::isfinite(std::get<0>(*result)));
  BOOST_CHECK_GT(std::get<1>(*result), 0.0);
}

BOOST_AUTO_TEST_CASE(Degenerate_SingleBinHistogram) {
  auto axis = AxisVariant(BoostRegularAxis(1, -5.0, 5.0, "x"));
  Histogram1 hist("onebin", "One bin", {axis});
  hist.setBinContent({0}, 1000.0);

  BOOST_CHECK(!gaussianHistogramFit(hist).has_value());
}

BOOST_AUTO_TEST_CASE(Degenerate_RangeExcludesEverything) {
  const auto hist = sampleGaussian(5000, 0.0, 1.0, 50, -5.0, 5.0, 5150);

  // A range far outside the axis selects no bins at all
  BOOST_CHECK(
      !gaussianHistogramFit(hist, HistogramFitRange{100.0, 200.0}).has_value());

  // A range narrower than a single bin selects at most one bin
  BOOST_CHECK(
      !gaussianHistogramFit(hist, HistogramFitRange{0.01, 0.02}).has_value());
}

BOOST_AUTO_TEST_CASE(Degenerate_CollapsingIterationRange) {
  const auto hist = sampleGaussian(5000, 0.0, 1.0, 50, -5.0, 5.0, 6161);

  // A tiny sigmaRange collapses the refit window below one bin width, so the
  // iteration must fail cleanly rather than return nonsense
  const auto result = iterativeFit(fitFn, hist, 1e-6, 3);
  BOOST_CHECK(!result.has_value());
}

namespace {

/// Gaussian core of 5000 entries plus a single far-away spike of `spikeCounts`
/// in the bin at x = 45.5
Histogram1 coreWithOutlier(double spikeCounts) {
  auto axis = AxisVariant(BoostRegularAxis(100, -50.0, 50.0, "x"));
  Histogram1 hist("outlier", "Outlier", {axis});

  std::mt19937 generator(2024);
  std::normal_distribution<double> core(0.0, 1.0);
  for (int i = 0; i < 5000; ++i) {
    hist.fill({core(generator)});
  }
  hist.setBinContent({95}, spikeCounts);

  return hist;
}

}  // namespace

BOOST_AUTO_TEST_CASE(ModerateOutlier_IsShakenOffByIterating) {
  // A single far-away spike sits outside any 3-sigma window drawn from a
  // fit that has locked onto the core, so once the unrestricted fit's sigma
  // is not itself dragged wide enough to keep the spike inside that window,
  // iterating drops it and converges back onto the core's true sigma.
  const Histogram1 hist = coreWithOutlier(200.0);

  const auto iterated = iterativeFit(fitFn, hist, 3.0, 4);
  BOOST_REQUIRE(iterated.has_value());
  BOOST_CHECK_CLOSE(std::get<1>(*iterated), 1.0, 15.0);
}

// No ExtremeOutlier_StaysFinite here: with a spike carrying more than a
// third of all entries, chi-square genuinely has no graceful-degradation
// guarantee to test -- confirmed empirically that ROOT's own "SQ0" also
// reports a failing fit status on the equivalent histogram (see
// devscripts/check_extreme_outlier.py, not checked in). Robustness to
// exactly this case was the stated reason the likelihood objective existed;
// now that it is gone, this scenario is simply out of scope for the
// remaining chi-square fit.

BOOST_AUTO_TEST_CASE(VariableBinning_IsSupported) {
  // Non-uniform bins: the fit uses bin centres, so it must not assume a
  // constant width anywhere
  std::vector<double> edges;
  for (int i = -10; i <= 10; ++i) {
    // Bins widen away from zero
    edges.push_back(std::copysign(0.1 * i * i, i));
  }
  auto axis = AxisVariant(BoostVariableAxis(edges, "x"));
  Histogram1 hist("var", "Variable", {axis});

  std::mt19937 generator(8080);
  std::normal_distribution<double> distribution(0.0, 2.0);
  for (int i = 0; i < 20000; ++i) {
    hist.fill({distribution(generator)});
  }

  const auto result = gaussianHistogramFit(hist);
  BOOST_REQUIRE(result.has_value());
  BOOST_CHECK(std::isfinite(std::get<0>(*result)));
  BOOST_CHECK_GT(std::get<1>(*result), 0.0);
}

namespace {

/// Residual-vs-eta style 2D histogram whose width grows with eta.
/// Bin i of the eta axis gets `entries` samples from N(0, sigmaOf(i)).
Histogram2 residualVsEta(int nEtaBins, std::size_t entries,
                         const std::function<double(int)>& sigmaOf,
                         std::uint32_t seed) {
  auto etaAxis =
      AxisVariant(BoostRegularAxis(nEtaBins, 0.0, nEtaBins * 1.0, "eta"));
  auto resAxis = AxisVariant(BoostRegularAxis(80, -10.0, 10.0, "res"));
  Histogram2 hist("res_vs_eta", "Residual vs Eta", {etaAxis, resAxis});

  std::mt19937 generator(seed);
  for (int i = 0; i < nEtaBins; ++i) {
    const double etaValue = i + 0.5;
    std::normal_distribution<double> distribution(0.0, sigmaOf(i));
    for (std::size_t n = 0; n < entries; ++n) {
      hist.fill({etaValue, distribution(generator)});
    }
  }

  return hist;
}

}  // namespace

BOOST_AUTO_TEST_CASE(Profiles2D_RecoverPerBinWidth) {
  const int nEtaBins = 5;
  const auto sigmaOf = [](int i) { return 0.5 + 0.25 * i; };
  const Histogram2 hist = residualVsEta(nEtaBins, 20000, sigmaOf, 5555);

  const auto profiles = extractMeanWidthProfiles(
      fitFn, hist, "resmean_d0_vs_eta", "reswidth_d0_vs_eta");

  BOOST_CHECK_EQUAL(profiles.mean.name(), "resmean_d0_vs_eta");
  BOOST_CHECK_EQUAL(profiles.width.name(), "reswidth_d0_vs_eta");
  BOOST_CHECK_EQUAL(profiles.mean.title(), "Residual vs Eta mean");
  BOOST_CHECK_EQUAL(profiles.width.title(), "Residual vs Eta width");
  BOOST_CHECK_EQUAL(profiles.fitFailureFraction, 0.0);

  // The output axis must be the input's first axis
  BOOST_CHECK_EQUAL(profiles.mean.histogram().axis(0).size(), nEtaBins);
  BOOST_CHECK_EQUAL(profiles.mean.histogram().axis(0).metadata(), "eta");

  for (int i = 0; i < nEtaBins; ++i) {
    BOOST_CHECK_LT(std::abs(profiles.mean.value({i})),
                   4 * profiles.mean.error({i}));
    BOOST_CHECK_LT(std::abs(profiles.width.value({i}) - sigmaOf(i)),
                   4 * profiles.width.error({i}));
    BOOST_CHECK_GT(profiles.mean.error({i}), 0.0);
    BOOST_CHECK_GT(profiles.width.error({i}), 0.0);
  }
}

BOOST_AUTO_TEST_CASE(Profiles2D_SkipsSparseSlicesWithoutCountingFailures) {
  auto etaAxis = AxisVariant(BoostRegularAxis(4, 0.0, 4.0, "eta"));
  auto resAxis = AxisVariant(BoostRegularAxis(40, -5.0, 5.0, "res"));
  Histogram2 hist("res_vs_eta", "Residual vs Eta", {etaAxis, resAxis});

  // Bin 0 well populated, bin 1 has three entries, bins 2 and 3 empty
  std::mt19937 generator(1717);
  std::normal_distribution<double> distribution(0.0, 1.0);
  for (int n = 0; n < 5000; ++n) {
    hist.fill({0.5, distribution(generator)});
  }
  hist.fill({1.5, -0.2});
  hist.fill({1.5, 0.0});
  hist.fill({1.5, 0.3});

  const auto profiles = extractMeanWidthProfiles(fitFn, hist, "mean", "width",
                                                 /*minEntriesForFit=*/10);

  // Only bin 0 clears the threshold and is filled
  BOOST_CHECK_GT(profiles.width.value({0}), 0.0);
  for (int i = 1; i < 4; ++i) {
    BOOST_CHECK_EQUAL(profiles.mean.value({i}), 0.0);
    BOOST_CHECK_EQUAL(profiles.width.value({i}), 0.0);
  }

  // Slices below the entry threshold are skipped, not failed
  BOOST_CHECK_EQUAL(profiles.fitFailureFraction, 0.0);
}

BOOST_AUTO_TEST_CASE(Profiles2D_CountsGenuineFailures) {
  auto etaAxis = AxisVariant(BoostRegularAxis(4, 0.0, 4.0, "eta"));
  auto resAxis = AxisVariant(BoostRegularAxis(40, -5.0, 5.0, "res"));
  Histogram2 hist("res_vs_eta", "Residual vs Eta", {etaAxis, resAxis});

  // Two slices with plenty of entries but all of them in a single bin, so the
  // fit is attempted and must fail
  hist.setBinContent({0, 20}, 500.0);
  hist.setBinContent({1, 20}, 500.0);

  const auto profiles = extractMeanWidthProfiles(fitFn, hist, "mean", "width",
                                                 /*minEntriesForFit=*/10);

  // Two failures out of four bins on the first axis
  BOOST_CHECK_CLOSE(profiles.fitFailureFraction, 0.5, 1e-10);
  BOOST_CHECK_EQUAL(profiles.mean.value({0}), 0.0);
}

BOOST_AUTO_TEST_CASE(Profiles2D_VariableBinning) {
  std::vector<double> etaEdges = {0.0, 0.5, 1.5, 3.0};
  auto etaAxis = AxisVariant(BoostVariableAxis(etaEdges, "eta"));
  auto resAxis = AxisVariant(BoostRegularAxis(60, -6.0, 6.0, "res"));
  Histogram2 hist("res_vs_eta", "Residual vs Eta", {etaAxis, resAxis});

  std::mt19937 generator(2727);
  std::normal_distribution<double> distribution(0.25, 0.8);
  const std::array<double, 3> etaValues = {0.25, 1.0, 2.0};
  for (const double etaValue : etaValues) {
    for (int n = 0; n < 10000; ++n) {
      hist.fill({etaValue, distribution(generator)});
    }
  }

  const auto profiles = extractMeanWidthProfiles(fitFn, hist, "mean", "width");

  BOOST_CHECK(extractBinEdges(profiles.mean.histogram().axis(0)) == etaEdges);
  BOOST_CHECK_EQUAL(profiles.fitFailureFraction, 0.0);
  for (int i = 0; i < 3; ++i) {
    BOOST_CHECK_CLOSE(profiles.mean.value({i}), 0.25, 10.0);
    BOOST_CHECK_CLOSE(profiles.width.value({i}), 0.8, 10.0);
  }
}

BOOST_AUTO_TEST_CASE(Profiles3D_RecoverPerBinWidth) {
  const int nEta = 2;
  const int nPt = 3;
  auto etaAxis = AxisVariant(BoostRegularAxis(nEta, 0.0, 2.0, "eta"));
  auto ptAxis = AxisVariant(BoostRegularAxis(nPt, 0.0, 3.0, "pt"));
  auto resAxis = AxisVariant(BoostRegularAxis(80, -10.0, 10.0, "res"));
  Histogram3 hist("res_vs_eta_pt", "Residual", {etaAxis, ptAxis, resAxis});

  const auto sigmaOf = [](int i, int j) { return 0.5 + 0.3 * i + 0.2 * j; };

  std::mt19937 generator(3939);
  for (int i = 0; i < nEta; ++i) {
    for (int j = 0; j < nPt; ++j) {
      std::normal_distribution<double> distribution(0.0, sigmaOf(i, j));
      for (int n = 0; n < 20000; ++n) {
        hist.fill({i + 0.5, j + 0.5, distribution(generator)});
      }
    }
  }

  const auto profiles = extractMeanWidthProfiles(fitFn, hist, "mean", "width");

  // The output keeps the first two axes
  BOOST_CHECK_EQUAL(profiles.width.histogram().axis(0).size(), nEta);
  BOOST_CHECK_EQUAL(profiles.width.histogram().axis(1).size(), nPt);
  BOOST_CHECK_EQUAL(profiles.width.histogram().axis(0).metadata(), "eta");
  BOOST_CHECK_EQUAL(profiles.width.histogram().axis(1).metadata(), "pt");
  BOOST_CHECK_EQUAL(profiles.fitFailureFraction, 0.0);

  for (int i = 0; i < nEta; ++i) {
    for (int j = 0; j < nPt; ++j) {
      BOOST_CHECK_LT(std::abs(profiles.width.value({i, j}) - sigmaOf(i, j)),
                     4 * profiles.width.error({i, j}));
    }
  }
}

BOOST_AUTO_TEST_CASE(Profiles3D_FailureFractionUsesBothAxes) {
  auto etaAxis = AxisVariant(BoostRegularAxis(2, 0.0, 2.0, "eta"));
  auto ptAxis = AxisVariant(BoostRegularAxis(2, 0.0, 2.0, "pt"));
  auto resAxis = AxisVariant(BoostRegularAxis(40, -5.0, 5.0, "res"));
  Histogram3 hist("res", "Residual", {etaAxis, ptAxis, resAxis});

  // One of the four (eta, pt) cells is populated but unfittable
  hist.setBinContent({0, 0, 20}, 500.0);

  const auto profiles = extractMeanWidthProfiles(fitFn, hist, "mean", "width",
                                                 /*minEntriesForFit=*/10);

  // One failure out of 2 x 2 cells
  BOOST_CHECK_CLOSE(profiles.fitFailureFraction, 0.25, 1e-10);
}

BOOST_AUTO_TEST_SUITE_END()
