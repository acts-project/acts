// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/GaussianHistogramFit.hpp"
#include "Acts/Utilities/Histogram.hpp"
#include "Acts/Utilities/IterativeFit.hpp"

#include <array>
#include <cmath>
#include <functional>
#include <random>
#include <vector>

using namespace Acts;
using namespace Acts::Experimental;

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
// The default fitter, i.e. chi-square -- this is what ROOT's "SQ0" produces
// and what most tests below exercise.
const GaussianHistogramFit fit;
// The likelihood fitter, used where a test specifically exercises the
// likelihood's handling of empty bins or outliers (see PROGRESS notes: this
// was the original motivation for offering it at all).
const GaussianHistogramFit likelihoodFit{GaussianHistogramFit::Config{
    GaussianHistogramFit::Objective::PoissonLikelihood}};
}  // namespace

BOOST_AUTO_TEST_CASE(ExactGaussian_RecoversParameters) {
  // A histogram filled with the analytic Gaussian shape rather than samples:
  // the likelihood maximum must sit essentially on the truth, so this isolates
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

  const auto result = fit.fit(hist);
  BOOST_REQUIRE(result.ok());
  BOOST_CHECK_CLOSE(result->mean, trueMean, 0.01);
  BOOST_CHECK_CLOSE(result->sigma, trueSigma, 0.01);
  BOOST_CHECK_GT(result->meanError, 0.0);
  BOOST_CHECK_GT(result->sigmaError, 0.0);
}

BOOST_AUTO_TEST_CASE(HighStatistics_RecoversTruth) {
  const double trueMean = -0.2;
  const double trueSigma = 1.3;
  const auto hist =
      sampleGaussian(100000, trueMean, trueSigma, 100, -8.0, 8.0, 12345);

  const auto result = fit.fit(hist);
  BOOST_REQUIRE(result.ok());

  // With 1e5 entries the statistical uncertainties are tiny, so require the
  // truth within a few of the reported errors
  BOOST_CHECK_LT(std::abs(result->mean - trueMean), 4 * result->meanError);
  BOOST_CHECK_LT(std::abs(result->sigma - trueSigma), 4 * result->sigmaError);

  // Sanity check that the reported errors are near the asymptotic expectation
  const double expectedMeanError = trueSigma / std::sqrt(100000.0);
  BOOST_CHECK_CLOSE(result->meanError, expectedMeanError, 20.0);
}

BOOST_AUTO_TEST_CASE(LowStatistics_StillFits) {
  const auto hist = sampleGaussian(60, 0.0, 1.0, 20, -5.0, 5.0, 777);

  const auto result = fit.fit(hist);
  BOOST_REQUIRE(result.ok());
  BOOST_CHECK_LT(std::abs(result->mean), 4 * result->meanError);
  BOOST_CHECK_LT(std::abs(result->sigma - 1.0), 4 * result->sigmaError);
}

BOOST_AUTO_TEST_CASE(SparseHistogram_EmptyBinsAreInformative) {
  // Fine binning with few entries leaves most bins empty. The likelihood
  // treats those bins as informative and lands close to the truth.
  const auto hist = sampleGaussian(80, 0.5, 0.4, 200, -5.0, 5.0, 4242);

  const auto result = likelihoodFit.fit(hist);
  BOOST_REQUIRE(result.ok());
  BOOST_CHECK_LT(std::abs(result->mean - 0.5), 4 * result->meanError);
  BOOST_CHECK_LT(std::abs(result->sigma - 0.4), 4 * result->sigmaError);
}

BOOST_AUTO_TEST_CASE(ChiSquareVsLikelihood_DifferOnSparseHistogram) {
  // Same histogram as above. Chi-square drops every empty bin from the sum,
  // so with this few entries and this fine a binning it is dominated by the
  // handful of populated bins that happen to have low counts (the Neyman
  // bias -- see ChiSquarePulls_ShowTheKnownNeymanBias); the likelihood, for
  // which empty bins are informative, recovers the truth substantially
  // better. This is the documented rationale for offering both objectives.
  const auto hist = sampleGaussian(80, 0.5, 0.4, 200, -5.0, 5.0, 4242);

  const auto chiSquareResult = fit.fit(hist);
  const auto likelihoodResult = likelihoodFit.fit(hist);

  BOOST_REQUIRE(chiSquareResult.ok());
  BOOST_REQUIRE(likelihoodResult.ok());

  const double chiSquareSigmaError = std::abs(chiSquareResult->sigma - 0.4);
  const double likelihoodSigmaError = std::abs(likelihoodResult->sigma - 0.4);
  BOOST_TEST_MESSAGE("chi2 sigma error = " << chiSquareSigmaError
                                           << ", likelihood sigma error = "
                                           << likelihoodSigmaError);
  BOOST_CHECK_GT(chiSquareSigmaError, 3 * likelihoodSigmaError);
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

PullStatistics pullStatistics(const GaussianHistogramFit& fitter,
                              double trueMean, double trueSigma,
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
    const auto result = fitter.fit(hist);
    if (!result.ok()) {
      continue;
    }
    ++successes;

    const double meanPull = (result->mean - trueMean) / result->meanError;
    const double sigmaPull = (result->sigma - trueSigma) / result->sigmaError;
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

BOOST_AUTO_TEST_CASE(LikelihoodPulls_AreStandardNormal) {
  // The pulls (fitted - true) / error must be distributed like N(0, 1). This
  // validates the curvature-based uncertainties independently of ROOT and
  // catches a systematically mis-scaled Hessian. The likelihood is unbiased
  // here even with many low-count tail bins, which is why it is used for this
  // check rather than chi-square (see ChiSquarePulls_ShowTheKnownNeymanBias).
  const auto stats =
      pullStatistics(likelihoodFit, 0.0, 1.0, 500, 2000, 60, -6.0, 6.0, 1000);

  BOOST_TEST_MESSAGE("mean pull:  mean = " << stats.meanPullMean
                                           << " rms = " << stats.meanPullRms);
  BOOST_TEST_MESSAGE("sigma pull: mean = " << stats.sigmaPullMean
                                           << " rms = " << stats.sigmaPullRms);

  // Every toy is comfortably populated, so none of them may fail
  BOOST_CHECK_EQUAL(stats.successes, 500u);

  // Bounds are loose enough that a correct implementation cannot fail by
  // chance, but tight enough to catch an error scale off by more than ~20%
  BOOST_CHECK_LT(std::abs(stats.meanPullMean), 0.15);
  BOOST_CHECK_GT(stats.meanPullRms, 0.85);
  BOOST_CHECK_LT(stats.meanPullRms, 1.2);
  BOOST_CHECK_LT(std::abs(stats.sigmaPullMean), 0.15);
  BOOST_CHECK_GT(stats.sigmaPullRms, 0.85);
  BOOST_CHECK_LT(stats.sigmaPullRms, 1.2);
}

BOOST_AUTO_TEST_CASE(ChiSquarePulls_ValidateTheFactorTwoCovarianceScaling) {
  // Same check as LikelihoodPulls_AreStandardNormal, for chi-square, with
  // binning coarse enough relative to sigma that per-bin counts are large
  // throughout. This does not eliminate the mean-pull bias chi-square is
  // known to have (see ChiSquarePulls_ShowTheKnownNeymanBias) -- it is a
  // small-sample effect in the *estimator* itself, not something more
  // statistics removes by shrinking the error alone -- so only the pull RMS
  // is checked here. A missing or wrong UP factor in the covariance would
  // show up as a gross RMS deviation (e.g. ~1.4x for a missing factor of 2),
  // clearly distinguishable from the residual bias.
  const auto stats =
      pullStatistics(fit, 0.0, 1.0, 500, 20000, 20, -3.0, 3.0, 2000);

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
  // "SQ0" and this objective both do) is known to bias fits low-count bins
  // dominate: a bin that fluctuates down gets an over-large weight (1/n), a
  // bin that fluctuates up an under-large one. Fine binning relative to sigma
  // means most bins sit far in the tail with only a handful of counts, so the
  // effect is large here. This is expected, ROOT-matching behaviour, not a
  // defect -- it is exactly why the likelihood objective exists.
  const auto stats =
      pullStatistics(fit, 0.0, 1.0, 500, 2000, 60, -6.0, 6.0, 1000);

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
  const auto restricted = fit.fit(hist, 3.0, 7.0);
  BOOST_REQUIRE(restricted.ok());
  BOOST_CHECK_CLOSE(restricted->mean, 5.0, 15.0);

  // The full fit sees the flat pedestal too and must give a wider sigma
  const auto full = fit.fit(hist);
  BOOST_REQUIRE(full.ok());
  BOOST_CHECK_GT(full->sigma, restricted->sigma);
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

  const auto single = fit.fit(hist);
  BOOST_REQUIRE(single.ok());

  const auto iterated = iterativeFit(fit, hist, 2.0, 4);
  BOOST_REQUIRE(iterated.ok());

  BOOST_CHECK_LT(std::abs(iterated->sigma - 1.0),
                 std::abs(single->sigma - 1.0));
  BOOST_CHECK_CLOSE(iterated->sigma, 1.0, 10.0);
}

BOOST_AUTO_TEST_CASE(IterativeFit_SingleIterationMatchesPlainFit) {
  const auto hist = sampleGaussian(5000, 0.1, 0.9, 50, -5.0, 5.0, 31337);

  const auto plain = fit.fit(hist);
  BOOST_REQUIRE(plain.ok());

  // iterations == 1 means only the unrestricted fit runs
  const auto once = iterativeFit(fit, hist, 3.0, 1);
  BOOST_REQUIRE(once.ok());
  BOOST_CHECK_CLOSE(once->mean, plain->mean, 1e-9);
  BOOST_CHECK_CLOSE(once->sigma, plain->sigma, 1e-9);

  // Values below 1 must not run fewer than the initial fit either
  const auto zero = iterativeFit(fit, hist, 3.0, 0);
  BOOST_REQUIRE(zero.ok());
  BOOST_CHECK_CLOSE(zero->mean, plain->mean, 1e-9);
}

BOOST_AUTO_TEST_CASE(Degenerate_EmptyHistogram) {
  auto axis = AxisVariant(BoostRegularAxis(20, -5.0, 5.0, "x"));
  const Histogram1 hist("empty", "Empty", {axis});

  BOOST_CHECK(!fit.fit(hist).ok());
  BOOST_CHECK(!iterativeFit(fit, hist, 3.0, 3).ok());
}

BOOST_AUTO_TEST_CASE(Degenerate_SingleFilledBin) {
  auto axis = AxisVariant(BoostRegularAxis(20, -5.0, 5.0, "x"));
  Histogram1 hist("spike", "Spike", {axis});
  hist.setBinContent({10}, 500.0);

  // One bin cannot constrain three parameters; sigma would run to zero
  BOOST_CHECK(!fit.fit(hist).ok());
}

BOOST_AUTO_TEST_CASE(Degenerate_TwoFilledBins) {
  auto axis = AxisVariant(BoostRegularAxis(20, -5.0, 5.0, "x"));
  Histogram1 hist("two", "Two", {axis});
  hist.setBinContent({9}, 100.0);
  hist.setBinContent({10}, 120.0);

  // Still fewer populated bins than free parameters
  BOOST_CHECK(!fit.fit(hist).ok());
}

BOOST_AUTO_TEST_CASE(Degenerate_ThreeFilledBinsSucceed) {
  auto axis = AxisVariant(BoostRegularAxis(20, -5.0, 5.0, "x"));
  Histogram1 hist("three", "Three", {axis});
  hist.setBinContent({9}, 80.0);
  hist.setBinContent({10}, 120.0);
  hist.setBinContent({11}, 70.0);

  // Exactly at the limit, this must work rather than fail
  const auto result = fit.fit(hist);
  BOOST_REQUIRE(result.ok());
  BOOST_CHECK(std::isfinite(result->mean));
  BOOST_CHECK_GT(result->sigma, 0.0);
}

BOOST_AUTO_TEST_CASE(Degenerate_SingleBinHistogram) {
  auto axis = AxisVariant(BoostRegularAxis(1, -5.0, 5.0, "x"));
  Histogram1 hist("onebin", "One bin", {axis});
  hist.setBinContent({0}, 1000.0);

  BOOST_CHECK(!fit.fit(hist).ok());
}

BOOST_AUTO_TEST_CASE(Degenerate_RangeExcludesEverything) {
  const auto hist = sampleGaussian(5000, 0.0, 1.0, 50, -5.0, 5.0, 5150);

  // A range far outside the axis selects no bins at all
  BOOST_CHECK(!fit.fit(hist, 100.0, 200.0).ok());

  // A range narrower than a single bin selects at most one bin
  BOOST_CHECK(!fit.fit(hist, 0.01, 0.02).ok());
}

BOOST_AUTO_TEST_CASE(Degenerate_CollapsingIterationRange) {
  const auto hist = sampleGaussian(5000, 0.0, 1.0, 50, -5.0, 5.0, 6161);

  // A tiny sigmaRange collapses the refit window below one bin width, so the
  // iteration must fail cleanly rather than return nonsense
  const auto result = iterativeFit(fit, hist, 1e-6, 3);
  BOOST_CHECK(!result.ok());
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
  // With a spike carrying 200 of 5200 entries the unrestricted likelihood fit
  // is pulled wide (sigma ~ 9), but that puts the spike outside the 3-sigma
  // window, so the next iteration drops it and converges back onto the core.
  // Uses the likelihood explicitly: chi-square weights a high-count bin like
  // the spike by 1/n, so it is far less sensitive to this kind of outlier in
  // the first place and the scenario does not exercise the same thing there.
  const Histogram1 hist = coreWithOutlier(200.0);

  const auto single = likelihoodFit.fit(hist);
  BOOST_REQUIRE(single.ok());
  BOOST_CHECK_GT(single->sigma, 4.0);

  const auto iterated = iterativeFit(likelihoodFit, hist, 3.0, 4);
  BOOST_REQUIRE(iterated.ok());
  BOOST_CHECK_CLOSE(iterated->sigma, 1.0, 15.0);
}

BOOST_AUTO_TEST_CASE(ExtremeOutlier_StaysFinite) {
  // A spike carrying 3000 of 8000 entries genuinely dominates the likelihood:
  // the wide solution has a far lower NLL than the core one, so sigma ~ 33 is
  // the correct maximum-likelihood answer and iterating cannot escape it (the
  // spike stays inside the 3-sigma window). What we require here is only that
  // the fit degrades gracefully - finite, positive sigma, finite errors -
  // rather than diverging or returning NaN.
  const Histogram1 hist = coreWithOutlier(3000.0);

  const auto checkDegradesGracefully = [](const auto& result) {
    BOOST_REQUIRE(result.ok());
    BOOST_CHECK(std::isfinite(result->mean));
    BOOST_CHECK(std::isfinite(result->sigma));
    BOOST_CHECK_GT(result->sigma, 0.0);
    BOOST_CHECK(std::isfinite(result->meanError));
    BOOST_CHECK_GT(result->meanError, 0.0);
    BOOST_CHECK(std::isfinite(result->sigmaError));
    BOOST_CHECK_GT(result->sigmaError, 0.0);
  };
  checkDegradesGracefully(likelihoodFit.fit(hist));
  checkDegradesGracefully(iterativeFit(likelihoodFit, hist, 3.0, 4));
}

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

  const auto result = fit.fit(hist);
  BOOST_REQUIRE(result.ok());
  BOOST_CHECK(std::isfinite(result->mean));
  BOOST_CHECK_GT(result->sigma, 0.0);
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

  const auto profiles = extractMeanWidthProfiles(fit, hist, "resmean_d0_vs_eta",
                                                 "reswidth_d0_vs_eta");

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

  const auto profiles = extractMeanWidthProfiles(fit, hist, "mean", "width",
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

  const auto profiles = extractMeanWidthProfiles(fit, hist, "mean", "width",
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

  const auto profiles = extractMeanWidthProfiles(fit, hist, "mean", "width");

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

  const auto profiles = extractMeanWidthProfiles(fit, hist, "mean", "width");

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

  const auto profiles = extractMeanWidthProfiles(fit, hist, "mean", "width",
                                                 /*minEntriesForFit=*/10);

  // One failure out of 2 x 2 cells
  BOOST_CHECK_CLOSE(profiles.fitFailureFraction, 0.25, 1e-10);
}

BOOST_AUTO_TEST_SUITE_END()
