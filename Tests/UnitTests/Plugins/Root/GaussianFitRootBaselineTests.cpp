// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Validates the ROOT-free Gaussian fit in Core against ROOT itself.
//
// ROOT is the reference for this algorithm, so the point of these tests is not
// that Core's fit is sane in isolation (GaussianFitTests.cpp covers that
// without ROOT) but that it agrees with what ROOT would have produced on the
// very same histogram.
//
// The reference is `TH1::Fit("gaus", "LSQ0")`, i.e. ROOT's *likelihood* fit.
// Note that `ActsPlugins::RootHistogramFit` still uses "SQ0", the chi-square
// default; switching it to "LSQ0" is a separate change. Until that lands,
// comparing through that entry point would compare two different objectives,
// so the reference fit is driven directly here.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/GaussianFit.hpp"
#include "Acts/Utilities/GaussianFitError.hpp"
#include "Acts/Utilities/Histogram.hpp"
#include "ActsPlugins/Root/HistogramConverter.hpp"
#include "ActsPlugins/Root/RootHistogramFit.hpp"

#include <cmath>
#include <optional>
#include <random>
#include <string>
#include <vector>

#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>

using namespace Acts;
using namespace Acts::Experimental;
using ActsPlugins::toRoot;

namespace {

const GaussianHistogramFit fit;

/// A named histogram to fit with both implementations
struct Scenario {
  std::string name;
  Histogram1 hist;
};

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

/// The scenarios both implementations are run over.
///
/// Deliberately includes the cases that motivated moving to a likelihood fit:
/// sparse histograms with many empty bins, and peaks contaminated by
/// background or outliers.
std::vector<Scenario> scenarios() {
  std::vector<Scenario> all;

  all.push_back({"high_stats",
                 sampled("high_stats", 100000, 0.0, 1.0, 100, -8.0, 8.0, 101)});
  all.push_back({"medium_stats",
                 sampled("medium_stats", 2000, 0.3, 0.7, 60, -5.0, 5.0, 102)});
  all.push_back(
      {"low_stats", sampled("low_stats", 120, -0.4, 1.1, 30, -6.0, 6.0, 103)});
  all.push_back(
      {"narrow_peak_coarse_bins", sampled("narrow_peak_coarse_bins", 20000, 0.0,
                                          0.35, 20, -5.0, 5.0, 104)});
  all.push_back(
      {"wide_peak_fine_bins",
       sampled("wide_peak_fine_bins", 20000, 0.0, 2.5, 200, -10.0, 10.0, 105)});
  all.push_back(
      {"sparse_many_empty_bins",
       sampled("sparse_many_empty_bins", 150, 0.5, 0.4, 200, -5.0, 5.0, 106)});
  all.push_back({"offset_peak",
                 sampled("offset_peak", 20000, 3.2, 0.6, 80, -8.0, 8.0, 107)});

  // Peak close to the axis edge, so one tail is truncated
  all.push_back(
      {"edge_peak", sampled("edge_peak", 20000, 4.2, 0.8, 60, -5.0, 5.0, 108)});

  // Gaussian core on a uniform background
  {
    Histogram1 hist = makeHistogram("core_plus_background", 100, -10.0, 10.0);
    std::mt19937 generator(109);
    std::normal_distribution<double> core(0.0, 1.0);
    std::uniform_real_distribution<double> flat(-10.0, 10.0);
    for (int i = 0; i < 20000; ++i) {
      hist.fill({core(generator)});
    }
    for (int i = 0; i < 3000; ++i) {
      hist.fill({flat(generator)});
    }
    all.push_back({"core_plus_background", hist});
  }

  // Gaussian core with a handful of isolated far-out entries
  {
    Histogram1 hist = makeHistogram("core_plus_outliers", 100, -20.0, 20.0);
    std::mt19937 generator(110);
    std::normal_distribution<double> core(0.2, 1.4);
    for (int i = 0; i < 20000; ++i) {
      hist.fill({core(generator)});
    }
    hist.fill({15.5});
    hist.fill({16.2});
    hist.fill({-14.8});
    all.push_back({"core_plus_outliers", hist});
  }

  // Variable bin widths
  {
    std::vector<double> edges;
    for (int i = -20; i <= 20; ++i) {
      edges.push_back(std::copysign(0.02 * i * i, i));
    }
    auto axis = AxisVariant(BoostVariableAxis(edges, "x"));
    Histogram1 hist("variable_bins", "variable_bins", {axis});
    std::mt19937 generator(111);
    std::normal_distribution<double> distribution(0.0, 2.0);
    for (int i = 0; i < 20000; ++i) {
      hist.fill({distribution(generator)});
    }
    all.push_back({"variable_bins", hist});
  }

  return all;
}

/// ROOT reference: a single likelihood fit over the full range
std::optional<GaussianFitResult> rootFit(TH1& hist) {
  TFitResultPtr result = hist.Fit("gaus", "LSQ0", nullptr);
  if (result.Get() == nullptr || result->Status() % 1000 != 0) {
    return std::nullopt;
  }

  return GaussianFitResult{result->Parameter(1), result->Parameter(2),
                           result->ParError(1), result->ParError(2)};
}

/// ROOT reference: the same iterative narrowing that `iterativeGaussianFit`
/// performs, driven through ROOT
std::optional<GaussianFitResult> rootIterativeFit(TH1& hist, double sigmaRange,
                                                  int iterations) {
  TFitResultPtr result = hist.Fit("gaus", "LSQ0", nullptr);
  if (result.Get() == nullptr || result->Status() % 1000 != 0) {
    return std::nullopt;
  }

  double mean = result->Parameter(1);
  double sigma = result->Parameter(2);

  for (int i = 0; i < iterations - 1; ++i) {
    const double xMin = mean - sigmaRange * sigma;
    const double xMax = mean + sigmaRange * sigma;

    result = hist.Fit("gaus", "LRSQ0", nullptr, xMin, xMax);
    if (result.Get() == nullptr || result->Status() % 1000 != 0) {
      return std::nullopt;
    }

    mean = result->Parameter(1);
    sigma = result->Parameter(2);
  }

  return GaussianFitResult{mean, sigma, result->ParError(1),
                           result->ParError(2)};
}

/// Require agreement to `relativeTolerance`, or to a fraction of the fitted
/// uncertainty, whichever is looser. The second leg matters for the low
/// statistics scenarios, where the two minimisers legitimately stop at
/// slightly different points on a shallow likelihood.
void checkAgrees(const std::string& what, double ours, double reference,
                 double error, double relativeTolerance,
                 double errorTolerance) {
  const double allowed =
      std::max(relativeTolerance * std::abs(reference), errorTolerance * error);
  const double difference = std::abs(ours - reference);

  BOOST_TEST_CONTEXT(what << ": ours = " << ours << ", ROOT = " << reference
                          << ", diff = " << difference
                          << ", allowed = " << allowed) {
    BOOST_CHECK_LE(difference, allowed);
  }
}

/// Require plain relative agreement. Used for the uncertainties, where there
/// is no second scale to fall back on.
void checkRelative(const std::string& what, double ours, double reference,
                   double relativeTolerance) {
  const double allowed = relativeTolerance * std::abs(reference);
  const double difference = std::abs(ours - reference);

  BOOST_TEST_CONTEXT(what << ": ours = " << ours << ", ROOT = " << reference
                          << ", diff = " << difference
                          << ", allowed = " << allowed) {
    BOOST_CHECK_LE(difference, allowed);
  }
}

}  // namespace

BOOST_AUTO_TEST_SUITE(GaussianFitRootBaselineSuite)

// Confirms the modelling convention Core assumes: ROOT's "gaus" amplitude is
// compared directly against the bin content, with no bin-width factor and no
// integration of the model over the bin. Core profiles the amplitude out
// analytically as A = N / sum_i g_i, so if that matches ROOT's fitted par[0]
// the conventions agree. This underpins every other comparison here.
BOOST_AUTO_TEST_CASE(AmplitudeConvention_MatchesRoot) {
  const Histogram1 hist =
      sampled("amplitude", 50000, 0.0, 1.0, 100, -8.0, 8.0, 201);
  const auto rootHist = toRoot(hist);

  TFitResultPtr result = rootHist->Fit("gaus", "LSQ0", nullptr);
  BOOST_REQUIRE(result.Get() != nullptr);
  BOOST_REQUIRE_EQUAL(result->Status() % 1000, 0);

  const double mean = result->Parameter(1);
  const double sigma = result->Parameter(2);

  // Profiled-amplitude estimate at ROOT's own best-fit mean and sigma
  const auto& axis = hist.histogram().axis(0);
  double total = 0;
  double shapeSum = 0;
  for (int i = 0; i < axis.size(); ++i) {
    const double centre = 0.5 * (axis.bin(i).lower() + axis.bin(i).upper());
    const double z = (centre - mean) / sigma;
    total += hist.binContent({i});
    shapeSum += std::exp(-0.5 * z * z);
  }
  const double profiledAmplitude = total / shapeSum;

  BOOST_TEST_MESSAGE("ROOT par[0] = " << result->Parameter(0)
                                      << ", profiled A = "
                                      << profiledAmplitude);
  BOOST_CHECK_CLOSE(profiledAmplitude, result->Parameter(0), 1.0);
}

BOOST_AUTO_TEST_CASE(SingleFit_AgreesWithRoot) {
  for (auto& scenario : scenarios()) {
    BOOST_TEST_CONTEXT("scenario " << scenario.name) {
      const auto rootHist = toRoot(scenario.hist);
      const auto reference = rootFit(*rootHist);
      const auto ours = fit.fit(scenario.hist);

      // Every scenario is populated well enough for both to succeed. Requiring
      // it outright rather than skipping keeps the comparison from silently
      // becoming a no-op if one side starts failing.
      BOOST_REQUIRE_MESSAGE(reference.has_value(), "ROOT failed to fit");
      BOOST_REQUIRE_MESSAGE(ours.ok(), "Core failed to fit");

      checkAgrees("mean", ours->mean, reference->mean, reference->meanError,
                  1e-3, 0.05);
      checkAgrees("sigma", ours->sigma, reference->sigma, reference->sigmaError,
                  1e-3, 0.05);
      checkRelative("meanError", ours->meanError, reference->meanError, 0.02);
      checkRelative("sigmaError", ours->sigmaError, reference->sigmaError,
                    0.02);
    }
  }
}

BOOST_AUTO_TEST_CASE(IterativeFit_AgreesWithRoot) {
  constexpr double sigmaRange = 3.0;
  constexpr int iterations = 3;

  for (auto& scenario : scenarios()) {
    BOOST_TEST_CONTEXT("scenario " << scenario.name) {
      const auto rootHist = toRoot(scenario.hist);
      const auto reference =
          rootIterativeFit(*rootHist, sigmaRange, iterations);
      const auto ours =
          iterativeFit(fit, scenario.hist, sigmaRange, iterations);

      BOOST_REQUIRE_MESSAGE(reference.has_value(), "ROOT failed to fit");
      BOOST_REQUIRE_MESSAGE(ours.ok(), "Core failed to fit");

      // Looser than the single fit: a small difference in the first iteration
      // shifts the refit window, which compounds over iterations
      checkAgrees("mean", ours->mean, reference->mean, reference->meanError,
                  5e-3, 0.15);
      checkAgrees("sigma", ours->sigma, reference->sigma, reference->sigmaError,
                  5e-3, 0.15);
      checkRelative("meanError", ours->meanError, reference->meanError, 0.05);
      checkRelative("sigmaError", ours->sigmaError, reference->sigmaError,
                    0.05);
    }
  }
}

BOOST_AUTO_TEST_CASE(RestrictedRange_AgreesWithRoot) {
  // Exercises the range restriction on its own, so a disagreement here points
  // at bin selection rather than at the minimiser
  const Histogram1 hist =
      sampled("restricted", 30000, 0.4, 1.2, 80, -8.0, 8.0, 202);
  const auto rootHist = toRoot(hist);

  for (const double halfWidth : {1.0, 2.0, 3.5}) {
    BOOST_TEST_CONTEXT("half width " << halfWidth) {
      const double xMin = 0.4 - halfWidth;
      const double xMax = 0.4 + halfWidth;

      TFitResultPtr result =
          rootHist->Fit("gaus", "LRSQ0", nullptr, xMin, xMax);
      BOOST_REQUIRE(result.Get() != nullptr);
      BOOST_REQUIRE_EQUAL(result->Status() % 1000, 0);

      const auto ours = fit.fit(hist, xMin, xMax);
      BOOST_REQUIRE(ours.ok());

      checkAgrees("mean", ours->mean, result->Parameter(1), result->ParError(1),
                  1e-3, 0.05);
      checkAgrees("sigma", ours->sigma, result->Parameter(2),
                  result->ParError(2), 1e-3, 0.05);
    }
  }
}

BOOST_AUTO_TEST_CASE(DegenerateInputs_NeitherSucceedsWrongly) {
  // Where Core refuses to fit, ROOT must not be producing a usable answer we
  // are throwing away. ROOT is far more willing to return a "converged" status
  // on nonsense, so this only asserts that Core declines.
  auto axis = AxisVariant(BoostRegularAxis(20, -5.0, 5.0, "x"));

  Histogram1 empty("empty", "empty", {axis});
  const auto emptyResult = fit.fit(empty);
  BOOST_CHECK(!emptyResult.ok());
  BOOST_CHECK_EQUAL(emptyResult.error(), GaussianFitError::EmptyRange);

  Histogram1 spike("spike", "spike", {axis});
  spike.setBinContent({10}, 500.0);
  const auto spikeResult = fit.fit(spike);
  BOOST_CHECK(!spikeResult.ok());
  BOOST_CHECK_EQUAL(spikeResult.error(), GaussianFitError::TooFewNonEmptyBins);

  Histogram1 two("two", "two", {axis});
  two.setBinContent({9}, 100.0);
  two.setBinContent({10}, 120.0);
  const auto twoResult = fit.fit(two);
  BOOST_CHECK(!twoResult.ok());
  BOOST_CHECK_EQUAL(twoResult.error(), GaussianFitError::TooFewNonEmptyBins);
}

// `ActsPlugins::RootHistogramFit` satisfies
// `Acts::Experimental::GaussianFitter` and is otherwise uncalled in the
// codebase, so this is the only place proving that the generic profile
// extraction really does work with a second, ROOT-backed fitter and not just
// with `GaussianHistogramFit`.
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

  const ActsPlugins::RootHistogramFit rootFitter;
  const auto profiles = Acts::Experimental::extractMeanWidthProfiles(
      rootFitter, hist, "mean", "width");

  BOOST_CHECK_EQUAL(profiles.fitFailureFraction, 0.0);
  for (int i = 0; i < nEtaBins; ++i) {
    BOOST_CHECK_CLOSE(profiles.width.value({i}), 0.5 + 0.3 * i, 5.0);
  }
}

BOOST_AUTO_TEST_CASE(ValueHistogram1D_ConvertsWithErrors) {
  std::vector<double> edges = {0.0, 1.0, 3.0, 7.0};
  auto axis = AxisVariant(BoostVariableAxis(edges, "eta"));
  ValueHistogram1 hist("resmean_d0_vs_eta", "Mean", {axis});

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

BOOST_AUTO_TEST_CASE(ValueHistogram2D_ConvertsWithErrors) {
  auto xAxis = AxisVariant(BoostRegularAxis(2, 0.0, 2.0, "eta"));
  auto yAxis = AxisVariant(BoostRegularAxis(3, 0.0, 3.0, "pt"));
  ValueHistogram2 hist("reswidth_d0_vs_eta_pt", "Width", {xAxis, yAxis});

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
