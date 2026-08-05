// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/GaussianHistogramFit.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/GaussianHistogramFitError.hpp"
#include "Acts/Utilities/detail/NumericalMinimization.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <vector>

namespace Acts::Experimental {

namespace {

/// The bins entering a single fit, flattened out of the histogram
struct FitBins {
  std::vector<double> centres;
  std::vector<double> counts;
  /// Sum of all counts
  double total = 0;
  /// Number of bins with non-zero content
  std::size_t nonEmpty = 0;
  /// Span between the first and last selected bin centre
  double span = 0;
};

/// Collect the bins whose centre lies within `[xMin, xMax]`
///
/// The closed interval and the use of bin centres match how ROOT restricts a
/// fit range, so that a range boundary lands on the same side of a bin here as
/// it does there.
FitBins selectBins(const Histogram1& hist, double xMin, double xMax) {
  const auto& axis = hist.histogram().axis(0);

  FitBins bins;
  for (int i = 0; i < axis.size(); ++i) {
    const double centre = 0.5 * (axis.bin(i).lower() + axis.bin(i).upper());
    if (centre < xMin || centre > xMax) {
      continue;
    }

    const double count = hist.binContent({i});
    bins.centres.push_back(centre);
    bins.counts.push_back(count);
    bins.total += count;
    if (count != 0) {
      ++bins.nonEmpty;
    }
  }

  if (!bins.centres.empty()) {
    bins.span = bins.centres.back() - bins.centres.front();
  }

  return bins;
}

/// Negative log-likelihood with the amplitude profiled out, up to a constant
///
/// Evaluates `F(m, s) = N log(sum_i g_i) + 0.5 sum_i n_i z_i^2`. The first term
/// is computed in a shifted form so that `sum_i g_i` cannot underflow when the
/// trial Gaussian sits far away from every bin.
double negLogLikelihood(const FitBins& bins, double mean, double sigma) {
  if (!(sigma > 0) || !std::isfinite(mean)) {
    return std::numeric_limits<double>::infinity();
  }

  // Smallest 0.5 * z^2 over all bins, used as the log-sum-exp shift
  double minHalfZSq = std::numeric_limits<double>::infinity();
  for (const double centre : bins.centres) {
    const double z = (centre - mean) / sigma;
    minHalfZSq = std::min(minHalfZSq, 0.5 * z * z);
  }
  if (!std::isfinite(minHalfZSq)) {
    return std::numeric_limits<double>::infinity();
  }

  double shiftedSum = 0;
  double weightedChiSq = 0;
  for (std::size_t i = 0; i < bins.centres.size(); ++i) {
    const double z = (bins.centres[i] - mean) / sigma;
    const double halfZSq = 0.5 * z * z;
    shiftedSum += std::exp(-(halfZSq - minHalfZSq));
    weightedChiSq += bins.counts[i] * halfZSq;
  }

  const double logSum = -minHalfZSq + std::log(shiftedSum);
  const double value = bins.total * logSum + weightedChiSq;

  return std::isfinite(value) ? value : std::numeric_limits<double>::infinity();
}

/// Starting values for the search, following ROOT's `H1InitGaus`: the
/// count-weighted mean and RMS of the selected bins.
///
/// @return `(mean, sigma)`, or `std::nullopt` if no sensible seed exists
std::optional<std::array<double, 2>> initialGuess(const FitBins& bins) {
  if (bins.total <= 0) {
    return std::nullopt;
  }

  double sumX = 0;
  double sumXSq = 0;
  for (std::size_t i = 0; i < bins.centres.size(); ++i) {
    // Negative contents would corrupt the moments; they cannot arise from
    // counting but can from a hand-filled histogram
    const double count = std::max(0.0, bins.counts[i]);
    sumX += count * bins.centres[i];
    sumXSq += count * bins.centres[i] * bins.centres[i];
  }

  const double mean = sumX / bins.total;
  const double variance = sumXSq / bins.total - mean * mean;

  // A vanishing or negative variance means the counts sit in a single bin, or
  // rounding has eaten the spread; fall back on a fraction of the fit range
  const double sigma = (variance > 0) ? std::sqrt(variance) : 0.25 * bins.span;
  if (!(sigma > 0) || !std::isfinite(mean)) {
    return std::nullopt;
  }

  return std::array<double, 2>{mean, sigma};
}

/// Nelder-Mead simplex search over `(mean, log sigma)`
///
/// The log parametrisation keeps sigma positive without a barrier term and
/// makes the simplex scale free.
///
/// @param relativeStep Initial simplex step, relative to the seed sigma
/// @return The minimising `(mean, sigma)`, or the minimiser's error if the
///         search did not converge within the iteration cap
Result<std::array<double, 2>> minimise(const FitBins& bins, double meanSeed,
                                       double sigmaSeed, double relativeStep) {
  const auto objective = [&bins](const Vector<2>& p) {
    return negLogLikelihood(bins, p(0), std::exp(p(1)));
  };

  const Vector<2> start{meanSeed, std::log(sigmaSeed)};
  const Vector<2> steps{relativeStep * sigmaSeed, relativeStep};

  const Result<Vector<2>> optimum =
      Acts::detail::nelderMead<2>(objective, start, steps);
  if (!optimum.ok()) {
    return Result<std::array<double, 2>>::failure(optimum.error());
  }

  return std::array<double, 2>{(*optimum)(0), std::exp((*optimum)(1))};
}

/// Parameter uncertainties from the curvature of the profiled likelihood
///
/// The covariance is the inverse Hessian of the negative log-likelihood. Taking
/// it from the amplitude-profiled objective is legitimate: the inverse Hessian
/// of a profile likelihood equals the corresponding block of the full
/// three-parameter covariance matrix, so these are directly comparable to
/// ROOT's `ParError` for mean and sigma.
///
/// @param relativeStep Finite-difference step, relative to the asymptotic
///                     parameter uncertainty
/// @return `(meanError, sigmaError)`, or the minimiser's error if the
///         curvature is not that of a genuine minimum
Result<std::array<double, 2>> parameterErrors(const FitBins& bins, double mean,
                                              double sigma,
                                              double relativeStep) {
  const auto objective = [&bins](const Vector<2>& p) {
    return negLogLikelihood(bins, p(0), p(1));
  };

  // Large enough to lift the second difference well clear of floating point
  // noise, small enough to stay inside the quadratic region around the minimum
  const Vector<2> point{mean, sigma};
  const Vector<2> steps{relativeStep * sigma / std::sqrt(bins.total),
                        relativeStep * sigma / std::sqrt(2 * bins.total)};

  const auto covariance =
      Acts::detail::numericalCovariance<2>(objective, point, steps);
  if (!covariance.ok()) {
    return Result<std::array<double, 2>>::failure(covariance.error());
  }

  const double varMean = (*covariance)(0, 0);
  const double varSigma = (*covariance)(1, 1);
  if (!(varMean > 0) || !(varSigma > 0)) {
    return Result<std::array<double, 2>>::failure(
        Acts::detail::NumericalMinimizationError::NotPositiveDefinite);
  }

  return std::array<double, 2>{std::sqrt(varMean), std::sqrt(varSigma)};
}

Result<GaussianHistogramFitResult> fitBins(
    const FitBins& bins, const GaussianHistogramFit::Config& config) {
  if (bins.total <= 0) {
    return Result<GaussianHistogramFitResult>::failure(
        GaussianHistogramFitError::EmptyRange);
  }
  if (bins.nonEmpty < config.minNonEmptyBins) {
    return Result<GaussianHistogramFitResult>::failure(
        GaussianHistogramFitError::TooFewNonEmptyBins);
  }

  const std::optional<std::array<double, 2>> seed = initialGuess(bins);
  if (!seed.has_value()) {
    return Result<GaussianHistogramFitResult>::failure(
        GaussianHistogramFitError::NoValidSeed);
  }

  const auto optimum =
      minimise(bins, (*seed)[0], (*seed)[1], config.relativeStep);
  if (!optimum.ok()) {
    return Result<GaussianHistogramFitResult>::failure(optimum.error());
  }

  const double mean = (*optimum)[0];
  const double sigma = (*optimum)[1];
  if (!std::isfinite(mean) || !std::isfinite(sigma) || sigma <= 0) {
    return Result<GaussianHistogramFitResult>::failure(
        GaussianHistogramFitError::NonFiniteParameters);
  }

  const auto errors = parameterErrors(bins, mean, sigma, config.relativeStep);
  if (!errors.ok()) {
    return Result<GaussianHistogramFitResult>::failure(errors.error());
  }

  return GaussianHistogramFitResult{mean, sigma, (*errors)[0], (*errors)[1]};
}

}  // namespace

Result<GaussianHistogramFitResult> GaussianHistogramFit::fit(
    const Histogram1& hist) const {
  const double infinity = std::numeric_limits<double>::infinity();
  return fitBins(selectBins(hist, -infinity, infinity), m_config);
}

Result<GaussianHistogramFitResult> GaussianHistogramFit::fit(
    const Histogram1& hist, double xMin, double xMax) const {
  return fitBins(selectBins(hist, xMin, xMax), m_config);
}

}  // namespace Acts::Experimental
