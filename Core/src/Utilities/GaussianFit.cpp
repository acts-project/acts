// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/GaussianFit.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

namespace Acts::Experimental {

namespace {

/// Minimum number of populated bins. The underlying model has three
/// parameters (amplitude, mean, sigma), so fewer bins cannot constrain it.
constexpr std::size_t s_minNonEmptyBins = 3;

/// Iteration cap for the simplex search; a well behaved 2D problem converges
/// in far fewer steps, so hitting this means something is wrong.
constexpr std::size_t s_maxSimplexIterations = 2000;

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
/// makes the simplex scale free. Two parameters means a three-vertex simplex,
/// which is cheap enough that the tolerance can be set very tight.
///
/// @return The minimising `(mean, sigma)`, or `std::nullopt` if the search did
///         not converge within the iteration cap
std::optional<std::array<double, 2>> minimise(const FitBins& bins,
                                              double meanSeed,
                                              double sigmaSeed) {
  // Standard Nelder-Mead coefficients
  constexpr double reflection = 1.0;
  constexpr double expansion = 2.0;
  constexpr double contraction = 0.5;
  constexpr double shrink = 0.5;

  using Point = std::array<double, 2>;
  const auto evaluate = [&bins](const Point& p) {
    return negLogLikelihood(bins, p[0], std::exp(p[1]));
  };

  // Initial simplex: displace each coordinate on its own natural scale
  const double logSigmaSeed = std::log(sigmaSeed);
  std::array<Point, 3> simplex = {
      Point{meanSeed, logSigmaSeed},
      Point{meanSeed + 0.1 * sigmaSeed, logSigmaSeed},
      Point{meanSeed, logSigmaSeed + 0.1}};
  std::array<double, 3> values = {evaluate(simplex[0]), evaluate(simplex[1]),
                                  evaluate(simplex[2])};
  if (std::ranges::none_of(values, [](double v) { return std::isfinite(v); })) {
    return std::nullopt;
  }

  const auto sortSimplex = [&simplex, &values]() {
    std::array<std::size_t, 3> order = {0, 1, 2};
    std::ranges::sort(order, [&values](std::size_t a, std::size_t b) {
      return values[a] < values[b];
    });
    const std::array<Point, 3> sortedSimplex = {
        simplex[order[0]], simplex[order[1]], simplex[order[2]]};
    const std::array<double, 3> sortedValues = {
        values[order[0]], values[order[1]], values[order[2]]};
    simplex = sortedSimplex;
    values = sortedValues;
  };

  for (std::size_t iteration = 0; iteration < s_maxSimplexIterations;
       ++iteration) {
    sortSimplex();

    // Converged once the simplex is tiny both in objective and in parameters.
    // Requiring both avoids stopping early on a flat stretch.
    const double valueSpread = std::abs(values[2] - values[0]);
    double parameterSpread = 0;
    for (std::size_t v = 1; v < simplex.size(); ++v) {
      for (std::size_t d = 0; d < 2; ++d) {
        parameterSpread =
            std::max(parameterSpread, std::abs(simplex[v][d] - simplex[0][d]));
      }
    }
    const double scale = std::max(1.0, std::abs(values[0]));
    if (valueSpread < 1e-12 * scale && parameterSpread < 1e-10) {
      return Point{simplex[0][0], std::exp(simplex[0][1])};
    }

    // Centroid of all but the worst vertex
    Point centroid = {0, 0};
    for (std::size_t v = 0; v + 1 < simplex.size(); ++v) {
      for (std::size_t d = 0; d < 2; ++d) {
        centroid[d] += simplex[v][d] / (simplex.size() - 1);
      }
    }

    const auto along = [&centroid, &simplex](double factor) {
      Point p{};
      for (std::size_t d = 0; d < 2; ++d) {
        p[d] = centroid[d] + factor * (centroid[d] - simplex.back()[d]);
      }
      return p;
    };

    const Point reflected = along(reflection);
    const double reflectedValue = evaluate(reflected);

    if (reflectedValue < values[0]) {
      // Better than the best so far: try to push further in that direction
      const Point expanded = along(expansion);
      const double expandedValue = evaluate(expanded);
      if (expandedValue < reflectedValue) {
        simplex.back() = expanded;
        values.back() = expandedValue;
      } else {
        simplex.back() = reflected;
        values.back() = reflectedValue;
      }
      continue;
    }

    if (reflectedValue < values[1]) {
      simplex.back() = reflected;
      values.back() = reflectedValue;
      continue;
    }

    const Point contracted = along(-contraction);
    const double contractedValue = evaluate(contracted);
    if (contractedValue < values.back()) {
      simplex.back() = contracted;
      values.back() = contractedValue;
      continue;
    }

    // Nothing helped: pull every vertex towards the best one
    for (std::size_t v = 1; v < simplex.size(); ++v) {
      for (std::size_t d = 0; d < 2; ++d) {
        simplex[v][d] =
            simplex[0][d] + shrink * (simplex[v][d] - simplex[0][d]);
      }
      values[v] = evaluate(simplex[v]);
    }
  }

  return std::nullopt;
}

/// Parameter uncertainties from the curvature of the profiled likelihood
///
/// The covariance is the inverse Hessian of the negative log-likelihood. Taking
/// it from the amplitude-profiled objective is legitimate: the inverse Hessian
/// of a profile likelihood equals the corresponding block of the full
/// three-parameter covariance matrix, so these are directly comparable to
/// ROOT's `ParError` for mean and sigma.
///
/// @return `(meanError, sigmaError)`, or `std::nullopt` if the curvature is not
///         that of a genuine minimum
std::optional<std::array<double, 2>> parameterErrors(const FitBins& bins,
                                                     double mean,
                                                     double sigma) {
  // Step a tenth of the asymptotic uncertainty: large enough to lift the second
  // difference well clear of floating point noise, small enough to stay inside
  // the quadratic region around the minimum
  const double meanStep = 0.1 * sigma / std::sqrt(bins.total);
  const double sigmaStep = 0.1 * sigma / std::sqrt(2 * bins.total);
  if (!(meanStep > 0) || !(sigmaStep > 0)) {
    return std::nullopt;
  }

  const auto f = [&bins](double m, double s) {
    return negLogLikelihood(bins, m, s);
  };

  const double centre = f(mean, sigma);
  const double dMeanPlus = f(mean + meanStep, sigma);
  const double dMeanMinus = f(mean - meanStep, sigma);
  const double dSigmaPlus = f(mean, sigma + sigmaStep);
  const double dSigmaMinus = f(mean, sigma - sigmaStep);

  const double hMeanMean =
      (dMeanPlus - 2 * centre + dMeanMinus) / (meanStep * meanStep);
  const double hSigmaSigma =
      (dSigmaPlus - 2 * centre + dSigmaMinus) / (sigmaStep * sigmaStep);
  const double hMixed = (f(mean + meanStep, sigma + sigmaStep) -
                         f(mean + meanStep, sigma - sigmaStep) -
                         f(mean - meanStep, sigma + sigmaStep) +
                         f(mean - meanStep, sigma - sigmaStep)) /
                        (4 * meanStep * sigmaStep);

  const double determinant = hMeanMean * hSigmaSigma - hMixed * hMixed;
  if (!std::isfinite(determinant) || determinant <= 0 || hMeanMean <= 0 ||
      hSigmaSigma <= 0) {
    // Not a positive definite Hessian, so not a minimum we can put errors on
    return std::nullopt;
  }

  // Inverse of a symmetric 2x2 matrix
  const double varMean = hSigmaSigma / determinant;
  const double varSigma = hMeanMean / determinant;
  if (!(varMean > 0) || !(varSigma > 0)) {
    return std::nullopt;
  }

  return std::array<double, 2>{std::sqrt(varMean), std::sqrt(varSigma)};
}

std::optional<GaussianFitResult> fitBins(const FitBins& bins) {
  if (bins.total <= 0 || bins.nonEmpty < s_minNonEmptyBins) {
    return std::nullopt;
  }

  const std::optional<std::array<double, 2>> seed = initialGuess(bins);
  if (!seed.has_value()) {
    return std::nullopt;
  }

  const std::optional<std::array<double, 2>> optimum =
      minimise(bins, (*seed)[0], (*seed)[1]);
  if (!optimum.has_value()) {
    return std::nullopt;
  }

  const double mean = (*optimum)[0];
  const double sigma = (*optimum)[1];
  if (!std::isfinite(mean) || !std::isfinite(sigma) || sigma <= 0) {
    return std::nullopt;
  }

  const std::optional<std::array<double, 2>> errors =
      parameterErrors(bins, mean, sigma);
  if (!errors.has_value()) {
    return std::nullopt;
  }

  return GaussianFitResult{mean, sigma, (*errors)[0], (*errors)[1]};
}

}  // namespace

std::optional<GaussianFitResult> gaussianFit(const Histogram1& hist) {
  const double infinity = std::numeric_limits<double>::infinity();
  return fitBins(selectBins(hist, -infinity, infinity));
}

std::optional<GaussianFitResult> gaussianFit(const Histogram1& hist,
                                             double xMin, double xMax) {
  return fitBins(selectBins(hist, xMin, xMax));
}

std::optional<GaussianFitResult> iterativeGaussianFit(const Histogram1& hist,
                                                      double sigmaRange,
                                                      int iterations,
                                                      const Logger& logger) {
  std::optional<GaussianFitResult> result = gaussianFit(hist);
  if (!result.has_value()) {
    ACTS_DEBUG("Failed to fit initial Gaussian to '" << hist.name() << "'");
    return std::nullopt;
  }

  for (int i = 0; i < iterations - 1; ++i) {
    const double xMin = result->mean - sigmaRange * result->sigma;
    const double xMax = result->mean + sigmaRange * result->sigma;

    const std::optional<GaussianFitResult> restricted =
        gaussianFit(hist, xMin, xMax);
    if (!restricted.has_value()) {
      ACTS_DEBUG("Failed to fit iteration " << i << " Gaussian to '"
                                            << hist.name() << "'");
      return std::nullopt;
    }

    result = restricted;
  }

  return result;
}

namespace {

/// Total content of a 1D histogram's in-range bins
///
/// Stands in for ROOT's `TH1::GetEntries()` on a projection, which likewise
/// amounts to the summed bin contents.
double sliceEntries(const Histogram1& slice) {
  const auto& axis = slice.histogram().axis(0);

  double entries = 0;
  for (int i = 0; i < axis.size(); ++i) {
    entries += slice.binContent({i});
  }

  return entries;
}

}  // namespace

MeanWidthProfiles1 extractMeanWidthProfiles(const Histogram2& hist2d,
                                            const std::string& meanName,
                                            const std::string& widthName,
                                            int minEntriesForFit,
                                            double sigmaRange, int iterations,
                                            const Logger& logger) {
  const auto& xAxis = hist2d.histogram().axis(0);
  const std::array<AxisVariant, 1> axes = {xAxis};

  MeanWidthProfiles1 profiles{
      ValueHistogram1(meanName, hist2d.title() + " mean", axes),
      ValueHistogram1(widthName, hist2d.title() + " width", axes), 0.0};

  int fitFailures = 0;
  for (int i = 0; i < xAxis.size(); ++i) {
    const Histogram1 slice = sliceLastAxis(hist2d, i);
    if (sliceEntries(slice) < minEntriesForFit) {
      continue;
    }

    const std::optional<GaussianFitResult> fit =
        iterativeGaussianFit(slice, sigmaRange, iterations, logger);
    if (!fit.has_value()) {
      ++fitFailures;
      continue;
    }

    profiles.mean.setBin({i}, fit->mean, fit->meanError);
    profiles.width.setBin({i}, fit->sigma, fit->sigmaError);
  }

  profiles.fitFailureFraction =
      (xAxis.size() > 0) ? static_cast<double>(fitFailures) / xAxis.size() : 0;

  return profiles;
}

MeanWidthProfiles2 extractMeanWidthProfiles(const Histogram3& hist3d,
                                            const std::string& meanName,
                                            const std::string& widthName,
                                            int minEntriesForFit,
                                            double sigmaRange, int iterations,
                                            const Logger& logger) {
  const auto& xAxis = hist3d.histogram().axis(0);
  const auto& yAxis = hist3d.histogram().axis(1);
  const std::array<AxisVariant, 2> axes = {xAxis, yAxis};

  MeanWidthProfiles2 profiles{
      ValueHistogram2(meanName, hist3d.title() + " mean", axes),
      ValueHistogram2(widthName, hist3d.title() + " width", axes), 0.0};

  int fitFailures = 0;
  for (int i = 0; i < xAxis.size(); ++i) {
    for (int j = 0; j < yAxis.size(); ++j) {
      const Histogram1 slice = sliceLastAxis(hist3d, i, j);
      if (sliceEntries(slice) < minEntriesForFit) {
        continue;
      }

      const std::optional<GaussianFitResult> fit =
          iterativeGaussianFit(slice, sigmaRange, iterations, logger);
      if (!fit.has_value()) {
        ++fitFailures;
        continue;
      }

      profiles.mean.setBin({i, j}, fit->mean, fit->meanError);
      profiles.width.setBin({i, j}, fit->sigma, fit->sigmaError);
    }
  }

  const int totalBins = xAxis.size() * yAxis.size();
  profiles.fitFailureFraction =
      (totalBins > 0) ? static_cast<double>(fitFailures) / totalBins : 0;

  return profiles;
}

}  // namespace Acts::Experimental
