// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/GaussianFitResult.hpp"
#include "Acts/Utilities/Histogram.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <array>
#include <concepts>
#include <limits>
#include <string>

namespace Acts::Experimental {

/// @brief Mean and width profiles extracted from a histogram of dimension
///        @c Dim + 1
///
/// @tparam Dim Number of dimensions of the profiled outer axes
template <std::size_t Dim>
struct MeanWidthProfiles {
  /// Fitted mean per bin of the outer axes
  ValueHistogram<Dim> mean;
  /// Fitted width (sigma) per bin of the outer axes
  ValueHistogram<Dim> width;
  /// Fraction of bins where a fit was attempted but failed
  double fitFailureFraction{};
};

/// Mean and width profiles extracted from a 2D histogram
using MeanWidthProfiles1 = MeanWidthProfiles<1>;
/// Mean and width profiles extracted from a 3D histogram
using MeanWidthProfiles2 = MeanWidthProfiles<2>;

/// A single Gaussian fit to a 1D histogram, optionally restricted to a range
///
/// This is the interface shared by @c Acts::Experimental::GaussianHistogramFit
/// (a Poisson-likelihood fit on a boost-histogram) and
/// @c ActsPlugins::RootHistogramFit (ROOT's `TH1::Fit`). Anything built on top
/// of a single fit -- iterative narrowing, mean/width profile extraction -- is
/// written once, generically, against this concept; see @c iterativeFit and
/// @c extractMeanWidthProfiles below.
template <typename F>
concept GaussianFitter = requires(const F& fitter, const Histogram1& hist,
                                  double xMin, double xMax) {
  { fitter.fit(hist) } -> std::same_as<Result<GaussianFitResult>>;
  { fitter.fit(hist, xMin, xMax) } -> std::same_as<Result<GaussianFitResult>>;
};

/// Fit a Gaussian repeatedly, narrowing the fit range around the peak
///
/// The first fit is unrestricted. Each subsequent fit is restricted to
/// @f$ m \pm \mathrm{sigmaRange} \cdot s @f$ using the previous iteration's
/// parameters. The returned uncertainties come from the final iteration.
///
/// @param fitter The single-range fitter to iterate
/// @param hist The histogram to fit
/// @param sigmaRange Half-width of the restricted range, in fitted sigmas
/// @param iterations Total number of fits, including the initial unrestricted
///                   one; values below 1 are treated as 1
/// @param logger Logger for diagnostics on failed iterations
/// @return The fit result, or an error if any iteration failed
template <GaussianFitter Fitter>
Result<GaussianFitResult> iterativeFit(
    const Fitter& fitter, const Histogram1& hist, double sigmaRange,
    int iterations, const Logger& logger = getDummyLogger()) {
  Result<GaussianFitResult> result = fitter.fit(hist);
  if (!result.ok()) {
    ACTS_DEBUG("Failed to fit initial Gaussian to '"
               << hist.name() << "': " << result.error().message());
    return result;
  }

  for (int i = 0; i < iterations - 1; ++i) {
    const double xMin = result->mean - sigmaRange * result->sigma;
    const double xMax = result->mean + sigmaRange * result->sigma;

    Result<GaussianFitResult> restricted = fitter.fit(hist, xMin, xMax);
    if (!restricted.ok()) {
      ACTS_DEBUG("Failed to fit iteration "
                 << i << " Gaussian to '" << hist.name()
                 << "': " << restricted.error().message());
      return restricted;
    }

    result = std::move(restricted);
  }

  return result;
}

/// Fit a Gaussian to every slice of a histogram along its last axis
///
/// For each bin of the outer axes (every axis but the last), the distribution
/// along the last axis is fitted with @c iterativeFit and the resulting mean
/// and sigma are stored, with their uncertainties, in the corresponding
/// output bin.
///
/// @param fitter The single-range fitter to use for every slice
/// @param hist The histogram to profile
/// @param meanName Name for the mean output histogram
/// @param widthName Name for the width output histogram
/// @param minEntriesForFit Slices with fewer entries are skipped
/// @param sigmaRange Half-width of the iterative refit range, in fitted sigmas
/// @param iterations Number of fits per slice, including the unrestricted one
/// @param logger Logger for diagnostics on failed fits
/// @return The mean and width profiles and the fit failure fraction
/// @note Skipped slices leave their output bins empty and do not count towards
///       @c fitFailureFraction, which reports only genuine fit failures.
template <std::size_t Dim, GaussianFitter Fitter>
MeanWidthProfiles<Dim - 1> extractMeanWidthProfiles(
    const Fitter& fitter, const Histogram<Dim>& hist,
    const std::string& meanName, const std::string& widthName,
    int minEntriesForFit = 5, double sigmaRange = 3.0, int iterations = 3,
    const Logger& logger = getDummyLogger()) {
  constexpr std::size_t OuterDim = Dim - 1;

  std::array<AxisVariant, OuterDim> axes{};
  std::array<int, OuterDim> outerSizes{};
  int totalOuterBins = 1;
  for (std::size_t d = 0; d < OuterDim; ++d) {
    axes[d] = hist.histogram().axis(d);
    outerSizes[d] = hist.histogram().axis(d).size();
    totalOuterBins *= outerSizes[d];
  }

  MeanWidthProfiles<OuterDim> profiles{
      ValueHistogram<OuterDim>(meanName, hist.title() + " mean", axes),
      ValueHistogram<OuterDim>(widthName, hist.title() + " width", axes), 0.0};

  // Unravel a flat outer index into per-axis indices, last outer axis
  // fastest, matching the nested-loop order of the original
  // per-dimension overloads
  const auto unravel = [&](int flat) {
    std::array<int, OuterDim> outerBins{};
    int remaining = flat;
    for (std::size_t d = OuterDim; d-- > 0;) {
      outerBins[d] = remaining % outerSizes[d];
      remaining /= outerSizes[d];
    }
    return outerBins;
  };

  int fitFailures = 0;
  for (int flat = 0; flat < totalOuterBins; ++flat) {
    const std::array<int, OuterDim> outerBins = unravel(flat);
    const Histogram1 slice = detail::sliceLastAxis<Dim>(hist, outerBins);
    if (totalContent(slice) < minEntriesForFit) {
      // Too few entries: skipped, does not count as a fit failure
      continue;
    }

    const Result<GaussianFitResult> result =
        iterativeFit(fitter, slice, sigmaRange, iterations, logger);
    if (!result.ok()) {
      ++fitFailures;
      continue;
    }

    profiles.mean.setBin(outerBins, result->mean, result->meanError);
    profiles.width.setBin(outerBins, result->sigma, result->sigmaError);
  }

  profiles.fitFailureFraction =
      (totalOuterBins > 0) ? static_cast<double>(fitFailures) / totalOuterBins
                           : 0;

  return profiles;
}

/// Fit a Gaussian to a 1D histogram by maximising the Poisson likelihood
///
/// The model is @f$ \mu_i = A \exp(-(x_i - m)^2 / (2 s^2)) @f$ with @f$ x_i @f$
/// the bin centre, compared directly against the bin content. This mirrors
/// ROOT's predefined `"gaus"` under `TH1::Fit(..., "L")`: no bin-width factor
/// and no integration of the model over the bin. See @c ActsPlugin::RootHistogramFit
/// for a ROOT-backed fit with the same interface.
///
/// The amplitude @f$ A @f$ is not searched over. For fixed @f$ (m, s) @f$ its
/// maximum-likelihood value is @f$ \hat{A} = N / \sum_i g_i @f$, which reduces
/// the problem to two parameters and makes it markedly more robust. Up to a
/// constant the minimised objective is
/// @f[
///   F(m, s) = N \log \sum_i g_i + \frac{1}{2} \sum_i n_i z_i^2, \quad
///   z_i = \frac{x_i - m}{s}, \quad g_i = e^{-z_i^2 / 2}.
/// @f]
///
/// Unlike a least-squares fit, empty bins carry information here (through
/// @f$ \sum_i g_i @f$), which is what makes the likelihood the better choice
/// for sparse histograms and histograms with outliers.
class GaussianHistogramFit {
 public:
  /// Configuration for @c GaussianHistogramFit
  struct Config {
    /// Minimum number of populated bins entering a single fit. The
    /// underlying model has three parameters (amplitude, mean, sigma), so
    /// fewer bins cannot constrain it.
    std::size_t minNonEmptyBins = 3;
    /// Fractional step, relative to the seed sigma, used both to size the
    /// initial Nelder-Mead simplex and to size the finite-difference step
    /// for the parameter-error Hessian. Kept as a single value so the two
    /// stay consistent with each other.
    double relativeStep = 0.1;
  };

  /// Construct with the default configuration
  GaussianHistogramFit() = default;

  /// Construct with the given configuration
  /// @param config The fit configuration
  explicit GaussianHistogramFit(Config config) : m_config(config) {}

  /// The fit configuration
  /// @return The configuration this instance was constructed with
  const Config& config() const { return m_config; }

  /// Fit a Gaussian to a 1D histogram
  ///
  /// @param hist The histogram to fit
  /// @return The fit result, or an error if the fit could not be performed
  /// @note Under- and overflow bins are ignored.
  Result<GaussianFitResult> fit(const Histogram1& hist) const;

  /// Fit a Gaussian to the bins of a 1D histogram whose centre lies in a range
  ///
  /// @param hist The histogram to fit
  /// @param xMin Lower end of the fit range (inclusive)
  /// @param xMax Upper end of the fit range (inclusive)
  /// @return The fit result, or an error if the fit could not be performed
  /// @remark Bin selection is by bin centre on a closed interval, matching how
  ///         ROOT restricts a fit range.
  /// @note See @c fit(const Histogram1&) for the fitted model.
  Result<GaussianFitResult> fit(const Histogram1& hist, double xMin,
                                double xMax) const;

 private:
  Config m_config{};
};

}  // namespace Acts::Experimental
