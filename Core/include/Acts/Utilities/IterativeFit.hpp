// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/GaussianHistogramFit.hpp"
#include "Acts/Utilities/Histogram.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <array>
#include <concepts>
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
concept GaussianHistogramFitter = requires(
    const F& fitter, const Histogram1& hist, double xMin, double xMax) {
  { fitter.fit(hist) } -> std::same_as<Result<GaussianHistogramFitResult>>;
  {
    fitter.fit(hist, xMin, xMax)
  } -> std::same_as<Result<GaussianHistogramFitResult>>;
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
template <GaussianHistogramFitter Fitter>
Result<GaussianHistogramFitResult> iterativeFit(
    const Fitter& fitter, const Histogram1& hist, double sigmaRange,
    int iterations, const Logger& logger = getDummyLogger()) {
  Result<GaussianHistogramFitResult> result = fitter.fit(hist);
  if (!result.ok()) {
    ACTS_DEBUG("Failed to fit initial Gaussian to '"
               << hist.name() << "': " << result.error().message());
    return result;
  }

  for (int i = 0; i < iterations - 1; ++i) {
    const double xMin = result->mean - sigmaRange * result->sigma;
    const double xMax = result->mean + sigmaRange * result->sigma;

    Result<GaussianHistogramFitResult> restricted =
        fitter.fit(hist, xMin, xMax);
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
template <std::size_t Dim, GaussianHistogramFitter Fitter>
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

    const Result<GaussianHistogramFitResult> result =
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

}  // namespace Acts::Experimental
