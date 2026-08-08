// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Histogram.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <array>
#include <functional>
#include <optional>
#include <string>
#include <tuple>
#include <utility>

namespace ActsExamples {

/// @brief Outcome of a Gaussian fit to a 1D histogram: `(mean, sigma,
///        meanError, sigmaError)`
using HistogramFitResult = std::tuple<double, double, double, double>;

/// @brief Fit range `[xMin, xMax]`, closed, selected by bin centre
using HistogramFitRange = std::pair<double, double>;

/// A single Gaussian fit to a 1D histogram, optionally restricted to a range
///
/// Any backend -- @c ActsExamples::gaussianHistogramFit, a callable
/// `ActsPlugins::RootHistogramFit`, or a Python callable -- can be adapted to
/// this signature.
using HistogramFitFunction = std::function<std::optional<HistogramFitResult>(
    const Acts::Experimental::Histogram1&, std::optional<HistogramFitRange>)>;

/// @brief Mean and width profiles extracted from a histogram of dimension
///        @c Dim + 1
///
/// @tparam Dim Number of dimensions of the profiled outer axes
template <std::size_t Dim>
struct MeanWidthProfiles {
  /// Fitted mean per bin of the outer axes
  Acts::Experimental::ValueHistogram<Dim> mean;
  /// Fitted width (sigma) per bin of the outer axes
  Acts::Experimental::ValueHistogram<Dim> width;
  /// Fraction of bins where a fit was attempted but failed
  double fitFailureFraction{};
};

/// Mean and width profiles extracted from a 2D histogram
using MeanWidthProfiles1 = MeanWidthProfiles<1>;
/// Mean and width profiles extracted from a 3D histogram
using MeanWidthProfiles2 = MeanWidthProfiles<2>;

/// Fit a Gaussian repeatedly, narrowing the fit range around the peak
///
/// The first fit is unrestricted. Each subsequent fit is restricted to
/// @f$ m \pm \mathrm{sigmaRange} \cdot s @f$ using the previous iteration's
/// parameters. The returned uncertainties come from the final iteration.
///
/// @param fitFn The single-range fit function to iterate
/// @param hist The histogram to fit
/// @param sigmaRange Half-width of the restricted range, in fitted sigmas
/// @param iterations Total number of fits, including the initial unrestricted
///                   one; values below 1 are treated as 1
/// @param logger Logger for diagnostics on failed iterations
/// @return The fit result, or @c std::nullopt if any iteration failed
std::optional<HistogramFitResult> iterativeFit(
    const HistogramFitFunction& fitFn,
    const Acts::Experimental::Histogram1& hist, double sigmaRange,
    int iterations, const Acts::Logger& logger = Acts::getDummyLogger());

/// Fit a Gaussian to every slice of a histogram along its last axis
///
/// For each bin of the outer axes (every axis but the last), the distribution
/// along the last axis is fitted with @c iterativeFit and the resulting mean
/// and sigma are stored, with their uncertainties, in the corresponding
/// output bin.
///
/// @param fitFn The single-range fit function to use for every slice
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
template <std::size_t Dim>
MeanWidthProfiles<Dim - 1> extractMeanWidthProfiles(
    const HistogramFitFunction& fitFn,
    const Acts::Experimental::Histogram<Dim>& hist, const std::string& meanName,
    const std::string& widthName, int minEntriesForFit = 5,
    double sigmaRange = 3.0, int iterations = 3,
    const Acts::Logger& logger = Acts::getDummyLogger()) {
  constexpr std::size_t OuterDim = Dim - 1;

  std::array<Acts::Experimental::AxisVariant, OuterDim> axes{};
  std::array<int, OuterDim> outerSizes{};
  int totalOuterBins = 1;
  for (std::size_t d = 0; d < OuterDim; ++d) {
    axes[d] = hist.histogram().axis(d);
    outerSizes[d] = hist.histogram().axis(d).size();
    totalOuterBins *= outerSizes[d];
  }

  MeanWidthProfiles<OuterDim> profiles{
      Acts::Experimental::ValueHistogram<OuterDim>(
          meanName, hist.title() + " mean", axes),
      Acts::Experimental::ValueHistogram<OuterDim>(
          widthName, hist.title() + " width", axes),
      0.0};

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
    const Acts::Experimental::Histogram1 slice =
        Acts::Experimental::detail::sliceLastAxis<Dim>(hist, outerBins);
    if (Acts::Experimental::totalContent(slice) < minEntriesForFit) {
      // Too few entries: skipped, does not count as a fit failure
      continue;
    }

    const std::optional<HistogramFitResult> result =
        iterativeFit(fitFn, slice, sigmaRange, iterations, logger);
    if (!result.has_value()) {
      ++fitFailures;
      continue;
    }

    const auto& [mean, sigma, meanError, sigmaError] = *result;
    profiles.mean.setBin(outerBins, mean, meanError);
    profiles.width.setBin(outerBins, sigma, sigmaError);
  }

  profiles.fitFailureFraction =
      (totalOuterBins > 0) ? static_cast<double>(fitFailures) / totalOuterBins
                           : 0;

  return profiles;
}

}  // namespace ActsExamples
