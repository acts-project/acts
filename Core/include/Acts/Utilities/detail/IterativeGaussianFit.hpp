// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/GaussianFitResult.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <concepts>
#include <limits>
#include <optional>
#include <string>

namespace Acts::detail {

/// A single Gaussian fit over a range `[xMin, xMax]`, as consumed by
/// @c iterativeGaussianFit
///
/// Passing `xMin = -infinity, xMax = +infinity` requests an unrestricted fit.
template <typename Callable>
concept RangedGaussianFitter =
    requires(const Callable& fit, double xMin, double xMax) {
      {
        fit(xMin, xMax)
      } -> std::same_as<std::optional<Acts::Experimental::GaussianFitResult>>;
    };

/// Fit a Gaussian repeatedly, narrowing the fit range around the peak
///
/// This is the part of the iterative Gaussian fit that does not depend on how
/// a single fit is actually performed, so it is shared between
/// @c Acts::Experimental::GaussianHistogramFit (a Poisson-likelihood fit on a
/// boost-histogram) and @c ActsPlugin::RootHistogramFit (ROOT's `TH1::Fit`).
///
/// The first fit is unrestricted. Each subsequent fit is restricted to
/// @f$ m \pm \mathrm{sigmaRange} \cdot s @f$ using the previous iteration's
/// parameters. The returned uncertainties come from the final iteration.
///
/// @param fit Single-shot fit over a range; called with
///            `(-infinity, +infinity)` for the initial, unrestricted fit
/// @param sigmaRange Half-width of the restricted range, in fitted sigmas
/// @param iterations Total number of fits, including the initial unrestricted
///                   one; values below 1 are treated as 1
/// @param logger Logger for diagnostics on failed iterations
/// @param name Name of the histogram being fitted, for diagnostics only
/// @return The fit result, or `std::nullopt` if any iteration failed
template <RangedGaussianFitter Callable>
std::optional<Acts::Experimental::GaussianFitResult> iterativeGaussianFit(
    Callable&& fit, double sigmaRange, int iterations, const Logger& logger,
    const std::string& name) {
  const double infinity = std::numeric_limits<double>::infinity();

  std::optional<Acts::Experimental::GaussianFitResult> result =
      fit(-infinity, infinity);
  if (!result.has_value()) {
    ACTS_DEBUG("Failed to fit initial Gaussian to '" << name << "'");
    return std::nullopt;
  }

  for (int i = 0; i < iterations - 1; ++i) {
    const double xMin = result->mean - sigmaRange * result->sigma;
    const double xMax = result->mean + sigmaRange * result->sigma;

    const std::optional<Acts::Experimental::GaussianFitResult> restricted =
        fit(xMin, xMax);
    if (!restricted.has_value()) {
      ACTS_DEBUG("Failed to fit iteration " << i << " Gaussian to '" << name
                                            << "'");
      return std::nullopt;
    }

    result = restricted;
  }

  return result;
}

/// Loop over the outer bins of a histogram, fitting a Gaussian to the slice
/// along its last axis and storing the mean/width in output bins
///
/// This is the part of mean/width profile extraction that does not depend on
/// the concrete histogram type, so it is shared between
/// @c Acts::Experimental::GaussianHistogramFit::extractMeanWidthProfiles
/// (boost-histogram slices) and
/// @c ActsPlugin::RootHistogramFit::extractMeanWidthProfiles (ROOT TH2F/TH3F
/// projections): both skip low-statistics slices, fit the rest, and record
/// either a mean/width pair or a fit failure.
///
/// @param totalOuterBins Number of outer bins to loop over
/// @param trySlice Given a flat outer-bin index, returns `std::nullopt` if
///                 the slice should be skipped (too few entries), otherwise
///                 the outcome of fitting that slice: a fit result on
///                 success, or `std::nullopt` if the fit itself failed
/// @param storeResult Given a flat outer-bin index and a successful fit
///                    result, stores it in the output histograms
/// @return The fraction of attempted (non-skipped) fits that failed
template <typename TrySlice, typename StoreResult>
double extractMeanWidthProfilesImpl(int totalOuterBins, TrySlice&& trySlice,
                                    StoreResult&& storeResult) {
  int fitFailures = 0;
  for (int flat = 0; flat < totalOuterBins; ++flat) {
    const std::optional<std::optional<Acts::Experimental::GaussianFitResult>>
        outcome = trySlice(flat);
    if (!outcome.has_value()) {
      // Skipped: too few entries, does not count as a fit failure
      continue;
    }
    if (!outcome->has_value()) {
      ++fitFailures;
      continue;
    }

    storeResult(flat, **outcome);
  }

  return (totalOuterBins > 0)
             ? static_cast<double>(fitFailures) / totalOuterBins
             : 0;
}

}  // namespace Acts::detail
