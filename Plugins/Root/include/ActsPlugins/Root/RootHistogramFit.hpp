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

#include <memory>
#include <optional>
#include <string>
#include <tuple>

class TH1;
class TH1F;
class TH2F;
class TH3F;

namespace ActsPlugins {

/// Fit a Gaussian to a ROOT histogram via `TH1::Fit`
///
/// Provides the same interface as `Acts::Experimental::GaussianHistogramFit`,
/// backed by ROOT's own least-squares `"gaus"` fitter instead of Core's
/// Poisson-likelihood optimiser. The iterative-narrowing and mean/width
/// profile extraction loops are shared with the Core implementation (see
/// `Acts::detail::iterativeGaussianFit` and
/// `Acts::detail::extractMeanWidthProfilesImpl`); only the single-range fit
/// itself is ROOT-specific.
class RootHistogramFit {
 public:
  RootHistogramFit() = default;

  /// Fit a Gaussian to a ROOT histogram
  ///
  /// @param hist The histogram to fit
  /// @return The fit result, or `std::nullopt` if the fit could not be performed
  std::optional<Acts::Experimental::GaussianFitResult> fit(TH1& hist) const;

  /// Fit a Gaussian to a ROOT histogram, restricted to a range
  ///
  /// @param hist The histogram to fit
  /// @param xMin Lower end of the fit range (inclusive)
  /// @param xMax Upper end of the fit range (inclusive)
  /// @return The fit result, or `std::nullopt` if the fit could not be performed
  std::optional<Acts::Experimental::GaussianFitResult> fit(TH1& hist,
                                                           double xMin,
                                                           double xMax) const;

  /// Fit a Gaussian repeatedly, narrowing the fit range around the peak
  ///
  /// @param hist The histogram to fit
  /// @param sigmaRange Half-width of the restricted range, in fitted sigmas
  /// @param iterations Total number of fits, including the initial unrestricted
  ///                   one; values below 1 are treated as 1
  /// @param logger Logger for diagnostics on failed iterations
  /// @return The fit result, or `std::nullopt` if any iteration failed
  std::optional<Acts::Experimental::GaussianFitResult> iterativeFit(
      TH1& hist, double sigmaRange, int iterations,
      const Acts::Logger& logger = Acts::getDummyLogger()) const;

  /// Extract 1D mean/width profiles from a 2D histogram
  ///
  /// For each X bin, the Y projection is fitted with @c iterativeFit and the
  /// resulting mean and sigma are stored, with their uncertainties, in the
  /// corresponding output bin.
  ///
  /// @param hist2d The input 2D histogram to analyze
  /// @param meanName The name for the output mean profile histogram
  /// @param widthName The name for the output width profile histogram
  /// @param minEntriesForFit Minimum number of entries in a projection to attempt a fit
  /// @param sigmaRange The range in sigma for the iterative Gaussian fit
  /// @param iterations The maximum number of iterations for the iterative Gaussian fit
  /// @param logger Logger for debug messages
  /// @return pair of unique pointers to the mean and width TH1F histograms and a fit failure fraction
  std::tuple<std::unique_ptr<TH1F>, std::unique_ptr<TH1F>, double>
  extractMeanWidthProfiles(
      const TH2F& hist2d, const std::string& meanName,
      const std::string& widthName, int minEntriesForFit = 5,
      double sigmaRange = 3.0, int iterations = 3,
      const Acts::Logger& logger = Acts::getDummyLogger()) const;

  /// Extract 2D mean/width profiles from a 3D histogram
  ///
  /// For each (X, Y) bin, the Z projection is fitted with @c iterativeFit and
  /// the resulting mean and sigma are stored, with their uncertainties, in the
  /// corresponding output bin.
  ///
  /// @param hist3d The input 3D histogram to analyze
  /// @param meanName The name for the output mean profile histogram
  /// @param widthName The name for the output width profile histogram
  /// @param minEntriesForFit Minimum number of entries in a projection to attempt a fit
  /// @param sigmaRange The range in sigma for the iterative Gaussian fit
  /// @param iterations The maximum number of iterations for the iterative Gaussian fit
  /// @param logger Logger for debug messages
  /// @return pair of unique pointers to the mean and width TH2F histograms and a fit failure fraction
  std::tuple<std::unique_ptr<TH2F>, std::unique_ptr<TH2F>, double>
  extractMeanWidthProfiles(
      const TH3F& hist3d, const std::string& meanName,
      const std::string& widthName, int minEntriesForFit = 5,
      double sigmaRange = 3.0, int iterations = 3,
      const Acts::Logger& logger = Acts::getDummyLogger()) const;
};

}  // namespace ActsPlugins
