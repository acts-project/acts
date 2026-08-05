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
#include "Acts/Utilities/Result.hpp"

#include <string>

namespace ActsPlugins {

/// Fit a Gaussian to a histogram via ROOT's `TH1::Fit`
///
/// Satisfies `Acts::Experimental::GaussianHistogramFitter`, so it can be used
/// with `Acts::Experimental::iterativeFit` and
/// `Acts::Experimental::extractMeanWidthProfiles` in place of
/// `Acts::Experimental::GaussianHistogramFit`, backed by ROOT's own `"gaus"`
/// fitter instead of Core's optimiser.
class RootHistogramFit {
 public:
  /// Configuration for @c RootHistogramFit
  struct Config {
    /// `TH1::Fit` option string, applied to both @c fit overloads. Must keep
    /// `"S"` (return a `TFitResult`, which the implementation reads) and
    /// `"0"` (do not draw); `"Q"` is strongly recommended to suppress ROOT's
    /// fit printout. The ranged overload adds `"R"` itself.
    ///
    /// Defaults to `"SQ0"`, ROOT's least-squares fit -- the counterpart of
    /// `Acts::Experimental::GaussianHistogramFit::Objective::ChiSquare`. Use
    /// `"LSQ0"` for the likelihood fit, the counterpart of
    /// `Objective::PoissonLikelihood`.
    std::string fitOptions = "SQ0";
  };

  RootHistogramFit() = default;

  /// Construct with the given configuration
  /// @param config The fit configuration
  explicit RootHistogramFit(Config config) : m_config(std::move(config)) {}

  /// The fit configuration
  /// @return The configuration this instance was constructed with
  const Config& config() const { return m_config; }

  /// Fit a Gaussian to a histogram
  ///
  /// @param hist The histogram to fit
  /// @return The fit result, or an error if the fit could not be performed
  Acts::Result<Acts::Experimental::GaussianHistogramFitResult> fit(
      const Acts::Experimental::Histogram1& hist) const;

  /// Fit a Gaussian to the bins of a histogram whose centre lies in a range
  ///
  /// @param hist The histogram to fit
  /// @param xMin Lower end of the fit range (inclusive)
  /// @param xMax Upper end of the fit range (inclusive)
  /// @return The fit result, or an error if the fit could not be performed
  Acts::Result<Acts::Experimental::GaussianHistogramFitResult> fit(
      const Acts::Experimental::Histogram1& hist, double xMin,
      double xMax) const;

 private:
  Config m_config{};
};

}  // namespace ActsPlugins
