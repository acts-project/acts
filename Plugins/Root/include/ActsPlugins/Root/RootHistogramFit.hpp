// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/GaussianFit.hpp"
#include "Acts/Utilities/Histogram.hpp"
#include "Acts/Utilities/Result.hpp"

namespace ActsPlugins {

/// Fit a Gaussian to a histogram via ROOT's `TH1::Fit`
///
/// Satisfies `Acts::Experimental::GaussianFitter`, so it can be used with
/// `Acts::Experimental::iterativeFit` and
/// `Acts::Experimental::extractMeanWidthProfiles` in place of
/// `Acts::Experimental::GaussianHistogramFit`, backed by ROOT's own
/// least-squares `"gaus"` fitter instead of Core's Poisson-likelihood
/// optimiser.
class RootHistogramFit {
 public:
  RootHistogramFit() = default;

  /// Fit a Gaussian to a histogram
  ///
  /// @param hist The histogram to fit
  /// @return The fit result, or an error if the fit could not be performed
  Acts::Result<Acts::Experimental::GaussianFitResult> fit(
      const Acts::Experimental::Histogram1& hist) const;

  /// Fit a Gaussian to the bins of a histogram whose centre lies in a range
  ///
  /// @param hist The histogram to fit
  /// @param xMin Lower end of the fit range (inclusive)
  /// @param xMax Upper end of the fit range (inclusive)
  /// @return The fit result, or an error if the fit could not be performed
  Acts::Result<Acts::Experimental::GaussianFitResult> fit(
      const Acts::Experimental::Histogram1& hist, double xMin,
      double xMax) const;
};

}  // namespace ActsPlugins
