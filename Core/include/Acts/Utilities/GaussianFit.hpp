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

#include <optional>

namespace Acts::Experimental {

/// @brief Outcome of a Gaussian fit to a 1D histogram
struct GaussianFitResult {
  /// Fitted mean of the Gaussian
  double mean{};
  /// Fitted standard deviation of the Gaussian
  double sigma{};
  /// Estimated uncertainty on @c mean
  double meanError{};
  /// Estimated uncertainty on @c sigma
  double sigmaError{};
};

/// Fit a Gaussian to a 1D histogram by maximising the Poisson likelihood
///
/// The model is @f$ \mu_i = A \exp(-(x_i - m)^2 / (2 s^2)) @f$ with @f$ x_i @f$
/// the bin centre, compared directly against the bin content. This mirrors
/// ROOT's predefined `"gaus"` under `TH1::Fit(..., "L")`: no bin-width factor
/// and no integration of the model over the bin.
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
///
/// @param hist The histogram to fit
/// @return The fit result, or `std::nullopt` if the fit could not be performed
/// @note Under- and overflow bins are ignored.
std::optional<GaussianFitResult> gaussianFit(const Histogram1& hist);

/// Fit a Gaussian to the bins of a 1D histogram whose centre lies in a range
///
/// @param hist The histogram to fit
/// @param xMin Lower end of the fit range (inclusive)
/// @param xMax Upper end of the fit range (inclusive)
/// @return The fit result, or `std::nullopt` if the fit could not be performed
/// @remark Bin selection is by bin centre on a closed interval, matching how
///         ROOT restricts a fit range.
/// @note See @c gaussianFit(const Histogram1&) for the fitted model.
std::optional<GaussianFitResult> gaussianFit(const Histogram1& hist,
                                             double xMin, double xMax);

/// Fit a Gaussian repeatedly, narrowing the fit range around the peak
///
/// The first fit uses the full histogram. Each subsequent fit is restricted to
/// @f$ m \pm \mathrm{sigmaRange} \cdot s @f$ using the previous iteration's
/// parameters. The returned uncertainties come from the final iteration.
///
/// @param hist The histogram to fit
/// @param sigmaRange Half-width of the restricted range, in fitted sigmas
/// @param iterations Total number of fits, including the initial unrestricted
///                   one; values below 1 are treated as 1
/// @param logger Logger for diagnostics on failed iterations
/// @return The fit result, or `std::nullopt` if any iteration failed
std::optional<GaussianFitResult> iterativeGaussianFit(
    const Histogram1& hist, double sigmaRange, int iterations,
    const Logger& logger = getDummyLogger());

}  // namespace Acts::Experimental
