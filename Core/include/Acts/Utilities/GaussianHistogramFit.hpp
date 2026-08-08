// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Histogram.hpp"
#include "Acts/Utilities/Result.hpp"

namespace Acts::Experimental {

/// @brief Outcome of a Gaussian fit to a 1D histogram
struct GaussianHistogramFitResult {
  /// Fitted mean of the Gaussian
  double mean{};
  /// Fitted standard deviation of the Gaussian
  double sigma{};
  /// Estimated uncertainty on @c mean
  double meanError{};
  /// Estimated uncertainty on @c sigma
  double sigmaError{};
};

/// Fit a Gaussian to a 1D histogram, with the model
/// @f$ \mu_i = A \exp(-(x_i - m)^2 / (2 s^2)) @f$, @f$ x_i @f$ the bin centre,
/// compared directly against the bin content: no bin-width factor and no
/// integration of the model over the bin. See @c ActsPlugin::RootHistogramFit
/// for a ROOT-backed fit with the same interface.
///
/// Two objectives are available, selected via @c Config::objective:
///
/// - @c Objective::ChiSquare (the default) mirrors ROOT's predefined
///   `"gaus"` under `TH1::Fit(..., "SQ0")`: a least-squares fit with bin
///   variance @f$ n_i @f$, i.e. Poisson counting errors. Bins with
///   @f$ n_i = 0 @f$ have zero error in this convention and are dropped from
///   the sum entirely, so empty bins carry no information.
/// - @c Objective::PoissonLikelihood mirrors `TH1::Fit(..., "LSQ0")`: a
///   Poisson maximum-likelihood fit. Empty bins do carry information here,
///   which makes it the better choice for sparse histograms and histograms
///   with outliers.
///
/// In both cases the amplitude @f$ A @f$ is not searched over; for fixed
/// @f$ (m, s) @f$ it has a closed form (see the objective-specific comments
/// in the implementation), which reduces the search to two parameters and
/// makes it markedly more robust.
class GaussianHistogramFit {
 public:
  /// The objective minimised by @c fit
  enum class Objective {
    /// Least squares with Poisson counting errors; drops empty bins
    ChiSquare,
    /// Poisson maximum likelihood; empty bins are informative
    PoissonLikelihood,
  };

  /// Configuration for @c GaussianHistogramFit
  struct Config {
    /// The objective to minimise
    Objective objective = Objective::ChiSquare;
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
  Result<GaussianHistogramFitResult> fit(const Histogram1& hist) const;

  /// Fit a Gaussian to the bins of a 1D histogram whose centre lies in a range
  ///
  /// @param hist The histogram to fit
  /// @param xMin Lower end of the fit range (inclusive)
  /// @param xMax Upper end of the fit range (inclusive)
  /// @return The fit result, or an error if the fit could not be performed
  /// @remark Bin selection is by bin centre on a closed interval, matching how
  ///         ROOT restricts a fit range.
  /// @note See @c fit(const Histogram1&) for the fitted model.
  Result<GaussianHistogramFitResult> fit(const Histogram1& hist, double xMin,
                                         double xMax) const;

 private:
  Config m_config{};
};

}  // namespace Acts::Experimental
