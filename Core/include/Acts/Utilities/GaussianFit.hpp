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
#include "Acts/Utilities/detail/IterativeGaussianFit.hpp"

#include <array>
#include <optional>
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

/// Configuration for @c GaussianHistogramFit
struct GaussianHistogramFitConfig {
  /// Minimum number of populated bins entering a single fit. The underlying
  /// model has three parameters (amplitude, mean, sigma), so fewer bins
  /// cannot constrain it.
  std::size_t minNonEmptyBins = 3;
  /// Fractional step, relative to the seed sigma, used both to size the
  /// initial Nelder-Mead simplex and to size the finite-difference step for
  /// the parameter-error Hessian. Kept as a single value so the two stay
  /// consistent with each other.
  double relativeStep = 0.1;
};

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
  /// Construct with the given configuration
  /// @param config The fit configuration
  explicit GaussianHistogramFit(GaussianHistogramFitConfig config = {})
      : m_config(config) {}

  /// The fit configuration
  /// @return The configuration this instance was constructed with
  const GaussianHistogramFitConfig& config() const { return m_config; }

  /// Fit a Gaussian to a 1D histogram
  ///
  /// @param hist The histogram to fit
  /// @return The fit result, or `std::nullopt` if the fit could not be performed
  /// @note Under- and overflow bins are ignored.
  std::optional<GaussianFitResult> fit(const Histogram1& hist) const;

  /// Fit a Gaussian to the bins of a 1D histogram whose centre lies in a range
  ///
  /// @param hist The histogram to fit
  /// @param xMin Lower end of the fit range (inclusive)
  /// @param xMax Upper end of the fit range (inclusive)
  /// @return The fit result, or `std::nullopt` if the fit could not be performed
  /// @remark Bin selection is by bin centre on a closed interval, matching how
  ///         ROOT restricts a fit range.
  /// @note See @c fit(const Histogram1&) for the fitted model.
  std::optional<GaussianFitResult> fit(const Histogram1& hist, double xMin,
                                       double xMax) const;

  /// Fit a Gaussian repeatedly, narrowing the fit range around the peak
  ///
  /// @param hist The histogram to fit
  /// @param sigmaRange Half-width of the restricted range, in fitted sigmas
  /// @param iterations Total number of fits, including the initial unrestricted
  ///                   one; values below 1 are treated as 1
  /// @param logger Logger for diagnostics on failed iterations
  /// @return The fit result, or `std::nullopt` if any iteration failed
  std::optional<GaussianFitResult> iterativeFit(
      const Histogram1& hist, double sigmaRange, int iterations,
      const Logger& logger = getDummyLogger()) const;

  /// Fit a Gaussian to every slice of a histogram along its last axis
  ///
  /// For each bin of the outer axes (every axis but the last), the distribution
  /// along the last axis is fitted with @c iterativeFit and the resulting mean
  /// and sigma are stored, with their uncertainties, in the corresponding
  /// output bin.
  ///
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
      const Histogram<Dim>& hist, const std::string& meanName,
      const std::string& widthName, int minEntriesForFit = 5,
      double sigmaRange = 3.0, int iterations = 3,
      const Logger& logger = getDummyLogger()) const {
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
        ValueHistogram<OuterDim>(widthName, hist.title() + " width", axes),
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

    const auto trySlice =
        [&](int flat) -> std::optional<std::optional<GaussianFitResult>> {
      const std::array<int, OuterDim> outerBins = unravel(flat);
      const Histogram1 slice = detail::sliceLastAxis<Dim>(hist, outerBins);
      if (totalContent(slice) < minEntriesForFit) {
        return std::nullopt;
      }
      return iterativeFit(slice, sigmaRange, iterations, logger);
    };

    const auto storeResult = [&](int flat, const GaussianFitResult& result) {
      const std::array<int, OuterDim> outerBins = unravel(flat);
      profiles.mean.setBin(outerBins, result.mean, result.meanError);
      profiles.width.setBin(outerBins, result.sigma, result.sigmaError);
    };

    profiles.fitFailureFraction = Acts::detail::extractMeanWidthProfilesImpl(
        totalOuterBins, trySlice, storeResult);

    return profiles;
  }

 private:
  GaussianHistogramFitConfig m_config;
};

}  // namespace Acts::Experimental
