// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Histogram.hpp"

#include <optional>
#include <string>
#include <tuple>
#include <utility>

namespace ActsPlugins {

/// Fit a Gaussian to a histogram via ROOT's `TH1::Fit`
///
/// `operator()`'s signature matches `ActsExamples::HistogramFitFunction`
/// exactly, so a `RootHistogramFit` instance can be used directly as one with
/// no adapter.
class RootHistogramFit {
 public:
  /// Configuration for @c RootHistogramFit
  struct Config {
    /// `TH1::Fit` option string. Must keep `"S"` (return a `TFitResult`) and
    /// `"0"` (do not draw); the ranged fit adds `"R"` itself. Defaults to
    /// ROOT's least-squares fit; use `"LSQ0"` for the likelihood fit instead.
    std::string fitOptions = "SQ0";
  };

  /// Outcome of the fit: `(mean, sigma, meanError, sigmaError)`
  using Result = std::tuple<double, double, double, double>;
  /// Fit range `[xMin, xMax]`, closed, selected by bin centre
  using Range = std::pair<double, double>;

  RootHistogramFit() = default;

  /// Construct with the given configuration
  /// @param config The fit configuration
  explicit RootHistogramFit(Config config) : m_config(std::move(config)) {}

  /// The fit configuration
  /// @return The configuration this instance was constructed with
  const Config& config() const { return m_config; }

  /// Fit a Gaussian to a histogram, optionally restricted to a range
  ///
  /// @param hist The histogram to fit
  /// @param range If set, only bins whose centre lies in `[range->first,
  ///              range->second]` enter the fit
  /// @return `(mean, sigma, meanError, sigmaError)`, or `std::nullopt` if the
  ///         fit could not be performed
  std::optional<Result> operator()(
      const Acts::Experimental::Histogram1& hist,
      std::optional<Range> range = std::nullopt) const;

 private:
  Config m_config{};
};

}  // namespace ActsPlugins
