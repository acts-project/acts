// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Histogram.hpp"
#include "Acts/Utilities/HistogramFit.hpp"

#include <optional>
#include <string>
#include <utility>

namespace ActsPlugins {

/// Fit a Gaussian to a histogram via ROOT's `TH1::Fit`
///
/// `operator()`'s signature matches `Acts::Experimental::HistogramFitFunction`
/// exactly, so a `RootHistogramFit` instance can be used directly as one with
/// no adapter.
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
    /// `ActsExamples::gaussianHistogramFit`. Use `"LSQ0"` for the likelihood
    /// fit instead.
    std::string fitOptions = "SQ0";
  };

  /// @c Acts::Experimental::HistogramFitResult, re-exported for convenience
  using Result = Acts::Experimental::HistogramFitResult;
  /// @c Acts::Experimental::HistogramFitRange, re-exported for convenience
  using Range = Acts::Experimental::HistogramFitRange;

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
