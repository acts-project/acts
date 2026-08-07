// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Histogram.hpp"

#include <functional>
#include <optional>
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
/// Any backend -- @c ActsExamples::gaussianHistogramFit (a chi-square fit on
/// a boost-histogram), @c ActsPlugins::RootHistogramFit (ROOT's `TH1::Fit`),
/// or a Python callable -- can be adapted to this signature.
/// `ActsPlugins::RootHistogramFit::fit` matches it exactly (both
/// its element types are plain `std::tuple`/`std::pair`, not shared with
/// this header, since the ROOT plugin must not depend on Examples), so it
/// needs no adapter beyond a capturing lambda for the member-function call.
using HistogramFitFunction = std::function<std::optional<HistogramFitResult>(
    const Acts::Experimental::Histogram1&, std::optional<HistogramFitRange>)>;

}  // namespace ActsExamples
