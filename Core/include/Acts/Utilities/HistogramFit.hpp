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

namespace Acts::Experimental {

/// @brief Outcome of a Gaussian fit to a 1D histogram: `(mean, sigma,
///        meanError, sigmaError)`
///
/// A plain @c std::tuple rather than a named struct so that a Python fit
/// backend can return an ordinary 4-tuple and have it convert automatically
/// via `pybind11/stl.h`, with no dedicated binding required.
using HistogramFitResult = std::tuple<double, double, double, double>;

/// @brief Fit range `[xMin, xMax]`, closed, selected by bin centre
using HistogramFitRange = std::pair<double, double>;

/// A single Gaussian fit to a 1D histogram, optionally restricted to a range
///
/// Any backend -- @c ActsExamples::gaussianHistogramFit, a callable
/// `ActsPlugins::RootHistogramFit`, or a Python callable -- can be adapted to
/// this signature. Living in Core rather than Examples or a plugin lets both
/// sides share the same vocabulary types without depending on each other.
using HistogramFitFunction = std::function<std::optional<HistogramFitResult>(
    const Histogram1&, std::optional<HistogramFitRange>)>;

}  // namespace Acts::Experimental
