// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Histogram.hpp"
#include "ActsExamples/Validation/HistogramFit.hpp"

#include <optional>

namespace ActsExamples {

/// Fit a Gaussian to a 1D histogram by chi-square minimisation
///
/// The model is @f$ \mu_i = A \exp(-(x_i - m)^2 / (2 s^2)) @f$, @f$ x_i @f$
/// the bin centre, compared directly against the bin content: no bin-width
/// factor and no integration of the model over the bin. Minimises
/// @f$ \sum_i (n_i - \mu_i)^2 / n_i @f$ over bins with @f$ n_i > 0 @f$ --
/// Neyman's chi-square with Poisson counting errors, matching ROOT's
/// predefined `"gaus"` under `TH1::Fit(..., "SQ0")`, which gives zero-content
/// bins zero error and drops them from the sum. See
/// @c ActsPlugins::RootHistogramFit for a ROOT-backed fit with the same
/// interface.
///
/// @param hist The histogram to fit
/// @param range If set, only bins whose centre lies in `[range->first,
///              range->second]` enter the fit; matches how ROOT restricts a
///              fit range. If unset, the fit uses every bin.
/// @return `(mean, sigma, meanError, sigmaError)`, or @c std::nullopt if the
///         fit could not be performed
/// @note Under- and overflow bins are always ignored.
std::optional<HistogramFitResult> gaussianHistogramFit(
    const Acts::Experimental::Histogram1& hist,
    std::optional<HistogramFitRange> range = std::nullopt);

}  // namespace ActsExamples
