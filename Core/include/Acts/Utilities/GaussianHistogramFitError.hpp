// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <system_error>
#include <type_traits>

namespace Acts::Experimental {

/// Error codes for a single Gaussian fit to a histogram
/// @ingroup errors
enum class GaussianHistogramFitError {
  // ensure all values are non-zero
  /// No bins in the fit range, or their contents summed to zero or less
  EmptyRange = 1,
  /// Fewer populated bins than the fitter's configured minimum
  TooFewNonEmptyBins,
  /// The initial mean/sigma guess could not be formed
  NoValidSeed,
  /// The minimiser returned a non-finite mean or a non-positive sigma
  NonFiniteParameters,
  /// The underlying fitter reported failure
  FitFailed,
};

/// Create error code from GaussianHistogramFitError
/// @param e The error code enum value
/// @return Standard error code
std::error_code make_error_code(GaussianHistogramFitError e);

}  // namespace Acts::Experimental

namespace std {
// register with STL
template <>
struct is_error_code_enum<Acts::Experimental::GaussianHistogramFitError>
    : std::true_type {};
}  // namespace std
