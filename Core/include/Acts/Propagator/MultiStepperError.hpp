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

namespace Acts {

/// Error codes for multi-stepper operations
/// @ingroup errors
enum class MultiStepperError {
  // ensure all values are non-zero
  /// Component is not on a surface
  ComponentNotOnSurface = 1,
  /// The global BoundState/CurvilinearState can only be computed if only one
  /// component exists
  StateOfMultipleComponentsRequested = 2,
  /// The average track has left the current volume
  AverageTrackLeftCurrentVolume = 3,
  /// Stepping error occurred in all components
  AllComponentsSteppingError = 4,
  /// The conversion to the bound state failed for all components
  AllComponentsConversionToBoundFailed = 5,
  /// The conversion to the bound state failed for some components
  SomeComponentsConversionToBoundFailed = 6
};

/// Create error code from MultiStepperError
/// @param e The error code enum value
/// @return Standard error code
std::error_code make_error_code(Acts::MultiStepperError e);

}  // namespace Acts

// register with STL
namespace std {
template <>
struct is_error_code_enum<Acts::MultiStepperError> : std::true_type {};
}  // namespace std
