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

/// Error codes for Eigen stepper operations
/// @ingroup errors
enum class EigenStepperError {
  // ensure all values are non-zero
  /// Step size fell below minimum threshold
  StepSizeStalled = 1,
  /// Step calculation was invalid
  StepInvalid,
  /// Step size adjustment exceeds maximum trials
  StepSizeAdjustmentFailed,
};

/// Create error code from EigenStepperError
/// @param e The error code enum value
/// @return Standard error code
std::error_code make_error_code(EigenStepperError e);

}  // namespace Acts

namespace std {
// register with STL
template <>
struct is_error_code_enum<Acts::EigenStepperError> : std::true_type {};
}  // namespace std
