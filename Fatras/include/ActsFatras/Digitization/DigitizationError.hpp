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

namespace ActsFatras {

/// Error codes for digitization operations
/// @ingroup errors
enum class DigitizationError {
  // ensure all values are non-zero
  /// Smeared out of surface bounds.
  SmearingOutOfRange = 1,
  /// Smearing error occurred.
  SmearingError,
  /// Surface undefined for this operation.
  UndefinedSurface,
  /// Surface mask could not be applied.
  MaskingError,
  /// Maximum number of retries exceeded.
  MaximumRetriesExceeded,
};

/// Create error code from DigitizationError
/// @param e Digitization error enum value
/// @return Error code corresponding to the error
std::error_code make_error_code(DigitizationError e);

}  // namespace ActsFatras

namespace std {
// register with STL
template <>
struct is_error_code_enum<ActsFatras::DigitizationError> : std::true_type {};
}  // namespace std
