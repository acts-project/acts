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

/// Recoverable errors in material covariance updates.
/// @ingroup errors
enum class MaterialInteractionError {
  /// Requested noise removal would produce an invalid covariance.
  InvalidCovarianceAfterNoiseRemoval = 1,
};

/// Convert a material-interaction error to a standard error code.
/// @param error Error enum value
/// @return Standard error code
std::error_code make_error_code(MaterialInteractionError error);

}  // namespace Acts

namespace std {
template <>
struct is_error_code_enum<Acts::MaterialInteractionError> : std::true_type {};
}  // namespace std
