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

/// Errors specific to the AtlasStepper coordinate model.
/// @ingroup errors
enum class AtlasStepperError {
  /// The stepper's hand-written derivatives require global-spherical
  /// coordinates.
  UnsupportedBoundParametrization = 1,
};

/// Convert an AtlasStepper error to a standard error code.
/// @param error Error enum value
/// @return Standard error code
std::error_code make_error_code(AtlasStepperError error);

}  // namespace Acts

namespace std {
template <>
struct is_error_code_enum<Acts::AtlasStepperError> : std::true_type {};
}  // namespace std
