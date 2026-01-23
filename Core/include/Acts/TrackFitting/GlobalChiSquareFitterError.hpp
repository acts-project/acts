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

/// @enum GlobalChiSquareFitterError
/// @ingroup errors
enum class GlobalChiSquareFitterError {
  // ensure all values are non-zero
  /// aMatrix is not invertible.
  AIsNotInvertible = 1,
  /// Did not converge in 'nUpdateMax' updates.
  DidNotConverge = 2,
  /// Not enough measurements.
  NotEnoughMeasurements = 3,
  /// Update pushed the parameters to a new volume.
  UpdatePushedToNewVolume = 4,

};

std::error_code make_error_code(
    Acts::Experimental::GlobalChiSquareFitterError e);

}  // namespace Acts::Experimental

namespace std {
// register with STL
template <>
struct is_error_code_enum<Acts::Experimental::GlobalChiSquareFitterError>
    : std::true_type {};
}  // namespace std
