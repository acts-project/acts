// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <system_error>
#include <type_traits>

namespace ActsFatras {
namespace detail {

enum class SimulationError {
  // ensure all values are non-zero
  eInvalidInputParticleId = 1,
};

/// Construct and error_code from the enum.
///
/// Must use snake_case naming for STL compatibility.
std::error_code make_error_code(SimulationError e);

}  // namespace detail
}  // namespace ActsFatras

// Register the error enum as STL-compatible.
namespace std {
template <>
struct is_error_code_enum<ActsFatras::detail::SimulationError>
    : std::true_type {};
}  // namespace std
