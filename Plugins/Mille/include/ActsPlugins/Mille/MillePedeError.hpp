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

namespace ActsPlugins::ActsToMille {
/// Error codes for Millepede run
/// @ingroup errors
enum class MillePedeError {
  InstallationNotFound = 1,  // no valid install found
  SteeringNotFound = 2,      // steering file not found
  SolverCrash = 3,           // solver crashed
  InvalidSolution = 4,       // solver finished but the solution is invalid
  SolutionNotReadable = 5,   // solution file could not be read
};

/// @cond
/// Create error code from MillePedeError
/// @param e The error code enum value
/// @return Standard error code
std::error_code make_error_code(ActsToMille::MillePedeError e);
/// @endcond

}  // namespace ActsPlugins::ActsToMille

namespace std {
// register with STL
template <>
struct is_error_code_enum<ActsPlugins::ActsToMille::MillePedeError>
    : std::true_type {};
}  // namespace std
