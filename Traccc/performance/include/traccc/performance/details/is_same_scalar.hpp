// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Library include(s).
#include "traccc/definitions/common.hpp"
#include "traccc/definitions/primitives.hpp"

namespace traccc::details {

/// Comparison of scalars, used in algorithmic code validation
///
/// This helper function can be used to decide if two scalar values are "the
/// same" within some specific uncertainty.
///
/// @param lhs The Left Hand Side of the comparison
/// @param rhs The Right Hand Side of the comparison
/// @param unc The uncertainty percentage expressed in the 0.0-1.0 range
/// @return @c true if the two values are "the same" within uncertainty,
///         @c false otherwise
///
bool is_same_scalar(scalar lhs, scalar rhs, scalar unc = float_epsilon);

}  // namespace traccc::details
