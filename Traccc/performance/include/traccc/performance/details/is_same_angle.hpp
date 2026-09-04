/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Library include(s).
#include "traccc/definitions/common.hpp"
#include "traccc/definitions/primitives.hpp"

namespace traccc::details {

/// Comparison of angles, used in algorithmic code validation
///
/// This helper function can be used to decide if two angle values are "the
/// same" within some specific uncertainty.
///
/// It is different from @c traccc::details:is_same_scalar in that it considers
/// values just above -Pi and just below Pi to be close to each other.
///
/// @param lhs The Left Hand Side of the comparison
/// @param rhs The Right Hand Side of the comparison
/// @param unc The uncertainty percentage expressed in the 0.0-1.0 range
/// @return @c true if the two values are "the same" within uncertainty,
///         @c false otherwise
///
bool is_same_angle(scalar lhs, scalar rhs, scalar unc = float_epsilon);

}  // namespace traccc::details
