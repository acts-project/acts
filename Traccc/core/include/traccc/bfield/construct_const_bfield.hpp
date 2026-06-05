/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/bfield/magnetic_field.hpp"
#include "traccc/definitions/primitives.hpp"

namespace traccc {

/// Construct a constant magnetic field object
///
/// @tparam scalar_t The scalar type to construct the field with
///
/// @param x The X component of the constant field
/// @param y The Y component of the constant field
/// @param z The Z component of the constant field
///
template <typename scalar_t>
magnetic_field construct_const_bfield(scalar_t x, scalar_t y, scalar_t z);

/// Construct a constant magnetic field object
///
/// @param v The 3-vector describing the constant field
///
magnetic_field construct_const_bfield(const vector3& v);

}  // namespace traccc

// Include the implementation.
#include "traccc/bfield/impl/construct_const_bfield.ipp"
