// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/utils/concepts.hpp"

namespace detray {

/// Components of a bound track parameters vector.
///
enum bound_indices : unsigned int {
  // Local position on the reference surface.
  // This is intentionally named different from the position components in
  // the other data vectors, to clarify that this is defined on a surface
  // while the others are defined in free space.
  e_bound_loc0 = 0u,
  e_bound_loc1 = 1u,
  // Direction angles
  e_bound_phi = 2u,
  e_bound_theta = 3u,
  // Global inverse-momentum-like parameter, i.e. q/p or 1/p
  // The naming is inconsistent for the case of neutral track parameters where
  // the value is interpreted as 1/p not as q/p. This is intentional to avoid
  // having multiple aliases for the same element and for lack of an
  // acceptable
  // common name.
  e_bound_qoverp = 4u,
  e_bound_time = 5u,
  // Last uninitialized value contains the total number of components
  e_bound_size,
};

/// Components of a free track parameters vector.
///
/// To be used to access components by named indices instead of just numbers.
/// This must be a regular `enum` and not a scoped `enum class` to allow
/// implicit conversion to an integer. The enum value are thus visible directly
/// in `namespace Acts` and are prefixed to avoid naming collisions.
enum free_indices : unsigned int {
  // Spatial position
  // The spatial position components must be stored as one continuous block.
  e_free_pos0 = 0u,
  e_free_pos1 = e_free_pos0 + 1u,
  e_free_pos2 = e_free_pos0 + 2u,
  // Time
  e_free_time = 3u,
  // (Unit) direction
  // The direction components must be stored as one continuous block.
  e_free_dir0 = 4u,
  e_free_dir1 = e_free_dir0 + 1u,
  e_free_dir2 = e_free_dir0 + 2u,
  // Global inverse-momentum-like parameter, i.e. q/p or 1/p
  // See BoundIndices for further information
  e_free_qoverp = 7u,
  // Last uninitialized value contains the total number of components
  e_free_size,
};

/// Vector type for free track parametrization
template <concepts::algebra algebra_t>
using free_vector = dmatrix<algebra_t, e_free_size, 1>;

/// Covariance matrix type for free track parametrization
template <concepts::algebra algebra_t>
using free_matrix = dmatrix<algebra_t, e_free_size, e_free_size>;

/// Vector type for bound track parametrization
template <concepts::algebra algebra_t>
using bound_vector = dmatrix<algebra_t, e_bound_size, 1>;

/// Covariance matrix type for bound track parametrization
template <concepts::algebra algebra_t>
using bound_matrix = dmatrix<algebra_t, e_bound_size, e_bound_size>;

/// Mapping from bound to free track parameters.
template <concepts::algebra algebra_t>
using bound_to_free_matrix = dmatrix<algebra_t, e_free_size, e_bound_size>;

/// Submatrix for bound to free Jacobians for d(x,y,z) or d(nx, ny, nz) over
/// d(l0, l1) or d(phi, theta)
template <concepts::algebra algebra_t>
using bound_to_free_jacobian_submatrix = dmatrix<algebra_t, 3, 2>;

/// Mapping from free to bound track parameters.
template <concepts::algebra algebra_t>
using free_to_bound_matrix = dmatrix<algebra_t, e_bound_size, e_free_size>;

/// Submatrix for free to bound Jacobians d(l0, l1) or d(phi, theta) over
/// d(x,y,z) or d(nx, ny, nz)
template <concepts::algebra algebra_t>
using free_to_bound_jacobian_submatrix = dmatrix<algebra_t, 2, 3>;

/// Mapping from free to path
template <concepts::algebra algebra_t>
using free_to_path_matrix = dmatrix<algebra_t, 1, e_free_size>;

/// Mapping from path to free
template <concepts::algebra algebra_t>
using path_to_free_matrix = dmatrix<algebra_t, e_free_size, 1>;

}  // namespace detray
