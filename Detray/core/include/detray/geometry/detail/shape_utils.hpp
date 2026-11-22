// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/containers.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/units.hpp"

// System include(s)
#include <cstdint>

namespace detray {

namespace detail {

/// Generate phi tolerance from distance tolerance
///
/// @param tol is the distance tolerance in mm
/// @param radius is the radius of the shape
///
/// @return the opening angle of a chord the size of tol (= 2*arcsin(c/(2r)))
/// using a small angle approximation
template <concepts::scalar scalar_t>
constexpr scalar_t phi_tolerance(scalar_t tol, scalar_t radius) {
  return radius > 0.f ? tol / radius : tol;
}
/// Add vertices on a full circle around the local z-axis
///
/// @tparam algebra_t the algebra type
///
/// @param radius the radius of the circle
/// @param z the position of the circle on the local z-axis
/// @param n_seg the number of segments (at least 4 are used)
/// @param vertices the vertices in local cartesian coordinates (output)
template <concepts::algebra algebra_t>
DETRAY_HOST inline void add_circle_vertices(
    const dscalar<algebra_t> radius, const dscalar<algebra_t> z,
    const dindex n_seg, dvector<dpoint3D<algebra_t>>& vertices) {
  using scalar_t = dscalar<algebra_t>;

  const dindex n{n_seg < 4u ? 4u : n_seg};
  const scalar_t step{2.f * constant<scalar_t>::pi / static_cast<scalar_t>(n)};
  for (dindex i = 0u; i < n; ++i) {
    const scalar_t phi{static_cast<scalar_t>(i) * step};
    vertices.push_back({radius * math::cos(phi), radius * math::sin(phi), z});
  }
}

// Result of an 'inside' check: First entry precise check, second entry with
// tolerance
template <typename bool_t>
using boundary_check_result = darray<bool_t, 2>;

}  // namespace detail

/// Address different types of check results
enum class check_type : std::uint_least8_t {
  e_precise = 0u,
  e_with_edge = 1u,
};

/// Get the result of a check with edge tolerance
template <check_type C, typename bool_t>
DETRAY_HOST_DEVICE constexpr bool_t get(
    const detail::boundary_check_result<bool_t> result) {
  if constexpr (C == check_type::e_precise) {
    return result[0];
  } else if constexpr (C == check_type::e_with_edge) {
    return result[1];
  } else {
    // Broadcast
    return bool_t(false);
  }
}

/// Get the result of a simple tolerance check
template <check_type C, typename bool_t>
DETRAY_HOST_DEVICE constexpr bool_t get(const bool_t& result) {
  return result;
}

}  // namespace detray
