// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/track_parametrization.hpp"
#include "detray/geometry/coordinates/concentric_cylindrical2D.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes/concentric_cylinder2D.hpp"
#include "detray/propagator/detail/jacobian.hpp"

namespace detray::detail {

/// @brief Specialization for 2D concentric cylindrical frames
template <concepts::algebra algebra_t>
struct jacobian<concentric_cylindrical2D<algebra_t>> {
  /// @name Type definitions for the struct
  /// @{
  using coordinate_frame = concentric_cylindrical2D<algebra_t>;

  using algebra_type = algebra_t;
  using transform3_type = dtransform3D<algebra_t>;
  using scalar_type = dscalar<algebra_t>;
  using point2_type = dpoint2D<algebra_t>;
  using point3_type = dpoint3D<algebra_t>;
  using vector3_type = dvector3D<algebra_t>;

  using bound_to_free_matrix_type = bound_to_free_matrix<algebra_t>;
  using free_to_bound_matrix_type = free_to_bound_matrix<algebra_t>;
  using free_to_path_matrix_type = free_to_path_matrix<algebra_t>;
  using bound_to_free_jacobian_submatrix_type =
      bound_to_free_jacobian_submatrix<algebra_type>;
  using free_to_bound_jacobian_submatrix_type =
      free_to_bound_jacobian_submatrix<algebra_type>;

  DETRAY_HOST_DEVICE static constexpr free_to_path_matrix_type path_derivative(
      const transform3_type & /*trf*/, const point3_type &pos,
      const vector3_type &dir, const vector3_type & /*dtds*/) {
    using mask_t = mask<concentric_cylinder2D, algebra_type, std::uint8_t>;
    constexpr mask_t dummy_mask{};

    const transform3_type identity{};

    const point2_type local =
        coordinate_frame::global_to_local(identity, pos, {});
    const vector3_type normal =
        coordinate_frame::normal(identity, local, dummy_mask);

    const vector3_type pos_term = (-1.f / vector::dot(normal, dir)) * normal;

    auto derivative{matrix::zero<free_to_path_matrix_type>()};
    getter::element(derivative, 0u, e_free_pos0) = pos_term[0];
    getter::element(derivative, 0u, e_free_pos1) = pos_term[1];
    getter::element(derivative, 0u, e_free_pos2) = pos_term[2];

    return derivative;
  }

  DETRAY_HOST_DEVICE
  static constexpr bound_to_free_jacobian_submatrix_type
  get_derivative_dpos_dloc(const transform3_type & /*trf*/,
                           const point3_type &pos,
                           const vector3_type & /*dir*/) {
    bound_to_free_jacobian_submatrix_type bound_pos_to_free_pos_derivative =
        matrix::zero<bound_to_free_jacobian_submatrix_type>();

    const scalar_type r{vector::perp(pos)};
    const scalar_type phi{vector::phi(pos)};

    // Assert the integrity of the local matrix wrt to the globally
    // defined track parameterization.
    static_assert(e_bound_loc1 == e_bound_loc0 + 1);
    static_assert(e_free_pos1 == e_free_pos0 + 1);
    static_assert(e_free_pos2 == e_free_pos0 + 2);

    constexpr unsigned int e_submatrix_bound_loc0 = 0u;
    constexpr unsigned int e_submatrix_bound_loc1 = 1u;
    constexpr unsigned int e_submatrix_free_pos0 = 0u;
    constexpr unsigned int e_submatrix_free_pos1 = 1u;
    constexpr unsigned int e_submatrix_free_pos2 = 2u;

    getter::element(bound_pos_to_free_pos_derivative, e_submatrix_free_pos0,
                    e_submatrix_bound_loc0) = -r * math::sin(phi);
    getter::element(bound_pos_to_free_pos_derivative, e_submatrix_free_pos1,
                    e_submatrix_bound_loc0) = r * math::cos(phi);
    getter::element(bound_pos_to_free_pos_derivative, e_submatrix_free_pos2,
                    e_submatrix_bound_loc1) = 1.f;

    return bound_pos_to_free_pos_derivative;
  }

  DETRAY_HOST_DEVICE
  static constexpr bound_to_free_jacobian_submatrix_type
  get_derivative_dpos_dangle(
      const transform3_type & /*unused*/, const point3_type & /*unused*/,
      const vector3_type & /*unused*/,
      const bound_to_free_jacobian_submatrix_type & /*unused*/) {
    return matrix::zero<bound_to_free_jacobian_submatrix_type>();
  }

  DETRAY_HOST_DEVICE
  static constexpr free_to_bound_jacobian_submatrix_type
  get_derivative_dloc_dpos(const transform3_type & /*unused*/,
                           const point3_type &pos,
                           const vector3_type & /*unused*/) {
    // Get d(loc0, loc1)/d(x,y,z)
    free_to_bound_jacobian_submatrix_type free_pos_to_bound_pos_derivative =
        matrix::zero<free_to_bound_jacobian_submatrix_type>();

    const scalar_type r_inv{1.f / vector::perp(pos)};
    const scalar_type phi{vector::phi(pos)};

    // Assert the integrity of the local matrix wrt to the globally
    // defined track parameterization.
    static_assert(e_bound_loc1 == e_bound_loc0 + 1);
    static_assert(e_free_pos1 == e_free_pos0 + 1);
    static_assert(e_free_pos2 == e_free_pos0 + 2);

    constexpr unsigned int e_submatrix_bound_loc0 = 0u;
    constexpr unsigned int e_submatrix_bound_loc1 = 1u;
    constexpr unsigned int e_submatrix_free_pos0 = 0u;
    constexpr unsigned int e_submatrix_free_pos1 = 1u;
    constexpr unsigned int e_submatrix_free_pos2 = 2u;

    getter::element(free_pos_to_bound_pos_derivative, e_submatrix_bound_loc0,
                    e_submatrix_free_pos0) = -r_inv * math::sin(phi);
    getter::element(free_pos_to_bound_pos_derivative, e_submatrix_bound_loc0,
                    e_submatrix_free_pos1) = r_inv * math::cos(phi);
    getter::element(free_pos_to_bound_pos_derivative, e_submatrix_bound_loc1,
                    e_submatrix_free_pos2) = 1.f;

    return free_pos_to_bound_pos_derivative;
  }
};

}  // namespace detray::detail
