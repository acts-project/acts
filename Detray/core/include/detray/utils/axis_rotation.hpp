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
#include "detray/definitions/math.hpp"

namespace detray {

/// @brief Helper struct to rotate a vector around a given axis and angle
/// counterclockwisely
template <concepts::algebra algebra_t>
struct axis_rotation {
 public:
  using algebra_type = algebra_t;
  using scalar_type = dscalar<algebra_t>;
  using vector3_type = dvector3D<algebra_t>;
  template <std::size_t ROWS, std::size_t COLS>
  using matrix_type = dmatrix<algebra_t, ROWS, COLS>;

  /// @brief Constructor for axis rotation
  ///
  /// @param axis rotation axis
  /// @param theta rotation angle
  DETRAY_HOST_DEVICE
  axis_rotation(const vector3_type& axis, const scalar_type theta) {
    // normalize the axis
    const vector3_type U = vector::normalize(axis);

    scalar_type cos_theta{math::cos(theta)};

    auto I = matrix::identity<matrix_type<3, 3>>();
    matrix_type<3, 3> axis_cross = matrix::cross_matrix(U);
    matrix_type<3, 3> axis_outer = matrix::outer_product(U, U);

    R = cos_theta * I + math::sin(theta) * axis_cross +
        (1.f - cos_theta) * axis_outer;
  }

  /// @param v vector to be rotated
  /// @returns Get the counterclockwisely-rotated vector
  template <typename vector3_t>
    requires(concepts::vector3D<vector3_t> ||
             concepts::column_matrix3D<vector3_t>)
  DETRAY_HOST_DEVICE vector3_t operator()(const vector3_t& v) const {
    return R * v;
  }

 private:
  /// Rotation matrix
  matrix_type<3, 3> R{matrix::identity<matrix_type<3, 3>>()};
};

/// @brief Helper struct to perform an euler rotation for a given vector
/// All rotation operations are counterclockwise
template <concepts::algebra algebra_t>
struct euler_rotation {
 public:
  using algebra_type = algebra_t;
  using scalar_type = dscalar<algebra_t>;
  using vector3_type = dvector3D<algebra_t>;

  // Following the z-x-z convention
  vector3_type x{1.0f, 0.f, 0.f};
  vector3_type z{0.f, 0.f, 1.f};
  scalar_type alpha{0.f};
  scalar_type beta{0.f};
  scalar_type gamma{0.f};

  /// @returns Get the new x and z axis
  DETRAY_HOST_DEVICE std::pair<vector3_type, vector3_type> operator()() const {
    // alpha around z axis
    axis_rotation<algebra_t> axis_rot_alpha(z, alpha);
    auto new_x = axis_rot_alpha(x);

    // beta around x' axis
    axis_rotation<algebra_t> axis_rot_beta(new_x, beta);
    const auto new_z = axis_rot_beta(z);

    axis_rotation<algebra_t> axis_rot_gamma(new_z, gamma);
    new_x = axis_rot_gamma(new_x);

    return std::make_pair(new_x, new_z);
  }

  /// @param v vector to be rotated
  /// @returns Get the rotated vector
  DETRAY_HOST_DEVICE vector3_type operator()(const vector3_type& v0) const {
    // alpha around z axis
    axis_rotation<algebra_t> axis_rot_alpha(z, alpha);
    const auto new_x = axis_rot_alpha(x);
    const auto v1 = axis_rot_alpha(v0);

    // beta around x' axis
    axis_rotation<algebra_t> axis_rot_beta(new_x, beta);
    const auto new_z = axis_rot_beta(z);
    const auto v2 = axis_rot_beta(v1);

    // gamma around z' axis
    axis_rotation<algebra_t> axis_rot_gamma(new_z, gamma);
    const auto v3 = axis_rot_gamma(v2);

    return v3;
  }
};

}  // namespace detray
