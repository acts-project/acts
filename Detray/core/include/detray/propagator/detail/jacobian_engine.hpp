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
#include "detray/definitions/track_parametrization.hpp"
#include "detray/propagator/detail/jacobian.hpp"
#include "detray/propagator/detail/jacobian_cartesian.hpp"
#include "detray/propagator/detail/jacobian_concentric_cylindrical.hpp"
#include "detray/propagator/detail/jacobian_cylindrical.hpp"
#include "detray/propagator/detail/jacobian_line.hpp"
#include "detray/propagator/detail/jacobian_polar.hpp"
#include "detray/tracks/detail/transform_track_parameters.hpp"

namespace detray::detail {

/// @brief Generate Jacobians
template <typename algebra_t>
struct jacobian_engine {
  /// @name Type definitions for the struct
  /// @{
  using algebra_type = algebra_t;
  using transform3_type = dtransform3D<algebra_type>;
  using scalar_type = dscalar<algebra_type>;
  using point3_type = dpoint3D<algebra_type>;
  using vector3_type = dvector3D<algebra_type>;

  using bound_to_free_matrix_type = bound_to_free_matrix<algebra_type>;
  using free_to_bound_matrix_type = free_to_bound_matrix<algebra_type>;
  using free_to_path_matrix_type = free_to_path_matrix<algebra_type>;
  using path_to_free_matrix_type = path_to_free_matrix<algebra_type>;
  using bound_to_free_jacobian_submatrix_type =
      bound_to_free_jacobian_submatrix<algebra_type>;
  using free_to_bound_jacobian_submatrix_type =
      free_to_bound_jacobian_submatrix<algebra_type>;
  /// @}

  template <typename frame_t, typename mask_t>
    requires concepts::point<typename frame_t::loc_point>
  DETRAY_HOST_DEVICE static constexpr bound_to_free_matrix_type
  bound_to_free_jacobian(
      const transform3_type& trf3, const mask_t& mask,
      const bound_parameters_vector<algebra_type>& bound_vec) {
    // Global position and direction
    const free_track_parameters<algebra_type> free_params =
        bound_to_free_vector(trf3, mask, bound_vec);

    const vector3_type pos = free_params.pos();
    const vector3_type dir = free_params.dir();

    bound_to_free_matrix_type jac_to_global =
        matrix::zero<bound_to_free_matrix_type>();

    const bound_to_free_jacobian_submatrix_type dpos_dloc =
        bound_to_free_jacobian_submatrix_dpos_dloc<frame_t>(trf3, pos, dir);

    getter::set_block(jac_to_global, dpos_dloc, e_free_pos0, e_bound_loc0);

    const bound_to_free_jacobian_submatrix_type ddir_dangle =
        bound_to_free_jacobian_submatrix_ddir_dangle(bound_vec);

    getter::set_block(jac_to_global, ddir_dangle, e_free_dir0, e_bound_phi);

    const bound_to_free_jacobian_submatrix_type dpos_dangle =
        bound_to_free_jacobian_submatrix_dpos_dangle<frame_t>(trf3, pos, dir,
                                                              ddir_dangle);

    getter::set_block(jac_to_global, dpos_dangle, e_free_pos0, e_bound_phi);

    // Set d(bound time)/d(free time) and d(qop)/d(qop)
    getter::element(jac_to_global, e_free_time, e_bound_time) = 1.f;
    getter::element(jac_to_global, e_free_qoverp, e_bound_qoverp) = 1.f;

    return jac_to_global;
  }

  DETRAY_HOST_DEVICE static constexpr bound_to_free_jacobian_submatrix_type
  bound_to_free_jacobian_submatrix_ddir_dangle(
      const bound_parameters_vector<algebra_type>& bound_vec) {
    bound_to_free_jacobian_submatrix_type bound_angle_to_free_dir_derivative =
        matrix::zero<bound_to_free_jacobian_submatrix_type>();

    // Get trigonometric values
    const scalar_type theta{bound_vec.theta()};
    const scalar_type phi{bound_vec.phi()};
    const scalar_type cos_theta{math::cos(theta)};
    const scalar_type sin_theta{math::sin(theta)};
    const scalar_type cos_phi{math::cos(phi)};
    const scalar_type sin_phi{math::sin(phi)};

    // Assert the integrity of the local matrix wrt to the globally
    // defined track parameterization.
    static_assert(e_bound_theta == e_bound_phi + 1);
    static_assert(e_free_dir1 == e_free_dir0 + 1);
    static_assert(e_free_dir2 == e_free_dir0 + 2);

    constexpr unsigned int e_submatrix_bound_phi = 0u;
    constexpr unsigned int e_submatrix_bound_theta = 1u;
    constexpr unsigned int e_submatrix_free_dir0 = 0u;
    constexpr unsigned int e_submatrix_free_dir1 = 1u;
    constexpr unsigned int e_submatrix_free_dir2 = 2u;

    // Set d(n_x,n_y,n_z)/d(phi, theta)
    getter::element(bound_angle_to_free_dir_derivative, e_submatrix_free_dir0,
                    e_submatrix_bound_phi) = -sin_theta * sin_phi;
    getter::element(bound_angle_to_free_dir_derivative, e_submatrix_free_dir0,
                    e_submatrix_bound_theta) = cos_theta * cos_phi;
    getter::element(bound_angle_to_free_dir_derivative, e_submatrix_free_dir1,
                    e_submatrix_bound_phi) = sin_theta * cos_phi;
    getter::element(bound_angle_to_free_dir_derivative, e_submatrix_free_dir1,
                    e_submatrix_bound_theta) = cos_theta * sin_phi;
    getter::element(bound_angle_to_free_dir_derivative, e_submatrix_free_dir2,
                    e_submatrix_bound_theta) = -sin_theta;

    return bound_angle_to_free_dir_derivative;
  }

  template <typename frame_t>
    requires concepts::point<typename frame_t::loc_point>
  DETRAY_HOST_DEVICE static constexpr bound_to_free_jacobian_submatrix_type
  bound_to_free_jacobian_submatrix_dpos_dloc(const transform3_type& trf3,
                                             const vector3_type& pos,
                                             const vector3_type& dir) {
    using jacobian_t = jacobian<frame_t>;

    // Get d(x,y,z)/d(loc0, loc1)
    return jacobian_t::get_derivative_dpos_dloc(trf3, pos, dir);
  }

  template <typename frame_t>
    requires concepts::point<typename frame_t::loc_point>
  DETRAY_HOST_DEVICE static constexpr bound_to_free_jacobian_submatrix_type
  bound_to_free_jacobian_submatrix_dpos_dangle(
      const transform3_type& trf3, const vector3_type& pos,
      const vector3_type& dir,
      const bound_to_free_jacobian_submatrix_type& ddir_dangle) {
    using jacobian_t = jacobian<frame_t>;

    // Get d(x,y,z)/d(loc0, loc1)
    return jacobian_t::get_derivative_dpos_dangle(trf3, pos, dir, ddir_dangle);
  }

  DETRAY_HOST_DEVICE static constexpr free_to_bound_jacobian_submatrix_type
  free_to_bound_jacobian_submatrix_dangle_ddir(const vector3_type& dir) {
    free_to_bound_jacobian_submatrix_type dangle_ddir =
        matrix::zero<free_to_bound_jacobian_submatrix_type>();

    const scalar_type theta{vector::theta(dir)};
    const scalar_type phi{vector::phi(dir)};

    const scalar_type cos_theta{math::cos(theta)};
    const scalar_type sin_theta{math::sin(theta)};
    const scalar_type cos_phi{math::cos(phi)};
    const scalar_type sin_phi{math::sin(phi)};

    // Assert the integrity of the local matrix wrt to the globally
    // defined track parameterization.
    static_assert(e_bound_theta == e_bound_phi + 1);
    static_assert(e_free_dir1 == e_free_dir0 + 1);
    static_assert(e_free_dir2 == e_free_dir0 + 2);

    constexpr unsigned int e_submatrix_bound_phi = 0u;
    constexpr unsigned int e_submatrix_bound_theta = 1u;
    constexpr unsigned int e_submatrix_free_dir0 = 0u;
    constexpr unsigned int e_submatrix_free_dir1 = 1u;
    constexpr unsigned int e_submatrix_free_dir2 = 2u;

    // Set d(phi, theta)/d(n_x, n_y, n_z)
    // @note This codes have a serious bug when theta is equal to zero...
    getter::element(dangle_ddir, e_submatrix_bound_phi, e_submatrix_free_dir0) =
        -sin_phi / sin_theta;
    getter::element(dangle_ddir, e_submatrix_bound_phi, e_submatrix_free_dir1) =
        cos_phi / sin_theta;
    getter::element(dangle_ddir, e_submatrix_bound_theta,
                    e_submatrix_free_dir0) = cos_phi * cos_theta;
    getter::element(dangle_ddir, e_submatrix_bound_theta,
                    e_submatrix_free_dir1) = sin_phi * cos_theta;
    getter::element(dangle_ddir, e_submatrix_bound_theta,
                    e_submatrix_free_dir2) = -sin_theta;

    return dangle_ddir;
  }

  template <typename frame_t>
    requires concepts::point<typename frame_t::loc_point>
  DETRAY_HOST_DEVICE static constexpr free_to_bound_matrix_type
  free_to_bound_jacobian(
      const transform3_type& trf3,
      const free_track_parameters<algebra_type>& free_params) {
    // Global position and direction
    const vector3_type pos = free_params.pos();
    const vector3_type dir = free_params.dir();

    // Declare jacobian for bound to free coordinate transform
    free_to_bound_matrix_type jac_to_local =
        matrix::zero<free_to_bound_matrix_type>();

    const free_to_bound_jacobian_submatrix_type dloc_dpos =
        free_to_bound_jacobian_submatrix_dloc_dpos<frame_t>(trf3, pos, dir);

    getter::set_block(jac_to_local, dloc_dpos, e_bound_loc0, e_free_pos0);

    const free_to_bound_jacobian_submatrix_type dangle_ddir =
        free_to_bound_jacobian_submatrix_dangle_ddir(dir);

    getter::set_block(jac_to_local, dangle_ddir, e_bound_phi, e_free_dir0);

    // Set d(free time)/d(bound time)
    getter::element(jac_to_local, e_bound_time, e_free_time) = 1.f;
    // Set d(Free Qop)/d(Bound Qop)
    getter::element(jac_to_local, e_bound_qoverp, e_free_qoverp) = 1.f;

    return jac_to_local;
  }

  template <typename frame_t>
    requires concepts::point<typename frame_t::loc_point>
  DETRAY_HOST_DEVICE static constexpr free_to_bound_jacobian_submatrix_type
  free_to_bound_jacobian_submatrix_dloc_dpos(const transform3_type& trf3,
                                             const vector3_type& pos,
                                             const vector3_type& dir) {
    using jacobian_t = jacobian<frame_t>;

    // Get d(loc0, loc1)/d(x, y, z)
    return jacobian_t::get_derivative_dloc_dpos(trf3, pos, dir);
  }

  template <typename frame_t>
    requires concepts::point<typename frame_t::loc_point>
  DETRAY_HOST_DEVICE static constexpr free_to_path_matrix_type
  free_to_path_derivative(const vector3_type& pos, const vector3_type& dir,
                          const vector3_type& dtds,
                          const transform3_type& trf3) {
    using jacobian_t = jacobian<frame_t>;

    return jacobian_t::path_derivative(trf3, pos, dir, dtds);
  }

  DETRAY_HOST_DEVICE static constexpr path_to_free_matrix_type
  path_to_free_derivative(const vector3_type& dir, const vector3_type& dtds,
                          const scalar_type dqopds) {
    path_to_free_matrix_type derivative =
        matrix::zero<path_to_free_matrix_type>();
    getter::element(derivative, e_free_pos0, 0u) = dir[0];
    getter::element(derivative, e_free_pos1, 0u) = dir[1];
    getter::element(derivative, e_free_pos2, 0u) = dir[2];
    getter::element(derivative, e_free_dir0, 0u) = dtds[0];
    getter::element(derivative, e_free_dir1, 0u) = dtds[1];
    getter::element(derivative, e_free_dir2, 0u) = dtds[2];
    getter::element(derivative, e_free_qoverp, 0u) = dqopds;

    return derivative;
  }

  template <typename frame_t>
    requires concepts::point<typename frame_t::loc_point>
  DETRAY_HOST_DEVICE static constexpr free_matrix<algebra_type> path_correction(
      const vector3_type& pos, const vector3_type& dir,
      const vector3_type& dtds, const scalar_type dqopds,
      const transform3_type& trf3) {
    return path_to_free_derivative(dir, dtds, dqopds) *
           free_to_path_derivative<frame_t>(pos, dir, dtds, trf3);
  }
};

}  // namespace detray::detail
