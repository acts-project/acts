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
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/geometry/coordinates/cartesian2D.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes/rectangle2D.hpp"
#include "detray/propagator/detail/jacobian_engine.hpp"
#include "detray/tracks/bound_track_parameters.hpp"
#include "detray/utils/unit_vectors.hpp"

namespace detray {

/// Frame and Jacobian data associated with a free track's curvilinear frame.
template <concepts::algebra algebra_t>
struct curvilinear_frame {
  /// Three-dimensional transform type for the curvilinear frame.
  using transform3_type = dtransform3D<algebra_t>;
  /// Three-dimensional vector type used by the transform.
  using vector3 = typename transform3_type::vector3;
  /// Helper for constructing the transverse unit vectors.
  using unit_vectors_type = unit_vectors<vector3>;
  /// Matrix type mapping bound parameter variations to free parameters.
  using bound_to_free_matrix_type = bound_to_free_matrix<algebra_t>;
  /// Vector holding the bound parameters at the frame's reference point.
  using bound_vector_type = bound_parameters_vector<algebra_t>;
  /// Engine that computes the bound-to-free Jacobian.
  using jacobian_engine_type = detail::jacobian_engine<algebra_t>;
  /// Type of the free track parameters used to initialize the frame.
  using free_track_parameters_type = free_track_parameters<algebra_t>;

  /// Build the frame and bound parameter vector from a valid free track.
  /// @param free_params Free track parameters defining the reference point.
  DETRAY_HOST_DEVICE
  explicit curvilinear_frame(const free_track_parameters_type& free_params) {
    assert(!free_params.is_invalid());

    const vector3 t = free_params.pos();
    const vector3 z = free_params.dir();
    const vector3 x = unit_vectors_type().make_curvilinear_unit_u(z);

    m_trf = transform3_type(t, z, x);
    m_bound_vec = detail::free_to_bound_vector<cartesian2D<algebra_t>>(
        m_trf, free_params);
  }

  /// Compute the Jacobian from bound curvilinear parameters to free parameters.
  /// @returns The bound-to-free Jacobian at the reference point.
  DETRAY_HOST_DEVICE
  bound_to_free_matrix_type bound_to_free_jacobian() const {
    return jacobian_engine_type()
        .template bound_to_free_jacobian<cartesian2D<algebra_t>>(
            m_trf, mask<rectangle2D, algebra_t>{}, m_bound_vec);
  }

  /// Transform placing the curvilinear frame at the track reference point.
  transform3_type m_trf{};
  /// Bound curvilinear parameters of the reference track.
  bound_vector_type m_bound_vec{};

};  // curvilinear frame

}  // namespace detray
