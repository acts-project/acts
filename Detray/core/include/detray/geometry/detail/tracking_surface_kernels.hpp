// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/track_parametrization.hpp"
#include "detray/geometry/detail/surface_kernels.hpp"
#include "detray/propagator/detail/jacobian_engine.hpp"
#include "detray/tracks/detail/transform_track_parameters.hpp"
#include "detray/tracks/tracks.hpp"

// System include(s)
#include <ostream>

namespace detray::detail {

/// Functors to be used in the @c tracking_surface class
template <typename algebra_t>
struct tracking_surface_kernels : public surface_kernels<algebra_t> {
  using vector3_type = dvector3D<algebra_t>;
  using transform3_type = dtransform3D<algebra_t>;
  using bound_param_vector_type = bound_parameters_vector<algebra_t>;
  using free_param_vector_type = free_parameters_vector<algebra_t>;
  using free_matrix_type = free_matrix<algebra_t>;
  using free_to_path_matrix_type = free_to_path_matrix<algebra_t>;

  /// A functor to get from a free to a bound vector
  struct free_to_bound_vector {
    // Visitor to the detector mask store that is called on the mask
    // collection that contains the mask (shape) type of the surface
    template <typename mask_group_t, typename index_t>
    DETRAY_HOST_DEVICE inline bound_param_vector_type operator()(
        const mask_group_t& /*mask_group*/, const index_t& /*index*/,
        const transform3_type& trf3,
        const free_param_vector_type& free_vec) const {
      using frame_t = typename mask_group_t::value_type::local_frame;

      return detail::free_to_bound_vector<frame_t>(trf3, free_vec);
    }
  };

  /// A functor to get from a bound to a free vector
  struct bound_to_free_vector {
    template <typename mask_group_t, concepts::index index_t>
    DETRAY_HOST_DEVICE inline free_param_vector_type operator()(
        const mask_group_t& mask_group, const index_t& index,
        const transform3_type& trf3,
        const bound_param_vector_type& bound_vec) const {
      return detail::bound_to_free_vector(trf3, mask_group[index], bound_vec);
    }

    template <typename mask_group_t, concepts::interval idx_range_t>
    DETRAY_HOST_DEVICE inline free_param_vector_type operator()(
        const mask_group_t& mask_group, const idx_range_t& idx_range,
        const transform3_type& trf3,
        const bound_param_vector_type& bound_vec) const {
      return detail::bound_to_free_vector(trf3, mask_group[idx_range.lower()],
                                          bound_vec);
    }
  };

  /// A functor to get the free-to-bound Jacobian
  struct free_to_bound_jacobian {
    template <typename mask_group_t, typename index_t>
    DETRAY_HOST_DEVICE inline auto operator()(
        const mask_group_t& /*mask_group*/, const index_t& /*index*/,
        const transform3_type& trf3,
        const free_param_vector_type& free_vec) const {
      using frame_t = typename mask_group_t::value_type::local_frame;

      return detail::jacobian_engine<
          algebra_t>::template free_to_bound_jacobian<frame_t>(trf3, free_vec);
    }
  };

  /// A functor to get the bound-to-free Jacobian
  struct bound_to_free_jacobian {
    template <typename mask_group_t, concepts::index index_t>
    DETRAY_HOST_DEVICE inline auto operator()(
        const mask_group_t& mask_group, const index_t& index,
        const transform3_type& trf3,
        const bound_param_vector_type& bound_vec) const {
      using frame_t = typename mask_group_t::value_type::local_frame;

      return detail::jacobian_engine<algebra_t>::
          template bound_to_free_jacobian<frame_t>(trf3, mask_group[index],
                                                   bound_vec);
    }

    template <typename mask_group_t, concepts::interval idx_range_t>
    DETRAY_HOST_DEVICE inline auto operator()(
        const mask_group_t& mask_group, const idx_range_t& idx_range,
        const transform3_type& trf3,
        const bound_param_vector_type& bound_vec) const {
      using frame_t = typename mask_group_t::value_type::local_frame;

      return detail::jacobian_engine<algebra_t>::
          template bound_to_free_jacobian<frame_t>(
              trf3, mask_group[idx_range.lower()], bound_vec);
    }
  };

  /// A functor to get the path correction
  struct path_correction {
    template <typename mask_group_t, typename index_t,
              concepts::scalar scalar_t>
    DETRAY_HOST_DEVICE inline free_matrix_type operator()(
        const mask_group_t& /*mask_group*/, const index_t& /*index*/,
        const transform3_type& trf3, const vector3_type& pos,
        const vector3_type& dir, const vector3_type& dtds,
        const scalar_t dqopds) const {
      using frame_t = typename mask_group_t::value_type::local_frame;

      return detail::jacobian_engine<algebra_t>::template path_correction<
          frame_t>(pos, dir, dtds, dqopds, trf3);
    }
  };

  /// A function object to get the free to path derivative
  struct free_to_path_derivative {
    template <typename mask_group_t, typename index_t>
    DETRAY_HOST_DEVICE inline free_to_path_matrix_type operator()(
        const mask_group_t& /*mask_group*/, const index_t& /*index*/,
        const transform3_type& trf3, const vector3_type& pos,
        const vector3_type& dir, const vector3_type& dtds) const {
      using frame_t = typename mask_group_t::value_type::local_frame;

      return detail::jacobian_engine<
          algebra_t>::template free_to_path_derivative<frame_t>(pos, dir, dtds,
                                                                trf3);
    }
  };
};
}  // namespace detray::detail
