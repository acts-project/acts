// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/track_parametrization.hpp"
#include "detray/tracks/tracks.hpp"

namespace detray::detail {

/// Transform a free track parameter vector to a bound track parameter vector
///
/// @param trf3 transform of the surface the bound vector should be defined on
/// @param free_param the free track parameters to be transformed
///
/// @returns the bound track parameter vector
template <typename local_frame_t>
DETRAY_HOST_DEVICE inline auto free_to_bound_vector(
    const dtransform3D<typename local_frame_t::algebra_type>& trf3,
    const free_parameters_vector<typename local_frame_t::algebra_type>&
        free_vec) {
  using algebra_t = typename local_frame_t::algebra_type;

  const auto pos = free_vec.pos();
  const auto dir = free_vec.dir();

  const auto bound_local = local_frame_t::global_to_local(trf3, pos, dir);

  return bound_parameters_vector<algebra_t>{bound_local, vector::phi(dir),
                                            vector::theta(dir), free_vec.qop(),
                                            free_vec.time()};
}

/// Transform a bound track parameter vector to a free track parameter vector
///
/// @param trf3 transform of the surface the bound parameters are defined on
/// @param mask the mask of the surface the bound parameters are defined on
/// @param bound_vec the bound track vector to be transformed
///
/// @returns the free track parameter vector
template <typename mask_t>
DETRAY_HOST_DEVICE inline auto bound_to_free_vector(
    const dtransform3D<typename mask_t::algebra_type>& trf3, const mask_t& mask,
    const bound_parameters_vector<typename mask_t::algebra_type>& bound_vec) {
  using algebra_t = typename mask_t::algebra_type;
  using local_frame_t = typename mask_t::local_frame;

  const auto bound_local = bound_vec.bound_local();
  const auto dir = bound_vec.dir();

  const auto pos = local_frame_t::local_to_global(trf3, mask, bound_local, dir);

  // The free vector constructor expects momentum and charge, so set the
  // values explicitly instead
  free_parameters_vector<algebra_t> free_vec{};

  free_vec.set_pos(pos);
  free_vec.set_time(bound_vec.time());
  free_vec.set_dir(dir);
  free_vec.set_qop(bound_vec.qop());

  return free_vec;
}

}  // namespace detray::detail
