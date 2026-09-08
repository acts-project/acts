/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/primitives.hpp"
#include "traccc/edm/track_parameters.hpp"

namespace traccc::device {

using sort_key = traccc::scalar;

/*
 * Spacing between two consecutive volumes in the sort key.
 *
 * Has to be larger than the largest value the intra-volume part of the key
 * can take, i.e. max(|theta - pi/2|) = pi/2.
 */
inline constexpr sort_key sort_key_volume_stride{2.f};

/*
 * Offset that is added to the key of a dead track, so that dead tracks sort
 * behind every live one.
 *
 * Has to be larger than the largest key a live track can have, i.e.
 * sort_key_volume_stride * (number of volumes) + pi/2. Detray encodes the
 * volume in a 12 bit field, so there are at most 4096 of them (see the mask
 * layout in detray's geometry/identifier.hpp).
 */
inline constexpr sort_key dead_track_sort_key_offset{10000.f};

static_assert(
    sort_key_volume_stride * 4096.f + 2.f < dead_track_sort_key_offset,
    "Dead tracks are no longer guaranteed to sort behind the live ones");

template <detray::concepts::algebra algebra_t>
TRACCC_HOST_DEVICE inline sort_key get_sort_key(
    const bound_track_parameters<algebra_t>& params) {
  // key = volume * stride + |theta - pi/2|
  return static_cast<sort_key>(params.surface_link().volume()) *
             sort_key_volume_stride +
         math::fabs(params.theta() - constant<traccc::scalar>::pi_2);
}

}  // namespace traccc::device
