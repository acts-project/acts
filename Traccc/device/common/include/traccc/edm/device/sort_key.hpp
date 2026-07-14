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

template <detray::concepts::algebra algebra_t>
TRACCC_HOST_DEVICE inline sort_key get_sort_key(
    const bound_track_parameters<algebra_t>& params) {
    // key = |theta - pi/2|
    return math::fabs(params.theta() - constant<traccc::scalar>::pi_2);
}

}  // namespace traccc::device
