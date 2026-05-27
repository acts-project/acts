/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/utils/detray_conversion.hpp"

namespace traccc::edm {

template <detray::concepts::algebra algebra_t, typename spacepoint_backend_t>
TRACCC_HOST_DEVICE detray::dpoint3D<algebra_t> get_spacepoint_global(
    const edm::spacepoint<spacepoint_backend_t>& sp) {

    return utils::to_dpoint3D<algebra_t>(sp.global());
}

}  // namespace traccc::edm
