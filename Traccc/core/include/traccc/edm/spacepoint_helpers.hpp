/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/spacepoint_collection.hpp"

// Detray include(s).
#include <detray/definitions/algebra.hpp>

namespace traccc::edm {

/// Get the global position of a spacepoint as a 3D point
///
/// @tparam algebra_t The algebra type used to describe the tracks
///
/// @param sp The spacepoint to extract the global position from
/// @param pos The 3D point to fill with the global position of the spacepoint
///
template <detray::concepts::algebra algebra_t, typename spacepoint_backend_t>
TRACCC_HOST_DEVICE detray::dpoint3D<algebra_t> get_spacepoint_global(
    const edm::spacepoint<spacepoint_backend_t>& sp);

}  // namespace traccc::edm

// Implementation include(s).
#include "traccc/edm/impl/spacepoint_helpers.ipp"
