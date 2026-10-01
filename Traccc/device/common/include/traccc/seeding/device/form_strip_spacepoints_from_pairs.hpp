/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/device/global_index.hpp"
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/spacepoint_collection.hpp"
#include "traccc/seeding/detail/strip_pair.hpp"
#include "traccc/seeding/strip_spacepoint_formation_data.hpp"

namespace traccc::device {

/// Form a spacepoint from a compatible pair using static strip geometry.
/// The pair carries measurement indices and explicit endpoint tolerances.
TRACCC_HOST_DEVICE inline void form_strip_spacepoints_from_pairs(
    global_index_t globalIndex,
    const edm::measurement_collection::const_view& measurements_view,
    const strip_pair_collection_types::const_view& pairs_view,
    const strip_measurement_surface_info_collection_types::const_view&
        surface_infos_view,
    const point3& beam_spot, edm::spacepoint_collection::view spacepoints_view);

}  // namespace traccc::device

// Include the implementation.
#include "traccc/seeding/device/impl/form_strip_spacepoints_from_pairs.ipp"
