/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <vecmem/containers/data/vector_view.hpp>

#include "traccc/definitions/qualifiers.hpp"
#include "traccc/device/global_index.hpp"
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/geometry/detector.hpp"
#include "traccc/seeding/detail/strip_pair.hpp"

namespace traccc::device {

template <detray::concepts::detector detector_t>
TRACCC_HOST_DEVICE inline void count_strip_pairs(
    global_index_t globalIndex,
    const detray::detector_view_t<detector_t> det_view,
    const edm::measurement_collection::const_view& measurements_view,
    const strip_measurement_surface_info_collection_types::const_view&
        surface_infos_view,
    const strip_pairing_rule_collection_types::const_view& rules_view,
    const point3& beam_spot, unsigned int& n_opposite_pairs,
    unsigned int& n_overlap_pairs);

}  // namespace traccc::device

#include "traccc/seeding/device/impl/count_strip_pairs.ipp"
