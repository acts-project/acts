/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/device/global_index.hpp"

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/spacepoint_collection.hpp"
#include "traccc/geometry/detector_conditions_description.hpp"
#include "traccc/geometry/detector_design_description.hpp"

namespace traccc::device {

/// Function for creating 3D spacepoints out of 2D measurements
///
/// @param[in] globalIndex          The index for the current thread
/// @param[in] det_view             A view type object of detector
/// @param[in] measurements_view    Collection of measurements
/// @param[out] spacepoints_view    Collection of spacepoints
///
template <typename detector_t>
TRACCC_HOST_DEVICE inline void form_spacepoints(
    global_index_t globalIndex, typename detector_t::view det_view,
    const edm::measurement_collection::const_view& measurements_view,
    edm::spacepoint_collection::view spacepoints_view);

}  // namespace traccc::device

// Include the implementation.
#include "traccc/seeding/device/impl/form_spacepoints.ipp"
