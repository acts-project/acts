/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/primitives.hpp"
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/silicon_cell_collection.hpp"
#include "traccc/edm/silicon_cluster_collection.hpp"
#include "traccc/geometry/detector_conditions_description.hpp"
#include "traccc/geometry/detector_design_description.hpp"

namespace traccc::details {

/// Get the local position of a cell on a module
///
/// @param cell      The cell to get the position of
/// @param det_descr The (silicon) detector description
/// @return The local position of the cell (upper bound) and optionality the
/// lower bound
///
template <typename TCell, typename TDesign>
TRACCC_HOST_DEVICE inline vector2 position_from_cell(
    const edm::silicon_cell<TCell>& cell,
    const traccc::detector_design_description_interface<TDesign>& module_dd);

}  // namespace traccc::details

// Include the implementation.
#include "traccc/clusterization/impl/measurement_creation.ipp"
