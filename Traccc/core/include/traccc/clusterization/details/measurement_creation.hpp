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

/// Function used for calculating the properties of the cluster during
/// measurement creation
///
/// @param[in] cluster      The silicon cluster to calculate the properties of
/// @param[in] cells        All silicon cells in the event
/// @param[in] det_descr    The detector segmentation description
/// @param[in] det_cond     The detector conditionas data
/// @param[out] mean        The mean position of the cluster/measurement
/// @param[out] var         The variation on the mean position of the
///                         cluster/measurement
/// @param[out] totalWeight The total weight of the cluster/measurement
///
template <typename T, typename TDesign>
TRACCC_HOST_DEVICE inline void calc_cluster_properties(
    const edm::silicon_cluster<T>& cluster,
    const edm::silicon_cell_collection::const_device& cells,
    const traccc::detector_design_description_interface<TDesign>& module_dd,
    point2& mean, point2& var, scalar& totalWeight);

/// Function used for calculating the properties of the cluster during
/// measurement creation
///
/// @param[out] measurement Measurement object to be filled
/// @param[in] cluster   The silicon cluster to turn into a measurement
/// @param[in] index     The index of the cluster/measurement in the collection
/// @param[in] cells     All silicon cells in the event
/// @param[in] det_descr Detector segmentation description
/// @param[in] det_cond  Detector conditions description
///
template <typename T1, typename T2>
TRACCC_HOST_DEVICE inline void fill_measurement(
    edm::measurement<T1>& measurement, const edm::silicon_cluster<T2>& cluster,
    unsigned int index, const edm::silicon_cell_collection::const_device& cells,
    const detector_design_description::const_device& det_descr,
    const detector_conditions_description::const_device& det_cond);

}  // namespace traccc::details

// Include the implementation.
#include "traccc/clusterization/impl/measurement_creation.ipp"
