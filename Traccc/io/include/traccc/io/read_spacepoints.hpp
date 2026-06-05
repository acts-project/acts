/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/io/data_format.hpp"

// Project include(s).
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/spacepoint_collection.hpp"
#include "traccc/geometry/detector.hpp"
#include "traccc/geometry/detector_conditions_description.hpp"
#include "traccc/geometry/detector_design_description.hpp"
#include "traccc/geometry/host_detector.hpp"

// System include(s).
#include <cstddef>
#include <string_view>

namespace traccc::io {

/// Read spacepoint data into memory
///
/// The file to read is selected according the naming conventions used in
/// our data.
///
/// @param[out] spacepoints The spacepoint collection to fill
/// @param[out] measurements The measurement collection to fill
/// @param[in]  event     The event ID to read in the spacepoints for
/// @param[in]  directory The directory holding the spacepoint data files
/// @param[in]  detector  detray detector
/// @param[in]  design_description The detector segmnentation information
/// @param[in]  conditions_description The detector conditions description
/// @param[in]  format    The format of the data files (to read)
///
void read_spacepoints(edm::spacepoint_collection::host& spacepoints,
                      edm::measurement_collection::host& measurements,
                      std::size_t event, std::string_view directory,
                      const traccc::host_detector* detector = nullptr,
                      const traccc::detector_design_description::host*
                          detector_design_description = nullptr,
                      const traccc::detector_conditions_description::host*
                          detector_conditions_description = nullptr,
                      data_format format = data_format::csv);

/// Read spacepoint data into memory
///
/// The file name is selected explicitly by the user.
///
/// @param[out] spacepoints The spacepoint collection to fill
/// @param[out] measurements The measurement collection to fill
/// @param[in]  hit_filename  The file to read the hit/spacepoint data from
/// @param[in]  meas_filename The file to read the measurement data from
/// @param[in]  meas_hit_map_filename The file to read the mapping from
///                                   measurements to hits from
/// @param[in]  detector  detray detector
/// @param[in]  format The format of the data files (to read)
///
void read_spacepoints(edm::spacepoint_collection::host& spacepoints,
                      edm::measurement_collection::host& measurements,
                      std::string_view hit_filename,
                      std::string_view meas_filename,
                      std::string_view meas_hit_map_filename,
                      const traccc::host_detector* detector = nullptr,
                      const traccc::detector_design_description::host*
                          detector_design_description = nullptr,
                      const traccc::detector_conditions_description::host*
                          detector_conditions_description = nullptr,
                      data_format format = data_format::csv);

}  // namespace traccc::io
