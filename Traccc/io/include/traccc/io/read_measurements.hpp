/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/geometry/detector_conditions_description.hpp"
#include "traccc/geometry/detector_design_description.hpp"
#include "traccc/io/data_format.hpp"

// Project include(s).
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/geometry/detector.hpp"
#include "traccc/geometry/host_detector.hpp"

// System include(s).
#include <cstddef>
#include <string_view>
#include <vector>

namespace traccc::io {

/// Read measurement data into memory
///
/// The file to read is selected according the naming conventions used in
/// our data.
///
/// @param[out] measurements The measurement collection to fill
/// @param[in]  event     The event ID to read in the measurements for
/// @param[in]  directory The directory holding the measurement data files
/// @param[in]  detector  detray detector
/// @param[in]  format    The format of the measurement data files (to read)
///
std::vector<measurement_id_type> read_measurements(
    edm::measurement_collection::host& measurements, std::size_t event,
    std::string_view directory, const traccc::host_detector* detector = nullptr,
    const traccc::detector_design_description::host* det_desc = nullptr,
    const traccc::detector_conditions_description::host* det_cond = nullptr,
    const bool sort_measurements = true, data_format format = data_format::csv);

/// Read measurement data into memory
///
/// The file name is selected explicitly by the user.
///
/// @param[out] measurements The measurement collection to fill
/// @param[in]  filename The file to read the measurement data from
/// @param[in]  detector  detray detector
/// @param[in]  format   The format of the measurement data files (to read)
///
std::vector<measurement_id_type> read_measurements(
    edm::measurement_collection::host& measurements, std::string_view filename,
    const traccc::host_detector* detector = nullptr,
    const traccc::detector_design_description::host* det_desc = nullptr,
    const traccc::detector_conditions_description::host* det_cond = nullptr,
    const bool sort_measurements = true, data_format format = data_format::csv);

}  // namespace traccc::io
