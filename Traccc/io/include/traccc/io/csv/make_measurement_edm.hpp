/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/definitions/primitives.hpp"
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/geometry/detector_conditions_description.hpp"
#include "traccc/geometry/detector_design_description.hpp"

// Local include(s).
#include "traccc/io/csv/measurement.hpp"

#pragma once

namespace traccc::io::csv {

/// Make measurement EDM from csv measurement
///
/// @param[in] csv_meas input csv measurement
/// @param[out] meas output measurement proxy to fill
/// @param[in] acts_to_detray_id Map for acts-to-detray geometry ID converision
///
void make_measurement_edm(
    const traccc::io::csv::measurement& csv_meas,
    edm::measurement_collection::host::proxy_type& meas,
    const std::map<geometry_id, geometry_id>* acts_to_detray_id,
    const traccc::detector_design_description::host* det_desc = nullptr,
    const std::map<geometry_id, std::size_t>*
        geometry_id_to_detector_description_index = nullptr);

}  // namespace traccc::io::csv
