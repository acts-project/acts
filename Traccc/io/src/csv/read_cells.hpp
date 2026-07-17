/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/edm/silicon_cell_collection.hpp"
#include "traccc/geometry/detector_conditions_description.hpp"
#include "traccc/geometry/detector_design_description.hpp"
#include "traccc/utils/logging.hpp"

// System include(s).
#include <string_view>

namespace traccc::io::csv {

/// Read cell information from a specific CSV file
///
/// @param[out] cells       The cell collection to fill
/// @param[in]  filename    The name of the file to read
/// @param[in]  dd          The detector description to point the cells at
/// @param[in]  deduplicate Whether to deduplicate the cells
/// @param[in]  use_acts_geometry_id Whether to treat the geometry ID as an
///                                  "Acts geometry ID", or a
///                                  "Detray geometry ID"
///
void read_cells(edm::silicon_cell_collection::host& cells,
                std::string_view filename,
                std::unique_ptr<const Logger> logger = getDummyLogger().clone(),
                const detector_conditions_description::host* det_cond = nullptr,
                bool deduplicate = true, bool use_acts_geometry_id = true);

}  // namespace traccc::io::csv
