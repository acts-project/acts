/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/io/data_format.hpp"

// Project include(s).
#include "traccc/edm/silicon_cell_collection.hpp"
#include "traccc/geometry/detector_conditions_description.hpp"
#include "traccc/geometry/detector_design_description.hpp"
#include "traccc/utils/logging.hpp"

// System include(s).
#include <cstddef>
#include <string_view>

namespace traccc::io {

/// Read cell data into memory
///
/// The file to read is selected according the naming conventions used in
/// our data.
///
/// @param[out] cells       The cell collection to fill
/// @param[in]  event       The event ID to read in the cells for
/// @param[in]  directory   The directory holding the cell data files
/// @param[in]  dd          The detector description to point the cells at
/// @param[in]  format      The format of the cell data files (to read)
/// @param[in]  deduplicate Whether to deduplicate the cells
/// @param[in]  use_acts_geometry_id Whether to treat the geometry ID as an
///                                  "Acts geometry ID", or a
///                                  "Detray geometry ID"
///
void read_cells(edm::silicon_cell_collection::host& cells, std::size_t event,
                std::string_view directory,
                std::unique_ptr<const Logger> logger = getDummyLogger().clone(),
                const detector_conditions_description::host* det_cond = nullptr,
                data_format format = data_format::csv, bool deduplicate = true,
                bool use_acts_geometry_id = true);

/// Read cell data into memory
///
/// The file name is selected explicitly by the user.
///
/// @param[out] cells       The cell collection to fill
/// @param[in]  filename    The name of the file to read
/// @param[in]  det_cond    The detector conditions description
/// @param[in]  format      The format of the cell data files (to read)
/// @param[in]  deduplicate Whether to deduplicate the cells
/// @param[in]  use_acts_geometry_id Whether to treat the geometry ID as an
///                                  "Acts geometry ID", or a
///                                  "Detray geometry ID"
///
void read_cells(edm::silicon_cell_collection::host& cells,
                std::string_view filename,
                std::unique_ptr<const Logger> logger = getDummyLogger().clone(),
                const detector_conditions_description::host* det_cond = nullptr,
                data_format format = data_format::csv, bool deduplicate = true,
                bool use_acts_geometry_id = true);

}  // namespace traccc::io
