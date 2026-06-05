/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Library include(s).
#include "traccc/io/data_format.hpp"

// Project include(s).
#include "traccc/geometry/detector_conditions_description.hpp"
#include "traccc/geometry/detector_design_description.hpp"

// System include(s).
#include <cstdint>
#include <map>
#include <string_view>

namespace traccc::io {

/// Populate a @c traccc::detector_design_description object from text files.
///
/// @param dd The detector design description object to set up.
/// @param geometry_file The path to the geometry description file.
/// @param digitization_file The path to the digitization configuration file.
/// @param conditions_file The path to the conditions configuration file.
/// @param geometry_format The format of the geometry description file.
/// @param digitization_format The format of the digitization configuration
///                            file.
/// @param conditions_format The format of the conditions configuration
///                           file.
///
void read_detector_description(
    detector_design_description::host& det_desc,
    detector_conditions_description::host& det_cond,
    std::string_view geometry_file, std::string_view digitization_file,
    std::string_view conditions_file,
    data_format geometry_format = data_format::json,
    data_format digitization_format = data_format::json,
    data_format conditions_format = data_format::json);

}  // namespace traccc::io
