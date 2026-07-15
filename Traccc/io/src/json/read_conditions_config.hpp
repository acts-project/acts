/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/io/conditions_config.hpp"

// System include(s).
#include <string_view>

namespace traccc::io::json {

/// Read the detector conditions configuration from a JSON input file
///
/// @param filename The name of the file to read the data from
/// @return An object describing the conditions configuration of the
///         detector
///

conditions_config read_conditions_config(std::string_view filename);

}  // namespace traccc::io::json
