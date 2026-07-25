/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/io/digitization_config.hpp"

// System include(s).
#include <string_view>

namespace traccc::io::json {

/// Write a digitization configuration to a file
///
/// @param filename The name of the file to write the data to
/// @param config The digitization configuration to write
///
void write_digitization_config(std::string_view filename,
                               const digitization_config& config);

}  // namespace traccc::io::json
