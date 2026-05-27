/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/io/data_format.hpp"

// Project include(s).
#include "traccc/io/digitization_config.hpp"

// System include(s).
#include <string_view>

namespace traccc::io {

/// Read the detector digitization configuration from an input file
///
/// @param filename The name of the file to read the data from
/// @param format The format of the input file
/// @return An object describing the digitization configuration of the
///         detector
///
digitization_config read_digitization_config(
    std::string_view filename, data_format format = data_format::json);

}  // namespace traccc::io
