/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/io/csv/dfe.hpp"
#include "traccc/io/csv/measurement_hit_id.hpp"

// System include(s).
#include <string_view>

namespace traccc::io::csv {

/// Set up a reader for measurement<->hit association information
///
/// @param filename The name of the file to read
/// @return An object that can read the specified CSV file
///
dfe::NamedTupleCsvReader<measurement_hit_id> make_measurement_hit_id_reader(
    std::string_view filename);

}  // namespace traccc::io::csv
