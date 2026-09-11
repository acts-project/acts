/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/io/csv/cell.hpp"
#include "traccc/io/csv/dfe.hpp"

// System include(s).
#include <string_view>

namespace traccc::io::csv {

/// Set up an object for reading a CSV file containing cell information
///
/// @param filename The name of the file to read
/// @return An object that can read the specified CSV file
///
dfe::NamedTupleCsvReader<cell> make_cell_reader(std::string_view filename);

}  // namespace traccc::io::csv
