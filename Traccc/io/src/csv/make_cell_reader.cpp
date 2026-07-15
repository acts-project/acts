/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/io/csv/make_cell_reader.hpp"

namespace traccc::io::csv {

dfe::NamedTupleCsvReader<cell> make_cell_reader(std::string_view filename) {

    return {filename.data(),
            {"geometry_id", "measurement_id", "cannel0", "channel1",
             "timestamp", "value"}};
}

}  // namespace traccc::io::csv
