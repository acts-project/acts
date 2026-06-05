/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/io/csv/make_measurement_hit_id_reader.hpp"

namespace traccc::io::csv {

dfe::NamedTupleCsvReader<measurement_hit_id> make_measurement_hit_id_reader(
    std::string_view filename) {

    return {filename.data(), {"measurement_id", "hit_id"}};
}

}  // namespace traccc::io::csv
