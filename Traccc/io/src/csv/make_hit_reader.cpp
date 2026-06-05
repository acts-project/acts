/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/io/csv/make_hit_reader.hpp"

namespace traccc::io::csv {

dfe::NamedTupleCsvReader<hit> make_hit_reader(std::string_view filename) {

    return {filename.data(),
            {"particle_id", "geometry_id", "tx", "ty", "tz", "tt", "tpx", "tpy",
             "tpz", "te", "deltapx", "deltapy", "deltapz", "deltae", "index"}};
}

}  // namespace traccc::io::csv
