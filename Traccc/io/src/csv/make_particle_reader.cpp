/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/io/csv/make_particle_reader.hpp"

namespace traccc::io::csv {

dfe::NamedTupleCsvReader<particle> make_particle_reader(
    std::string_view filename) {

    return {filename.data(),
            {"particle_id", "particle_type", "process", "vx", "vy", "vz", "vt",
             "px", "py", "pz", "m", "q"}};
}

}  // namespace traccc::io::csv
