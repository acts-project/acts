/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/io/csv/make_surface_reader.hpp"

namespace traccc::io::csv {

dfe::NamedTupleCsvReader<surface> make_surface_reader(
    std::string_view filename) {

    return {filename.data(),
            {"geometry_id", "cx", "cy", "cz", "rot_xu", "rot_xv", "rot_xw",
             "rot_zu", "rot_zv", "rot_zw"}};
}

}  // namespace traccc::io::csv
