/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "read_surfaces.hpp"

#include "traccc/io/csv/dfe.hpp"
#include "traccc/io/csv/make_surface_reader.hpp"

namespace traccc::io::csv {

std::map<geometry_id, transform3> read_surfaces(std::string_view filename) {

    // Construct the surface reader object.
    auto reader = make_surface_reader(filename);

    // Fill an std::map with entries from the file.
    std::map<geometry_id, transform3> result;
    surface surf;
    while (reader.read(surf)) {
        result.insert(
            {surf.geometry_id,
             transform3{vector3{surf.cx, surf.cy, surf.cz},
                        vector3{surf.rot_xw, surf.rot_yw, surf.rot_zw},
                        vector3{surf.rot_xu, surf.rot_yu, surf.rot_zu}}});
    }
    return result;
}

}  // namespace traccc::io::csv
