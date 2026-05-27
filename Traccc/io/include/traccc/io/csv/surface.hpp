/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "traccc/io/csv/dfe.hpp"

// System include(s).
#include <cstdint>

namespace traccc::io::csv {

/// Type to read information into about the detector modules/surfaces
struct surface {

    uint64_t geometry_id = 0;
    float cx = 0.f, cy = 0.f, cz = 0.f;
    float rot_xu = 0.f, rot_xv = 0.f, rot_xw = 0.f;
    float rot_yu = 0.f, rot_yv = 0.f, rot_yw = 0.f;
    float rot_zu = 0.f, rot_zv = 0.f, rot_zw = 0.f;

    // geometry_id,hit_id,channel0,channel1,timestamp,value
    DFE_NAMEDTUPLE(surface, geometry_id, cx, cy, cz, rot_xu, rot_xv, rot_xw,
                   rot_yu, rot_yv, rot_yw, rot_zu, rot_zv, rot_zw);

};  // struct surface

}  // namespace traccc::io::csv
