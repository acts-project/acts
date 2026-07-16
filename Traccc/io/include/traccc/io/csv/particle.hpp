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

/// Type used in reading CSV particle data into memory
struct particle {

    uint64_t particle_id = 0;
    int particle_type = 0;
    int process = 0;
    float vx = 0.f;
    float vy = 0.f;
    float vz = 0.f;
    float vt = 0.f;
    float px = 0.f;
    float py = 0.f;
    float pz = 0.f;
    float m = 0.f;
    float q = 0.f;

    DFE_NAMEDTUPLE(particle, particle_id, particle_type, process, vx, vy, vz,
                   vt, px, py, pz, m, q);
};

}  // namespace traccc::io::csv
