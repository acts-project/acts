/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "traccc/definitions/common.hpp"

namespace traccc {

struct truth_matching_config {
    float pT_min = 0.5f * traccc::unit<float>::GeV;

    float z_min = -500.f * traccc::unit<float>::mm;
    float z_max = 500.f * traccc::unit<float>::mm;
    float r_max = 200.f * traccc::unit<float>::mm;
    float eta_max = 3.f;
    int process_id = -1;

    unsigned int min_track_candidates = 3;
};

}  // namespace traccc
