/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "traccc/definitions/common.hpp"

namespace traccc {

struct track_matching_config {
    float matching_ratio = 0.5f;

    bool double_matching = true;
};

}  // namespace traccc
