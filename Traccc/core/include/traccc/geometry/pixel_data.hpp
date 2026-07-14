/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/primitives.hpp"
#include "traccc/definitions/qualifiers.hpp"

namespace traccc {

/// A very basic pixel segmentation with
/// a minimum corner and ptich x/y
///
/// No checking on out of bounds done
struct pixel_data {

    scalar min_corner_x = 0.f;
    scalar min_corner_y = 0.f;
    scalar pitch_x = 1.f;
    scalar pitch_y = 1.f;
    char dimension = 2;

    TRACCC_HOST_DEVICE
    vector2 get_pitch() const { return {pitch_x, pitch_y}; };
};

}  // namespace traccc
