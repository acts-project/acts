/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Turn off exessive detector building logs (only works with gcc!)
// clang-format off
#include <detray/utils/quiet_log_start.hpp>
// clang-format on

// Detray include(s).
#include <detray/test/common/build_telescope_detector.hpp>
#include <detray/test/common/build_toy_detector.hpp>
#include <detray/test/common/build_wire_chamber.hpp>

// clang-format off
#include <detray/utils/quiet_log_end.hpp>
// clang-format on

// useful for telescope detector creation
#include <detray/geometry/mask.hpp>
#include <detray/geometry/shapes/rectangle2D.hpp>
#include <detray/tracks/ray.hpp>
