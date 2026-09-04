/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/alpaka/utils/queue.hpp"
#include "traccc/bfield/magnetic_field.hpp"

namespace traccc::alpaka {

/// Create a magnetic field usable on the active device
///
/// @param bfield The magnetic field to be copied
//
magnetic_field make_magnetic_field(const magnetic_field& bfield,
                                   const queue& queue);

}  // namespace traccc::alpaka
