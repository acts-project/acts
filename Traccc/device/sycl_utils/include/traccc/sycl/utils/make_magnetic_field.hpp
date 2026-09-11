/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/sycl/utils/queue_wrapper.hpp"

// Project include(s).
#include "traccc/bfield/magnetic_field.hpp"

namespace traccc::sycl {

/// Create a magnetic field usable on the selected SYCL device
///
/// @param bfield The magnetic field to be copied
/// @param queue The SYCL queue to use for the operation
/// @return A copy of the magnetic field that can be used on the selected SYCL
///         device
///
magnetic_field make_magnetic_field(const magnetic_field& bfield,
                                   queue_wrapper& queue);

}  // namespace traccc::sycl
