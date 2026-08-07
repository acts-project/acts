/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/bfield/magnetic_field.hpp"

namespace traccc::hip {

/// Storage method for inhomogeneous magnetic fields
enum class magnetic_field_storage {
  global_memory  ///< Store the magnetic field in global device memory
};

/// Create a magnetic field usable on the active HIP device
///
/// @param bfield The magnetic field to be copied
/// @param storage The storage method to use for the magnetic field
/// @return A copy of the magnetic field that can be used on the active HIP
///         device
///
magnetic_field make_magnetic_field(
    const magnetic_field& bfield,
    magnetic_field_storage storage = magnetic_field_storage::global_memory);

}  // namespace traccc::hip
