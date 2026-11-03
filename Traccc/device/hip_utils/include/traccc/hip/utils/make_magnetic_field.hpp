// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

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
