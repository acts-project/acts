// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

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
