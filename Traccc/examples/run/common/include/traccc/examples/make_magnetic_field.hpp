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
#include "traccc/options/magnetic_field.hpp"

namespace traccc::details {

/// Create a magnetic field object based on the provided options
///
/// @param opts The command line options for the magnetic field
/// @return A magnetic field object configured according to the options
///
magnetic_field make_magnetic_field(const opts::magnetic_field& opts);

}  // namespace traccc::details
