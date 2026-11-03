// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <string>

namespace traccc::alpaka {

/// Function that prints the current device information to the console.
/// Included as part of the traccc::alpaka namespace, to avoid having to include
/// alpaka headers in any users of the library.
std::string get_device_info();

}  // namespace traccc::alpaka
