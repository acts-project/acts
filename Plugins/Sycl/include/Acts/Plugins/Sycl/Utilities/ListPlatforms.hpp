// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

namespace Acts::Sycl {
/// @brief This function allows us to list available SYCL platforms and devices.
///
/// Available platforms and devices only include previously linked targets by
/// CMake, which can optionally be altered by environment variable SYCL_BE.
void listPlatforms();

}  // namespace Acts::Sycl
