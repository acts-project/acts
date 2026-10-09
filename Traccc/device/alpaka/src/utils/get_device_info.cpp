// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "utils.hpp"

// Project include(s).
#include "traccc/alpaka/utils/get_device_info.hpp"

namespace traccc::alpaka {

std::string get_device_info() {
  int device = 0;
  auto devAcc = ::alpaka::getDevByIdx(::alpaka::Platform<Acc>{}, 0u);
  return std::string("Using Alpaka device: " + ::alpaka::getName(devAcc) +
                     " [id: " + std::to_string(device) + "] ");
}

}  // namespace traccc::alpaka
