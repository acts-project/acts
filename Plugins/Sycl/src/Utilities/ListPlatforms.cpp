// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// System include(s)
#include <iostream>

// SYCL plugin include(s)
#include "Acts/Plugins/Sycl/Utilities/ListPlatforms.hpp"

// SYCL include
#include <CL/sycl.hpp>

namespace Acts::Sycl {
/// @brief This function allows us to list available SYCL platforms and devices.
void listPlatforms() {
  for (const sycl::platform& platform : sycl::platform::get_platforms()) {
    // Print some information about the platform.
    std::cout << "============ Platform ============" << std::endl;
    std::cout << " Name   : " << platform.get_info<sycl::info::platform::name>()
              << std::endl;
    std::cout << " Vendor : "
              << platform.get_info<sycl::info::platform::vendor>() << std::endl;
    std::cout << " Version: "
              << platform.get_info<sycl::info::platform::version>()
              << std::endl;

    // Loop over all devices available from this platform.
    for (const sycl::device& device : platform.get_devices()) {
      // Print some information about the device.
      std::cout << "------------- Device -------------" << std::endl;
      std::cout << " Name   : " << device.get_info<sycl::info::device::name>()
                << std::endl;
      std::cout << " Vendor : " << device.get_info<sycl::info::device::vendor>()
                << std::endl;
      std::cout << " Version: "
                << device.get_info<sycl::info::device::version>() << std::endl;
    }
  }
}
}  // namespace Acts::Sycl
