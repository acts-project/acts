// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// SYCL plugin include(s)
#include "Acts/Plugins/Sycl/Utilities/DeviceSelector.hpp"

// SYCL include
#include <CL/sycl.hpp>

namespace Acts::Sycl {
DeviceSelector::DeviceSelector(const std::string& deviceName)
    : m_defaultSelector(cl::sycl::default_selector()),
      m_deviceName(deviceName) {}

int DeviceSelector::operator()(const cl::sycl::device& d) const {
  // Under no circumstances do we accept any NVidia OpenCL devices.
  const std::string vendor = d.get_info<cl::sycl::info::device::vendor>();
  const std::string version = d.get_info<cl::sycl::info::device::version>();
  if ((vendor.find("NVIDIA") != std::string::npos) &&
      (version.find("OpenCL") != std::string::npos)) {
    return -1;
  }

  // If the user provided a substring of the device name, look for that device
  if (!m_deviceName.empty()) {
    if (d.get_info<cl::sycl::info::device::name>().find(m_deviceName) !=
        std::string::npos) {
      return 1;
    } else {
      return -1;
    }
  }

  // Otherwise return the value defined by the default selector
  return m_defaultSelector(d);
};
}  // namespace Acts::Sycl
