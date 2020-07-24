// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// System include(s).
#include <cassert>

namespace Acts {
namespace Cuda {

template <typename T>
void copyToDevice(DeviceMatrix<T>& device, const HostMatrix<T>& host) {

  // Some security check(s).
  assert(host.nRows() == device.nRows());
  assert(host.nCols() == device.nCols());

  // Perform the copy.
  device.copyFrom(host.getPtr(), device.size(), 0);
  return;
}

template <typename T>
void copyToDevice(DeviceMatrix<T>& device, const HostMatrix<T>& host,
                  cudaStream_t stream) {

  // Some security check(s).
  assert(host.nRows() == device.nRows());
  assert(host.nCols() == device.nCols());

  // Perform the copy.
  device.copyFrom(host.getPtr(), device.size(), 0, stream);
  return;
}

template <typename T>
void copyToHost(HostMatrix<T>& host, const DeviceMatrix<T>& device) {

  // Some security check(s).
  assert(host.nRows() == device.nRows());
  assert(host.nCols() == device.nCols());

  // Perform the copy.
  host.copyFrom(device.getPtr(), host.size(), 0);
  return;
}

template <typename T>
void copyToHost(HostMatrix<T>& host, const DeviceMatrix<T>& device,
                cudaStream_t stream) {

  // Some security check(s).
  assert(host.nRows() == device.nRows());
  assert(host.nCols() == device.nCols());

  // Perform the copy.
  host.copyFrom(device.getPtr(), host.size(), 0, stream);
  return;
}

template <typename T>
void copyToDevice(DeviceVector<T>& device, const HostVector<T>& host) {

  // Some security check(s).
  assert(host.size() == device.size());

  // Perform the copy.
  device.copyFrom(host.getPtr(), device.size(), 0);
  return;
}

template <typename T>
void copyToDevice(DeviceVector<T>& device, const HostVector<T>& host,
                  cudaStream_t stream) {

  // Some security check(s).
  assert(host.size() == device.size());

  // Perform the copy.
  device.copyFrom(host.getPtr(), device.size(), 0, stream);
  return;
}

template <typename T>
void copyToHost(HostVector<T>& host, const DeviceVector<T>& device) {

  // Some security check(s).
  assert(host.size() == device.size());

  // Perform the copy.
  host.copyFrom(device.getPtr(), host.size(), 0);
  return;
}

template <typename T>
void copyToHost(HostVector<T>& host, const DeviceVector<T>& device,
                cudaStream_t stream) {

  // Some security check(s).
  assert(host.size() == device.size());

  // Perform the copy.
  host.copyFrom(device.getPtr(), host.size(), 0, stream);
  return;
}

} // namespace Cuda
} // namespace Acts
