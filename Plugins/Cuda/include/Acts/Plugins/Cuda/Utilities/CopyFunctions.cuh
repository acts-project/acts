// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// CUDA plugin include(s).
#include "Acts/Plugins/Cuda/Utilities/DeviceMatrix.cuh"
#include "Acts/Plugins/Cuda/Utilities/DeviceVector.cuh"
#include "Acts/Plugins/Cuda/Utilities/HostMatrix.cuh"
#include "Acts/Plugins/Cuda/Utilities/HostVector.cuh"

// CUDA include(s).
#include <cuda.h>

namespace Acts {
namespace Cuda {

  /// Copy the contents of a matrix from the host to the/a device
  template <typename T>
  void copyToDevice(DeviceMatrix<T>& device, const HostMatrix<T>& host);
  /// Copy the contents of a matrix from the host to the/a device asynchronously
  template <typename T>
  void copyToDevice(DeviceMatrix<T>& device, const HostMatrix<T>& host,
                    cudaStream_t stream);

  /// Copy the contents of a matrix from the/a device to the host
  template <typename T>
  void copyToHost(HostMatrix<T>& host, const DeviceMatrix<T>& device);
  /// Copy the contents of a matrix from the/a device to the host asynchronously
  template <typename T>
  void copyToHost(HostMatrix<T>& host, const DeviceMatrix<T>& device,
                  cudaStream_t stream);

  /// Copy the contents of a vector from the host to the/a device
  template <typename T>
  void copyToDevice(DeviceVector<T>& device, const HostVector<T>& host);
  /// Copy the contents of a vector from the host to the/a device asynchronously
  template <typename T>
  void copyToDevice(DeviceVector<T>& device, const HostVector<T>& host,
                    cudaStream_t stream);

  /// Copy the contents of a vector from the/a device to the host
  template <typename T>
  void copyToHost(HostVector<T>& host, const DeviceVector<T>& device);
  /// Copy the contents of a vector from the/a device to the host asynchronously
  template <typename T>
  void copyToHost(HostVector<T>& host, const DeviceVector<T>& device,
                  cudaStream_t stream);

} // namespace Cuda
} // namespace Acts

// Include the template implementation.
#include "CopyFunctions.ipp"
