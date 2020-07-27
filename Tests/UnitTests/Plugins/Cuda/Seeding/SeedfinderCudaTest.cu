// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Local include(s).
#include "SeedfinderCudaTest.hpp"

// CUDA plugin include(s).
#include "Acts/Plugins/Cuda/Utilities/ErrorCheck.cuh"

// CUDA include(s).
#include <cuda_runtime.h>

// System include(s).
#include <iostream>

void setupCudaDevice(int deviceID, int& maxThreadsPerBlock) {

  // Force the test to use the selected device.
  ACTS_CUDA_ERROR_CHECK(cudaSetDevice(deviceID));

  // Get the properties of the selected device.
  cudaDeviceProp prop;
  ACTS_CUDA_ERROR_CHECK(cudaGetDeviceProperties(&prop, deviceID));
  std::cout << "GPU Device " << deviceID << ": \"" << prop.name
            << "\" with compute capability " << prop.major << "." << prop.minor
            << std::endl;

  // Set up the configuration object based on the device's properties.
  maxThreadsPerBlock = prop.maxThreadsPerBlock;
  return;
}
