// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// CUDA plugin include(s).
#include "Acts/Plugins/Cuda/Utilities/CopyFunctions.hpp"
#include "Acts/Plugins/Cuda/Utilities/DeviceMatrix.hpp"
#include "Acts/Plugins/Cuda/Utilities/ErrorCheck.cuh"
#include "Acts/Plugins/Cuda/Utilities/HostMatrix.hpp"

// Boost include(s).
#include <boost/test/unit_test.hpp>

// System include(s).
#include <cmath>

namespace Acts {
namespace Cuda {
namespace Test {

/// Simple kernel performing a multiplication.
__global__ void matrixTransform(const int size, const float* input,
                                float* output) {
  // Get the index to work on.
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= size) {
    return;
  }

  // Perform the transformation.
  output[i] = input[i] * 2.0f;
  return;
}

BOOST_AUTO_TEST_SUITE(Utilities)
BOOST_AUTO_TEST_CASE(MatrixCopy) {
  // The matrix size to use in the test.
  static constexpr int MATRIX_SIZE = 1000;

  // Create an input matrix.
  HostMatrix<float> inputHost(MATRIX_SIZE, MATRIX_SIZE);
  BOOST_TEST_REQUIRE(inputHost.size() == MATRIX_SIZE * MATRIX_SIZE);
  for (int i = 0; i < MATRIX_SIZE; ++i) {
    for (int j = 0; j < MATRIX_SIZE; ++j) {
      inputHost.set(i, j, M_PI);
    }
  }

  // Create an output matrix.
  HostMatrix<float> outputHost(MATRIX_SIZE, MATRIX_SIZE);
  BOOST_TEST_REQUIRE(outputHost.size() == MATRIX_SIZE * MATRIX_SIZE);

  // Copy the input to the/a device, and run a simple kernel on it.
  DeviceMatrix<float> inputDevice(MATRIX_SIZE, MATRIX_SIZE);
  copyToDevice(inputDevice, inputHost);
  DeviceMatrix<float> outputDevice(MATRIX_SIZE, MATRIX_SIZE);
  static constexpr int blockSize = 256;
  static const int numBlocks = (inputDevice.size() + blockSize - 1) / blockSize;
  matrixTransform<<<numBlocks, blockSize>>>(
      inputDevice.size(), inputDevice.getPtr(), outputDevice.getPtr());
  ACTS_CUDA_ERROR_CHECK(cudaDeviceSynchronize());

  // Copy the result back to the host, and check it.
  copyToHost(outputHost, outputDevice);
  float maxDeviation = 0.0f;
  static constexpr float EXPECTED_RESULT = M_PI * 2.0f;
  for (int i = 0; i < MATRIX_SIZE; ++i) {
    for (int j = 0; j < MATRIX_SIZE; ++j) {
      maxDeviation = std::max(maxDeviation,
                              std::abs(outputHost.get(i, j) - EXPECTED_RESULT));
    }
  }
  BOOST_TEST_REQUIRE(maxDeviation < 0.001);
}
BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test
}  // namespace Cuda
}  // namespace Acts
