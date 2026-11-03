// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <cuda_runtime.h>
#include <gtest/gtest.h>

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);

  int r = RUN_ALL_TESTS();

  cudaError_t cErr = cudaDeviceReset();

  if (cErr == cudaSuccess || cErr == cudaErrorNoDevice ||
      cErr == cudaErrorInsufficientDriver) {
    return r;
  } else {
    return static_cast<int>(cErr);
  }
}
