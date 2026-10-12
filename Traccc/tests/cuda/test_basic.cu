// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <gtest/gtest.h>

TEST(CUDABasic, DeviceCount) {
  int nDev = -1;

  ASSERT_EQ(cudaGetDeviceCount(&nDev), cudaSuccess);

  ASSERT_GE(nDev, 1);
}

TEST(CUDABasic, Memory) {
  void *ptr = nullptr;

  ASSERT_EQ(cudaMalloc(&ptr, 1024), cudaSuccess);

  ASSERT_NE(ptr, nullptr);

  ASSERT_EQ(cudaFree(ptr), cudaSuccess);
}

__global__ void testKernel(int *output) {
  *output = 0x0BADF00D;  // This test sponsored by R1.
}

TEST(CUDABasic, LaunchKernel) {
  int *ptr = nullptr;
  int val = 0;

  ASSERT_EQ(cudaMalloc(&ptr, sizeof(int)), cudaSuccess);

  ASSERT_NE(ptr, nullptr);

  testKernel<<<1, 1>>>(ptr);

  ASSERT_EQ(cudaPeekAtLastError(), cudaSuccess);

  ASSERT_EQ(cudaMemcpy(&val, ptr, sizeof(int), cudaMemcpyDeviceToHost),
            cudaSuccess);

  ASSERT_EQ(val, 0x0BADF00D);

  ASSERT_EQ(cudaFree(ptr), cudaSuccess);
}
