// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Built with NDEBUG defined: the device-side assert below only fires if
// ACTS_FORCE_ASSERTIONS re-enabled it.

#include <cassert>
#include <cstdio>

#include <cuda_runtime.h>

__global__ void assertKernel() {
  assert(false && "forced assertion");
}

int main() {
  int nDevices = 0;
  if (cudaGetDeviceCount(&nDevices) != cudaSuccess || nDevices == 0) {
    // Reported as skipped by CTest
    return 77;
  }
  assertKernel<<<1, 1>>>();
  if (cudaDeviceSynchronize() == cudaErrorAssert) {
    // The device runtime has printed the assertion message by now
    return 1;
  }
  std::puts("ACTS_ASSERTION_NOT_TRIGGERED");
  return 0;
}
