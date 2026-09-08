// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "detray/definitions/detail/cuda_definitions.hpp"

// Detray test include(s)
#include "known_substructure_matrix_cuda_kernel.hpp"

namespace detray {

__global__ void ksm_kernel(double *out) {
  ksm_cuda_test::compute(out);
}

void ksm_run_on_device(double *out) {
  // A single thread is enough: this checks that the operations compile and run
  // on device and agree with the host, not that they are parallel.
  ksm_kernel<<<1, 1>>>(out);

  // cuda error check
  DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
  DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

}  // namespace detray
