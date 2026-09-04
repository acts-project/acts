// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "detector_construction.hpp"
#include "detray/definitions/detail/cuda_definitions.hpp"
#include "detray/utils/logging.hpp"

namespace detray::tutorial {

/// Kernel that runs the entire propagation loop
__global__ void print_kernel(
    typename detray::tutorial::detector_host_t::view_type det_data) {
  int gid = threadIdx.x + blockIdx.x * blockDim.x;

  if (gid > 0) {
    return;
  }

  // Setup of the device-side detector
  detray::tutorial::detector_device_t det(det_data);

  DETRAY_INFO_DEVICE("Number of volumes: %d", det.volumes().size());
  DETRAY_INFO_DEVICE("Number of transforms: %d", det.transform_store().size());
  DETRAY_INFO_DEVICE("First translation: {%f,%f,%f}",
                     det.transform_store().at(0).translation()[0],
                     det.transform_store().at(0).translation()[1],
                     det.transform_store().at(0).translation()[2]);
  DETRAY_INFO_DEVICE("Number of rectangles: %d",
                     det.mask_store().get<mask_id::e_rectangle2D>().size());
  DETRAY_INFO_DEVICE("Number of trapezoids: %d",
                     det.mask_store().get<mask_id::e_trapezoid2D>().size());
  DETRAY_INFO_DEVICE("Number of portal discs: %d",
                     det.mask_store().get<mask_id::e_ring2D>().size());
  DETRAY_INFO_DEVICE(
      "Number of portal cylinders: %d",
      det.mask_store().get<mask_id::e_concentric_cylinder2D>().size());
  DETRAY_INFO_DEVICE(
      "Number of portal collections: %d",
      det.accelerator_store().get<acc_id::e_surface_brute_force>().size());
  DETRAY_INFO_DEVICE(
      "Number of disc grids: %d",
      det.accelerator_store().get<acc_id::e_surface_ring2D_grid>().size());
  DETRAY_INFO_DEVICE("Number of cylinder grids: %d",
                     det.accelerator_store()
                         .get<acc_id::e_surface_concentric_cylinder2D_grid>()
                         .size());
}

void print(typename detray::tutorial::detector_host_t::view_type det_data) {
  // run the tutorial kernel
  print_kernel<<<1, 1>>>(det_data);

  // cuda error check
  DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
  DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

}  // namespace detray::tutorial
