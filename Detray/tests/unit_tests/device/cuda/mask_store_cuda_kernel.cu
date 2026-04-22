// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "detray/definitions/detail/cuda_definitions.hpp"

// Detray test include(s)
#include "mask_store_cuda_kernel.hpp"

namespace detray {

/// test kernel function to fill the output vector with is_inside function
/// return values
__global__ void mask_test_kernel(
    typename host_store_type::view_type store_data,
    vecmem::data::vector_view<point3> input_point3_data,
    vecmem::data::jagged_vector_view<int> output_data) {
  // get mask store
  device_store_type store(store_data);

  // get mask objects
  vecmem::device_vector<point3> input_point3(input_point3_data);
  vecmem::jagged_device_vector<int> output_device(output_data);

  const auto& rectangle_mask = store.get<mask_id::e_rectangle2D>()[0];
  const auto& trapezoid_mask = store.get<mask_id::e_trapezoid2D>()[0];
  const auto& ring_mask = store.get<mask_id::e_ring2D>()[0];
  const auto& cylinder_mask = store.get<mask_id::e_cylinder2D>()[0];
  const auto& annulus_mask = store.get<mask_id::e_annulus2D>()[0];

  // get device results from is_inside function
  for (int i = 0; i < n_points; i++) {
    output_device[0].push_back(rectangle_mask.is_inside(input_point3[i]));
    output_device[1].push_back(trapezoid_mask.is_inside(input_point3[i]));
    output_device[2].push_back(ring_mask.is_inside(input_point3[i]));
    output_device[3].push_back(cylinder_mask.is_inside(input_point3[i]));
    output_device[4].push_back(annulus_mask.is_inside(input_point3[i]));
  }
}

void mask_test(typename host_store_type::view_type store_data,
               vecmem::data::vector_view<point3> input_point3_data,
               vecmem::data::jagged_vector_view<int> output_data) {
  int block_dim = 1;
  int thread_dim = 1;

  // run the test kernel
  mask_test_kernel<<<block_dim, thread_dim>>>(store_data, input_point3_data,
                                              output_data);

  // cuda error check
  DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
  DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

}  // namespace detray
