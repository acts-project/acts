// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "detray/definitions/detail/cuda_definitions.hpp"
#include "detray/utils/ranges.hpp"

// Detray test include(s)
#include "transform_store_cuda_kernel.hpp"

namespace detray {

__global__ void transform_test_kernel(
    vecmem::data::vector_view<detray::point3> input_data,
    typename host_transform_store_t::view_type store_data,
    vecmem::data::vector_view<detray::point3> output_data) {
  typename device_transform_store_t::context_type ctx0;
  device_transform_store_t store(store_data);

  vecmem::device_vector<detray::point3> input(input_data);
  vecmem::device_vector<detray::point3> output(output_data);

  auto range = detray::ranges::subrange(store.get(ctx0),
                                        dindex_range{0u, store.size(ctx0)});
  output[threadIdx.x] = range[threadIdx.x].point_to_global(input[threadIdx.x]);
}

void transform_test(vecmem::data::vector_view<detray::point3> input_data,
                    typename host_transform_store_t::view_type store_data,
                    vecmem::data::vector_view<detray::point3> output_data,
                    std::size_t n_transforms) {
  int block_dim = 1;
  int thread_dim(n_transforms);

  // run the kernel
  transform_test_kernel<<<block_dim, thread_dim>>>(input_data, store_data,
                                                   output_data);

  // cuda error check
  DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
  DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}
}  // namespace detray
