// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "detray/definitions/detail/cuda_definitions.hpp"
#include "detray/utils/ranges.hpp"

// Detray test include(s)
#include "sf_finders_grid_cuda_kernel.hpp"

// Vecmem include(s)
#include <vecmem/containers/device_vector.hpp>

namespace detray {

//----------------------------------------------------
//  test function for grid data with replace populator
//----------------------------------------------------

/// cuda kernel for grid_replace_test
__global__ void grid_replace_test_kernel(
    host_grid3_single::view_type grid_view) {
  // Let's try building the grid object
  device_grid3_single g3_device(grid_view);

  // Get axes on the device-side
  const auto& axis_x = g3_device.template get_axis<axis::label::e_x>();
  const auto& axis_y = g3_device.template get_axis<axis::label::e_y>();
  const auto& axis_z = g3_device.template get_axis<axis::label::e_z>();

  dindex gid = g3_device.serialize(
      detray::axis::multi_bin<3>{threadIdx.x, threadIdx.y, threadIdx.z});

  point3 tp{axis_x.min() + gid * axis_x.bin_width(),
            axis_y.min() + gid * axis_y.bin_width(),
            axis_z.min() + gid * axis_z.bin_width()};

  // replace the bin elements
  g3_device.template populate<replace<>>(gid, std::move(tp));
}

/// grid_replace_test implementation
void grid_replace_test(host_grid3_single::view_type grid_view,
                       std::size_t dim_x, std::size_t dim_y,
                       std::size_t dim_z) {
  int n_blocks = 1;
  dim3 n_threads(dim_x, dim_y, dim_z);

  // run the kernel
  grid_replace_test_kernel<<<n_blocks, n_threads>>>(grid_view);

  // cuda error check
  DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
  DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

/// cuda kernel for grid_replace_ci_test
__global__ void grid_replace_ci_test_kernel(
    host_grid2_single_ci::view_type grid_view) {
  // Let's try building the grid object
  device_grid2_single_ci g2_device(grid_view);

  // Get axes on the device-side
  const auto& axis_r = g2_device.template get_axis<axis::label::e_r>();
  const auto& axis_phi = g2_device.template get_axis<axis::label::e_phi>();

  auto gid = threadIdx.x + threadIdx.y * blockDim.x;

  point3 tp{axis_r.min() + gid * axis_r.bin_width(threadIdx.x),
            axis_phi.min() + gid * axis_phi.bin_width(), 0.5f};

  // replace the bin elements
  g2_device.template populate<replace<>>(gid, std::move(tp));
}

// test function for replace populator with circular and irregular axis
void grid_replace_ci_test(host_grid2_single_ci::view_type grid_view,
                          std::size_t dim_x, std::size_t dim_y) {
  int n_blocks = 1;
  dim3 n_threads(dim_x, dim_y);

  // run the kernel
  grid_replace_ci_test_kernel<<<n_blocks, n_threads>>>(grid_view);

  // cuda error check
  DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
  DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

//----------------------------------------------------
// test function for grid data with complete populator
//----------------------------------------------------

// cuda kernel for grid_complete_test
__global__ void grid_complete_kernel(host_grid2_array::view_type grid_view) {
  // Let's try building the grid object
  device_grid2_array g2_device(grid_view);

  // Get axes on the device-side
  const auto& axis_r = g2_device.template get_axis<axis::label::e_r>();
  const auto& axis_phi = g2_device.template get_axis<axis::label::e_phi>();

  auto gid = threadIdx.x + threadIdx.y * blockDim.x;
  auto tp = point3{axis_r.min() + gid * axis_r.bin_width(),
                   axis_phi.min() + gid * axis_phi.bin_width(), 0.5f};

  g2_device.template populate<complete<>>(gid, std::move(tp));
}

// grid_complete_test implementation
void grid_complete_test(host_grid2_array::view_type grid_view,
                        std::size_t dim_x, std::size_t dim_y) {
  int block_dim = 1;
  dim3 thread_dim(dim_x, dim_y);

  // run the kernel
  grid_complete_kernel<<<block_dim, thread_dim>>>(grid_view);

  // cuda error check
  DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
  DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

//--------------------------------------------------
// test function for grid data with attach populator
//--------------------------------------------------

// cuda kernel for grid_attach_test
__global__ void grid_attach_kernel(host_grid2_array::view_type grid_view) {
  // Let's try building the grid object
  device_grid2_array g2_device(grid_view);

  // Get axes on the device-side
  const auto& axis_r = g2_device.template get_axis<axis::label::e_r>();
  const auto& axis_phi = g2_device.template get_axis<axis::label::e_phi>();

  auto width_r = axis_r.m_binning.bin_width();
  auto width_phi = axis_phi.m_binning.bin_width();

  auto gid = threadIdx.x + threadIdx.y * blockDim.x;
  auto tp = point3{axis_r.min() + gid * width_r,
                   axis_phi.min() + gid * width_phi, 0.5f};

  g2_device.template populate<attach<>>(gid, std::move(tp));
}

// grid_attach_test implementation
void grid_attach_test(host_grid2_array::view_type grid_view, std::size_t dim_x,
                      std::size_t dim_y) {
  int block_dim = 1;
  dim3 thread_dim(dim_x, dim_y);

  // run the kernel
  grid_attach_kernel<<<block_dim, thread_dim>>>(grid_view);

  // cuda error check
  DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
  DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

// cuda kernel for grid_dynamic_attach_test
__global__ void grid_dynamic_attach_kernel(
    host_grid2_dynamic_array::view_type grid_view) {
  // Let's try building the grid object
  device_grid2_dynamic_array g2_device(grid_view);

  // Get axes on the device-side
  const auto& axis_r = g2_device.template get_axis<axis::label::e_r>();
  const auto& axis_phi = g2_device.template get_axis<axis::label::e_phi>();

  auto width_r = axis_r.m_binning.bin_width();
  auto width_phi = axis_phi.m_binning.bin_width();

  auto gid = threadIdx.x + threadIdx.y * blockDim.x;
  auto tp = point3{axis_r.min() + gid * width_r,
                   axis_phi.min() + gid * width_phi, 0.5f};

  g2_device.template populate<attach<>>(gid, std::move(tp));
}

// grid_dynamic_attach_test implementation
void grid_dynamic_attach_test(host_grid2_dynamic_array::view_type grid_view,
                              std::size_t dim_x, std::size_t dim_y) {
  int block_dim = 1;
  dim3 thread_dim(dim_x, dim_y);

  // run the kernel
  grid_dynamic_attach_kernel<<<block_dim, thread_dim>>>(grid_view);

  // cuda error check
  DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
  DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

//----------------------------------------------------
// Device side grid reader for debugging
//----------------------------------------------------

// cuda kernel for attach_read_test
template <typename device_grid_t, typename view_t>
__global__ void print_grid_kernel(view_t grid_view) {
  // Let's try building the grid object
  device_grid_t g2_device(grid_view);

  axis::multi_bin<device_grid_t::dim> mbin;
  if constexpr (device_grid_t::dim == 2) {
    mbin = {threadIdx.x, threadIdx.y};
  } else {
    mbin = {threadIdx.x, threadIdx.y, threadIdx.z};
  }

  for (auto& pt : g2_device.bin(mbin)) {
    printf("[%f %f %f]\n", pt[0], pt[1], pt[2]);
  }
}

// grid_attach_read_test implementation
template <typename device_grid_t, typename view_t, typename... I>
void print_grid(view_t grid_view, I... dims) {
  int block_dim = 1;
  dim3 thread_dim(dims...);

  // run the kernel
  print_grid_kernel<device_grid_t, view_t>
      <<<block_dim, thread_dim>>>(grid_view);

  // cuda error check
  DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
  DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

// Explicit instantioations
template void print_grid<device_grid3_single>(host_grid3_single::view_type,
                                              dindex, dindex, dindex);

template void print_grid<device_grid2_single_ci>(
    host_grid2_single_ci::view_type, dindex, dindex);

template void print_grid<device_grid2_array>(host_grid2_array::view_type,
                                             dindex, dindex);

template void print_grid<device_grid2_dynamic_array>(
    device_grid2_dynamic_array::view_type, dindex, dindex);

//---------------------------------------
//  test function for collection of grids
//---------------------------------------

/// cuda kernel for grid_collection_test
__global__ void grid_collection_test_kernel(
    grid_collection<n_own_host_grid3_array>::view_type grid_coll_view,
    vecmem::data::vector_view<dindex> n_bins_view,
    vecmem::data::vector_view<std::array<dindex, 3>> result_bins_view) {
  // Let's try building the grid object
  grid_collection<n_own_device_grid3_array> device_coll(grid_coll_view);
  vecmem::device_vector<dindex> n_bins(n_bins_view);
  vecmem::device_vector<std::array<dindex, 3>> result_bins(result_bins_view);

  // test the grid axes of the second grid in the collection
  if (threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0) {
    const auto& axis_r =
        device_coll[blockIdx.x].template get_axis<axis::label::e_r>();
    const auto& axis_phi =
        device_coll[blockIdx.x].template get_axis<axis::label::e_phi>();
    const auto& axis_z =
        device_coll[blockIdx.x].template get_axis<axis::label::e_z>();

    n_bins[0 + blockIdx.x * 3] = axis_r.nbins();
    n_bins[1 + blockIdx.x * 3] = axis_phi.nbins();
    n_bins[2 + blockIdx.x * 3] = axis_z.nbins();
  }

  // Read the entire grid content
  int gid = threadIdx.z * blockDim.y * blockDim.x + threadIdx.y * blockDim.x +
            threadIdx.x;
  if (gid < device_coll[blockIdx.x].nbins()) {
    for (const auto [i, bin_entry] :
         detray::views::enumerate(device_coll[blockIdx.x].bin(gid))) {
      result_bins[gid + device_coll.offsets()[blockIdx.x]][i] = bin_entry;
    }
  }
}

/// grid_collection_test implementation
void grid_collection_test(
    grid_collection<n_own_host_grid3_array>::view_type grid_coll_view,
    vecmem::data::vector_view<dindex> n_bins_view,
    vecmem::data::vector_view<std::array<dindex, 3>> result_bins_view,
    std::size_t n_grids, std::size_t dim_x, std::size_t dim_y,
    std::size_t dim_z) {
  int n_blocks = n_grids;
  dim3 n_threads(dim_x, dim_y, dim_z);

  // run the kernel
  grid_collection_test_kernel<<<n_blocks, n_threads>>>(
      grid_coll_view, n_bins_view, result_bins_view);

  // cuda error check
  DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
  DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}
}  // namespace detray
