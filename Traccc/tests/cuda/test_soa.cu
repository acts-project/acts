// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <cstdint>

#include <gtest/gtest.h>
#include <traccc/utils/soa.hpp>
#include <traccc/utils/soa_types/soa_bound_track_parameters.hpp>
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/host_memory_resource.hpp>

using test_layout = traccc::soa_layout<std::uint32_t, float[3], std::uint32_t>;
using test_view = traccc::soa_view<test_layout>;
using test_buffer = traccc::soa_buffer<test_layout>;

__global__ void fillSoAKernel(test_view v) {
  traccc::soa_device<test_layout> d(v);

  unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;

  if (tid < v.size()) {
    d.at<0>(tid) = 1u;
    d.at<1, 0>(tid) = static_cast<float>(tid);
    d.at<1, 1>(tid) = static_cast<float>(tid) + 0.25f;
    d.at<1, 2>(tid) = static_cast<float>(tid) + 0.5f;
    d.at<2>(tid) = 3u;
  }
}

__global__ void testSoAKernel(const test_view v, std::uint32_t* total) {
  const traccc::soa_device<test_layout> d(v);

  __shared__ std::uint32_t block_total;

  if (threadIdx.x == 0) {
    block_total = 0;
  }

  __syncthreads();

  unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;

  std::uint32_t warp_total = d.at<0>(tid) + d.at<2>(tid);

  for (int i = 16; i >= 1; i /= 2) {
    warp_total += __shfl_xor_sync(0xffffffff, warp_total, i, 32);
  }

  if (threadIdx.x % 32 == 0) {
    atomicAdd(&block_total, warp_total);
  }

  __syncthreads();

  if (threadIdx.x == 0) {
    atomicAdd(total, block_total);
  }
}

TEST(CUDASoA, FillAndReduce) {
  vecmem::cuda::device_memory_resource mr;

  std::uint32_t n = 1024u * 1024u;

  test_buffer b(n, mr);

  fillSoAKernel<<<n / 1024u, 1024u>>>(b);

  ASSERT_EQ(cudaPeekAtLastError(), cudaSuccess);
  ASSERT_EQ(cudaDeviceSynchronize(), cudaSuccess);

  std::uint32_t host_result;
  std::uint32_t* dev_result = nullptr;

  ASSERT_EQ(cudaMalloc(&dev_result, sizeof(std::uint32_t)), cudaSuccess);
  ASSERT_EQ(cudaMemset(dev_result, 0, sizeof(std::uint32_t)), cudaSuccess);

  testSoAKernel<<<n / 1024u, 1024u>>>(b, dev_result);

  ASSERT_EQ(cudaPeekAtLastError(), cudaSuccess);
  ASSERT_EQ(cudaDeviceSynchronize(), cudaSuccess);

  ASSERT_EQ(cudaMemcpy(&host_result, dev_result, sizeof(std::uint32_t),
                       cudaMemcpyDeviceToHost),
            cudaSuccess);

  // Every row contributes 1 from the first column and 3 from the last.
  ASSERT_EQ(host_result, 4u * n);

  ASSERT_EQ(cudaFree(dev_result), cudaSuccess);
}

TEST(CUDASoA, ArrayColumns) {
  vecmem::cuda::device_memory_resource device_mr;
  vecmem::host_memory_resource host_mr;

  std::uint32_t n = 1000u;

  test_buffer device_buffer(n, device_mr);
  test_buffer host_buffer(n, host_mr);

  fillSoAKernel<<<(n + 255u) / 256u, 256u>>>(device_buffer);

  ASSERT_EQ(cudaPeekAtLastError(), cudaSuccess);
  ASSERT_EQ(cudaDeviceSynchronize(), cudaSuccess);

  // Two buffers of the same size share one layout, so the raw bytes can be
  // copied as they are.
  ASSERT_EQ(host_buffer.capacity(), device_buffer.capacity());
  const std::size_t bytes = test_layout::total_size *
                            static_cast<std::size_t>(host_buffer.capacity());
  ASSERT_EQ(cudaMemcpy(host_buffer.ptr(), device_buffer.ptr(), bytes,
                       cudaMemcpyDeviceToHost),
            cudaSuccess);

  const traccc::soa_device<test_layout> d(host_buffer);

  for (unsigned int i = 0u; i < n; ++i) {
    ASSERT_EQ(d.at<0>(i), 1u);
    ASSERT_EQ((d.at<1, 0>(i)), static_cast<float>(i));
    ASSERT_EQ((d.at<1, 1>(i)), static_cast<float>(i) + 0.25f);
    ASSERT_EQ((d.at<1, 2>(i)), static_cast<float>(i) + 0.5f);
    ASSERT_EQ(d.at<2>(i), 3u);
  }
}

__global__ void fillBoundSoA(
    traccc::soa_bound_track_parameters_view<traccc::default_algebra> view) {
  const unsigned int row = blockIdx.x * blockDim.x + threadIdx.x;
  if (row >= view.size())
    return;
  traccc::bound_track_parameters<> value;
  for (unsigned int i = 0; i < 6; ++i) {
    value[i] = static_cast<float>(row * 64 + i);
    for (unsigned int j = 0; j < 6; ++j) {
      detray::getter::element(value.covariance(), i, j) =
          static_cast<float>(row * 64 + 6 + i * 6 + j);
    }
  }
  value.set_surface_link(detray::geometry::identifier{}.set_index(row));
  traccc::soa_bound_track_parameters_device<traccc::default_algebra>(view).set(
      row, value);
}

__global__ void readBoundSoA(
    traccc::soa_bound_track_parameters_const_view<traccc::default_algebra> view,
    unsigned int* failures) {
  const unsigned int row = blockIdx.x * blockDim.x + threadIdx.x;
  if (row >= view.size())
    return;
  const auto value =
      traccc::soa_bound_track_parameters_const_device<traccc::default_algebra>(
          view)
          .get(row);
  bool ok = value.surface_link().index() == row;
  for (unsigned int i = 0; i < 6; ++i) {
    ok &= value[i] == static_cast<float>(row * 64 + i);
    for (unsigned int j = 0; j < 6; ++j) {
      ok &= detray::getter::element(value.covariance(), i, j) ==
            static_cast<float>(row * 64 + 6 + i * 6 + j);
    }
  }
  if (!ok)
    atomicAdd(failures, 1u);
}

TEST(CUDASoA, BoundParametersAcrossTileBoundaries) {
  vecmem::cuda::device_memory_resource mr;
  for (const unsigned int n : {1u, 7u, 8u, 9u, 31u, 32u, 33u, 65u, 1000u}) {
    traccc::soa_bound_track_parameters_buffer<traccc::default_algebra> buffer(
        n, mr);
    unsigned int* failures = nullptr;
    ASSERT_EQ(cudaMalloc(&failures, sizeof(unsigned int)), cudaSuccess);
    ASSERT_EQ(cudaMemset(failures, 0, sizeof(unsigned int)), cudaSuccess);
    fillBoundSoA<<<(n + 63) / 64, 64>>>(buffer);
    ASSERT_EQ(cudaPeekAtLastError(), cudaSuccess);
    readBoundSoA<<<(n + 63) / 64, 64>>>(buffer, failures);
    ASSERT_EQ(cudaPeekAtLastError(), cudaSuccess);
    ASSERT_EQ(cudaDeviceSynchronize(), cudaSuccess);
    unsigned int host_failures = 0;
    ASSERT_EQ(cudaMemcpy(&host_failures, failures, sizeof(unsigned int),
                         cudaMemcpyDeviceToHost),
              cudaSuccess);
    EXPECT_EQ(host_failures, 0u) << "row count " << n;
    ASSERT_EQ(cudaFree(failures), cudaSuccess);
  }
}
