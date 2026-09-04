// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "detray/definitions/detail/cuda_definitions.hpp"

// Detray test include(s)
#include "utils_ranges_cuda_kernel.hpp"

namespace detray {

//
// single
//
__global__ void single_kernel(const dindex value, dindex* result) {
  // single view should only add the value 'i' once
  for (auto i : detray::views::single(value)) {
    *result += i;
  }
}

void test_single(const dindex value, dindex& check) {
  dindex* result{nullptr};
  cudaMallocManaged(&result, sizeof(dindex));
  *result = 0u;

  // run the kernel
  single_kernel<<<1, 1>>>(value, result);

  // cuda error check
  DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
  DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());

  check = *result;
  cudaFree(result);
}

//
// pointer
//
__global__ void pointer_kernel(const dindex value, dindex* result) {
  // pointer view should only add the value 'i' once
  for (auto i : detray::views::pointer(value)) {
    *result += i;
  }
}

void test_pointer(const dindex value, dindex& check) {
  dindex* result{nullptr};
  cudaMallocManaged(&result, sizeof(dindex));
  *result = 0u;

  // run the kernel
  pointer_kernel<<<1, 1>>>(value, result);

  // cuda error check
  DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
  DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());

  check = *result;
  cudaFree(result);
}

//
// iota
//
__global__ void iota_kernel(const darray<dindex, 2> range,
                            vecmem::data::vector_view<dindex> check_data) {
  vecmem::device_vector<dindex> check(check_data);

  for (auto i : detray::views::iota(range)) {
    check.push_back(i);
  }
}

void test_iota(const darray<dindex, 2> range,
               vecmem::data::vector_view<dindex> check_data) {
  // run the kernel
  iota_kernel<<<1, 1>>>(range, check_data);

  // cuda error check
  DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
  DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

//
// cartesian product
//
__global__ void cartesian_product_kernel(
    const darray<dindex, 2> range1, const darray<dindex, 2> range2,
    const darray<dindex, 2> range3,
    vecmem::data::vector_view<std::tuple<dindex, dindex, dindex>> check_data) {
  vecmem::device_vector<std::tuple<dindex, dindex, dindex>> check(check_data);

  auto seq1 = detray::views::iota(range1);
  auto seq2 = detray::views::iota(range2);
  auto seq3 = detray::views::iota(range3);

  for (const auto [i, j, k] : detray::views::cartesian_product(
           std::move(seq1), std::move(seq2), std::move(seq3))) {
    check.emplace_back(i, j, k);
  }
}

void test_cartesian_product(
    const darray<dindex, 2> range1, const darray<dindex, 2> range2,
    const darray<dindex, 2> range3,
    vecmem::data::vector_view<std::tuple<dindex, dindex, dindex>> check_data) {
  // run the kernel
  cartesian_product_kernel<<<1, 1>>>(range1, range2, range3, check_data);

  // cuda error check
  DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
  DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

//
// enumerate
//
__global__ void enumerate_kernel(
    vecmem::data::vector_view<uint_holder> seq_data,
    vecmem::data::vector_view<dindex> check_idx_data,
    vecmem::data::vector_view<dindex> check_value_data) {
  vecmem::device_vector<uint_holder> seq(seq_data);
  vecmem::device_vector<dindex> check_idx(check_idx_data);
  vecmem::device_vector<dindex> check_value(check_value_data);

  for (auto [i, v] : detray::views::enumerate(seq)) {
    check_idx.push_back(i);
    check_value.push_back(v.ui);
  }
}

void test_enumerate(vecmem::data::vector_view<uint_holder> seq_data,
                    vecmem::data::vector_view<dindex> check_idx_data,
                    vecmem::data::vector_view<dindex> check_value_data) {
  // run the kernel
  enumerate_kernel<<<1, 1>>>(seq_data, check_idx_data, check_value_data);

  // cuda error check
  DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
  DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

//
// pick
//
__global__ void pick_kernel(
    vecmem::data::vector_view<uint_holder> seq_data,
    vecmem::data::vector_view<dindex> idx_data,
    vecmem::data::vector_view<dindex> check_idx_data,
    vecmem::data::vector_view<dindex> check_value_data) {
  vecmem::device_vector<uint_holder> seq(seq_data);
  vecmem::device_vector<dindex> idx(idx_data);
  vecmem::device_vector<dindex> check_idx(check_idx_data);
  vecmem::device_vector<dindex> check_value(check_value_data);

  for (auto [i, v] : detray::views::pick(seq, idx)) {
    check_idx.push_back(i);
    check_value.push_back(v.ui);
  }
}

void test_pick(vecmem::data::vector_view<uint_holder> seq_data,
               vecmem::data::vector_view<dindex> idx_data,
               vecmem::data::vector_view<dindex> check_idx_data,
               vecmem::data::vector_view<dindex> check_value_data) {
  // run the kernel
  pick_kernel<<<1, 1>>>(seq_data, idx_data, check_idx_data, check_value_data);

  // cuda error check
  DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
  DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

//
// join
//
__global__ void join_kernel(
    vecmem::data::vector_view<uint_holder> seq_data_1,
    vecmem::data::vector_view<uint_holder> seq_data_2,
    vecmem::data::vector_view<dindex> check_value_data) {
  vecmem::device_vector<uint_holder> seq_1(seq_data_1);
  vecmem::device_vector<uint_holder> seq_2(seq_data_2);
  vecmem::device_vector<dindex> check_value(check_value_data);
  std::array<vecmem::device_vector<uint_holder>, 2> vectors{seq_1, seq_2};

  for (auto v : detray::views::join(vectors)) {
    check_value.push_back(v.ui);
  }
}

void test_join(vecmem::data::vector_view<uint_holder> seq_data_1,
               vecmem::data::vector_view<uint_holder> seq_data_2,
               vecmem::data::vector_view<dindex> check_value_data) {
  // run the kernel
  join_kernel<<<1, 1>>>(seq_data_1, seq_data_2, check_value_data);

  // cuda error check
  DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
  DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

//
// static_join
//
__global__ void static_join_kernel(
    vecmem::data::vector_view<uint_holder> seq_data_1,
    vecmem::data::vector_view<uint_holder> seq_data_2,
    vecmem::data::vector_view<dindex> check_value_data) {
  vecmem::device_vector<uint_holder> seq_1(seq_data_1);
  vecmem::device_vector<uint_holder> seq_2(seq_data_2);
  vecmem::device_vector<dindex> check_value(check_value_data);

  for (auto v : detray::views::static_join(seq_1, seq_2)) {
    check_value.push_back(v.ui);
  }
}

void test_static_join(vecmem::data::vector_view<uint_holder> seq_data_1,
                      vecmem::data::vector_view<uint_holder> seq_data_2,
                      vecmem::data::vector_view<dindex> check_value_data) {
  // run the kernel
  static_join_kernel<<<1, 1>>>(seq_data_1, seq_data_2, check_value_data);

  // cuda error check
  DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
  DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

//
// subrange
//
__global__ void subrange_kernel(vecmem::data::vector_view<int> seq_data,
                                vecmem::data::vector_view<int> check_value_data,
                                const std::size_t begin,
                                const std::size_t end) {
  vecmem::device_vector<int> seq(seq_data);
  vecmem::device_vector<int> check(check_value_data);

  for (const auto& v :
       detray::ranges::subrange(seq, std::array<std::size_t, 2>{begin, end})) {
    check.push_back(v);
  }
}

void test_subrange(vecmem::data::vector_view<int> seq_data,
                   vecmem::data::vector_view<int> check_value_data,
                   const std::size_t begin, const std::size_t end) {
  // run the kernel
  subrange_kernel<<<1, 1>>>(seq_data, check_value_data, begin, end);

  // cuda error check
  DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
  DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

}  // namespace detray
