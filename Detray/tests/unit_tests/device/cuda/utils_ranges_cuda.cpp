// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Detray test include(s)
#include "utils_ranges_cuda_kernel.hpp"

// Vecmem include(s)
#include <vecmem/memory/cuda/managed_memory_resource.hpp>
#include <vecmem/utils/cuda/copy.hpp>

// GTest include(s)
#include <gtest/gtest.h>

using namespace detray;

// This tests the single value view
TEST(utils_ranges_cuda, single) {
  dindex value{251u};
  dindex check{std::numeric_limits<dindex>::max()};

  // Run the test code
  test_single(value, check);

  // Check result value
  ASSERT_EQ(value, check);
}

// This tests the non-owning single value view
TEST(utils_ranges_cuda, pointer) {
  dindex value{251u};
  dindex check{std::numeric_limits<dindex>::max()};

  // Run the test code
  test_pointer(value, check);

  // Check result value
  ASSERT_EQ(value, check);
}

// This tests the iota range generator
TEST(utils_ranges_cuda, iota) {
  // Helper object for performing memory copies.
  vecmem::cuda::copy copy;

  // memory resource
  vecmem::cuda::managed_memory_resource managed_resource;

  // Input reference vector for test
  vecmem::vector<dindex> reference({2u, 3u, 4u, 5u, 6u});

  // Const range for enumeration test
  const darray<dindex, 2> range = {2u, 7u};

  // Output vector buffer for enumeration test
  vecmem::data::vector_buffer<dindex> check_buffer(
      static_cast<vecmem::data::vector_buffer<dindex>::size_type>(range[1] -
                                                                  range[0]),
      managed_resource, vecmem::data::buffer_type::resizable);
  copy.setup(check_buffer)->wait();

  // Run test function
  test_iota(range, check_buffer);

  // Copy vector buffer to output vector
  vecmem::vector<dindex> check{&managed_resource};
  copy(check_buffer, check)->wait();

  // Check the result
  ASSERT_EQ(check, reference);
}

// This tests the cartesian product range adaptor
TEST(utils_ranges_cuda, cartesian_product) {
  // Helper object for performing memory copies.
  vecmem::cuda::copy copy;

  // memory resource
  vecmem::cuda::managed_memory_resource managed_resource;

  // Const range for enumeration test
  const darray<dindex, 2> range1 = {2u, 7u};
  const darray<dindex, 2> range2 = {1u, 10u};
  const darray<dindex, 2> range3 = {3u, 4u};

  auto seq1 = detray::views::iota(range1);
  auto seq2 = detray::views::iota(range2);
  auto seq3 = detray::views::iota(range3);

  dindex size{seq1.size() * seq2.size() * seq3.size()};

  // Input reference vector for test
  vecmem::vector<std::tuple<dindex, dindex, dindex>> result;
  for (auto i : seq1) {
    for (auto j : seq2) {
      for (auto k : seq3) {
        result.emplace_back(i, j, k);
      }
    }
  }

  // Output vector buffer for enumeration test
  using buffer_t =
      vecmem::data::vector_buffer<std::tuple<dindex, dindex, dindex>>;
  buffer_t check_buffer(static_cast<buffer_t::size_type>(size),
                        managed_resource, vecmem::data::buffer_type::resizable);
  copy.setup(check_buffer)->wait();

  // Run test function
  test_cartesian_product(range1, range2, range3, check_buffer);

  // Copy vector buffer to output vector
  vecmem::vector<std::tuple<dindex, dindex, dindex>> check{&managed_resource};
  copy(check_buffer, check)->wait();

  // Check the result
  ASSERT_EQ(result.size(), check.size());

  for (std::size_t r = 0; r < check.size(); ++r) {
    const auto [i, j, k] = check[r];
    const auto [l, m, n] = result[r];

    ASSERT_EQ(i, l);
    ASSERT_EQ(j, m);
    ASSERT_EQ(k, n);
  }
}

// This tests the convenience enumeration function
TEST(utils_ranges_cuda, enumerate) {
  // Helper object for performing memory copies.
  vecmem::cuda::copy copy;

  // memory resource
  vecmem::cuda::managed_memory_resource managed_resource;

  // Input vector sequence for test
  vecmem::vector<uint_holder> seq({{0u}, {1u}, {2u}, {3u}, {4u}, {5u}},
                                  &managed_resource);

  // Get vector_data object
  auto seq_data = vecmem::get_data(seq);

  // Output vector buffer for enumeration test
  vecmem::data::vector_buffer<dindex> idx_buffer(
      static_cast<vecmem::data::vector_buffer<dindex>::size_type>(seq.size()),
      managed_resource, vecmem::data::buffer_type::resizable);
  copy.setup(idx_buffer)->wait();

  vecmem::data::vector_buffer<dindex> value_buffer(
      static_cast<vecmem::data::vector_buffer<dindex>::size_type>(seq.size()),
      managed_resource, vecmem::data::buffer_type::resizable);
  copy.setup(value_buffer)->wait();

  // Run test function
  test_enumerate(seq_data, idx_buffer, value_buffer);

  // Copy vector buffer to output vector
  vecmem::vector<dindex> idx_vec{&managed_resource};
  copy(idx_buffer, idx_vec)->wait();

  vecmem::vector<dindex> value_vec{&managed_resource};
  copy(value_buffer, value_vec)->wait();

  // Check the result
  for (std::size_t i = 0u; i < idx_vec.size(); i++) {
    ASSERT_EQ(idx_vec[i], value_vec[i]);
  }
}

// This tests the convenience pick function
TEST(utils_ranges_cuda, pick) {
  // Helper object for performing memory copies.
  vecmem::cuda::copy copy;

  // memory resource
  vecmem::cuda::managed_memory_resource managed_resource;

  // Input vector sequence for test
  vecmem::vector<uint_holder> seq({{0u}, {1u}, {2u}, {3u}, {4u}, {5u}},
                                  &managed_resource);
  // Input index sequence for test
  vecmem::vector<dindex> idx({0u, 2u, 4u, 5u}, &managed_resource);

  // Get vector_data object
  auto seq_data = vecmem::get_data(seq);
  auto idx_data = vecmem::get_data(idx);

  // Output vector buffer for enumeration test
  vecmem::data::vector_buffer<dindex> idx_buffer(
      static_cast<vecmem::data::vector_buffer<dindex>::size_type>(seq.size()),
      managed_resource, vecmem::data::buffer_type::resizable);
  copy.setup(idx_buffer)->wait();

  vecmem::data::vector_buffer<dindex> value_buffer(
      static_cast<vecmem::data::vector_buffer<dindex>::size_type>(seq.size()),
      managed_resource, vecmem::data::buffer_type::resizable);
  copy.setup(value_buffer)->wait();

  // Run test function
  test_pick(seq_data, idx_data, idx_buffer, value_buffer);

  // Copy vector buffer to output vector
  vecmem::vector<dindex> idx_vec{&managed_resource};
  copy(idx_buffer, idx_vec)->wait();

  vecmem::vector<dindex> value_vec{&managed_resource};
  copy(value_buffer, value_vec)->wait();

  // Check the result
  for (std::size_t i = 0u; i < idx_vec.size(); i++) {
    ASSERT_EQ(idx[i], idx_vec[i]);
    ASSERT_EQ(seq[idx_vec[i]].ui, value_vec[i]);
  }
}

// This tests the detray join range
TEST(utils_ranges_cuda, join) {
  // Helper object for performing memory copies.
  vecmem::cuda::copy copy;

  // memory resource
  vecmem::cuda::managed_memory_resource managed_resource;

  // Input vector sequence for test
  vecmem::vector<uint_holder> seq_1({{0u}, {1u}, {2u}, {3u}, {4u}, {5u}},
                                    &managed_resource);
  vecmem::vector<uint_holder> seq_2({{2u}, {0u}, {9u}, {4u}, {15u}},
                                    &managed_resource);
  // Get vector_data object
  auto seq_data_1 = vecmem::get_data(seq_1);
  auto seq_data_2 = vecmem::get_data(seq_2);

  // Output vector buffer for join test
  vecmem::data::vector_buffer<dindex> value_buffer(
      static_cast<vecmem::data::vector_buffer<dindex>::size_type>(seq_1.size() +
                                                                  seq_2.size()),
      managed_resource, vecmem::data::buffer_type::resizable);
  copy.setup(value_buffer)->wait();

  // Run test function
  test_join(seq_data_1, seq_data_2, value_buffer);

  // Copy vector buffer to output vector
  vecmem::vector<dindex> value_vec{&managed_resource};
  copy(value_buffer, value_vec)->wait();

  // First sequence
  for (std::size_t i = 0u; i < seq_1.size(); i++) {
    ASSERT_EQ(seq_1[i].ui, value_vec[i]);
  }
  // Second sequence
  for (std::size_t i = 0u; i < seq_2.size(); i++) {
    ASSERT_EQ(seq_2[i].ui, value_vec[i + seq_1.size()]);
  }
}

// This tests the convenience enumeration function
TEST(utils_ranges_cuda, static_join) {
  // Helper object for performing memory copies.
  vecmem::cuda::copy copy;

  // memory resource
  vecmem::cuda::managed_memory_resource managed_resource;

  // Input vector sequence for test
  vecmem::vector<uint_holder> seq_1({{0u}, {1u}, {2u}, {3u}, {4u}, {5u}},
                                    &managed_resource);
  vecmem::vector<uint_holder> seq_2({{2u}, {0u}, {9u}, {4u}, {15u}},
                                    &managed_resource);
  // Get vector_data object
  auto seq_data_1 = vecmem::get_data(seq_1);
  auto seq_data_2 = vecmem::get_data(seq_2);

  // Output vector buffer for static_join test
  vecmem::data::vector_buffer<dindex> value_buffer(
      static_cast<vecmem::data::vector_buffer<dindex>::size_type>(seq_1.size() +
                                                                  seq_2.size()),
      managed_resource, vecmem::data::buffer_type::resizable);
  copy.setup(value_buffer)->wait();

  // Run test function
  test_static_join(seq_data_1, seq_data_2, value_buffer);

  // Copy vector buffer to output vector
  vecmem::vector<dindex> value_vec{&managed_resource};
  copy(value_buffer, value_vec)->wait();

  // First sequence
  for (std::size_t i = 0u; i < seq_1.size(); i++) {
    ASSERT_EQ(seq_1[i].ui, value_vec[i]);
  }
  // Second sequence
  for (std::size_t i = 0u; i < seq_1.size(); i++) {
    ASSERT_EQ(seq_2[i].ui, value_vec[i + seq_1.size()]);
  }
}

// This tests the subrange view
TEST(utils_ranges_cuda, subrange) {
  // Helper object for performing memory copies.
  vecmem::cuda::copy copy;

  // memory resource
  vecmem::cuda::managed_memory_resource managed_resource;

  // Input vector sequence for test
  vecmem::vector<int> seq({0, 1, 2, 3, 4, 5}, &managed_resource);
  // Get vector_data object
  auto seq_data = vecmem::get_data(seq);

  // Begin and end index for iteration
  const std::size_t begin{1u};
  const std::size_t end{4u};

  // Output vector buffer for iteration test
  vecmem::data::vector_buffer<int> check_buffer(
      static_cast<vecmem::data::vector_buffer<int>::size_type>(end - begin),
      managed_resource, vecmem::data::buffer_type::resizable);
  copy.setup(check_buffer)->wait();

  // Run test function
  test_subrange(seq_data, check_buffer, begin, end);

  // Copy vector buffer to output vector
  vecmem::vector<int> check{&managed_resource};
  copy(check_buffer, check)->wait();

  // Check the result
  ASSERT_EQ(check[0], 1);
  ASSERT_EQ(check[1], 2);
  ASSERT_EQ(check[2], 3);
}
