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

// Detray include(s)
#include "container_cuda_kernel.hpp"

// Vecmem include(s)
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/cuda/managed_memory_resource.hpp>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/utils/cuda/copy.hpp>

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <algorithm>
#include <iterator>
#include <numeric>

using namespace detray;

constexpr double tol{1e-6};

/// Test the access to the single store by summing the contained values
TEST(container_cuda, single_store) {
  // Vecmem memory resources
  vecmem::host_memory_resource host_mr;
  vecmem::cuda::managed_memory_resource mng_mr;
  vecmem::cuda::device_memory_resource dev_mr;
  vecmem::cuda::copy cpy;

  // Create single store(s)
  single_store_t store(host_mr);
  single_store_t mng_store(mng_mr);

  // Base store function check
  EXPECT_TRUE(store.empty());
  EXPECT_EQ(store.size(), 0u);
  EXPECT_TRUE(mng_store.empty());
  EXPECT_EQ(mng_store.size(), 0u);

  // Test the managed memory allocation
  geometry_context ctx{};
  mng_store.reserve(4, ctx);
  mng_store.emplace_back(ctx, 1.);
  mng_store.push_back(2., ctx);
  mng_store.insert(vecmem::vector<double>{10.5, 7.6}, ctx);
  store.append(mng_store, ctx);

  EXPECT_FALSE(store.empty());
  EXPECT_EQ(store.size(), 4u);
  EXPECT_FALSE(mng_store.empty());
  EXPECT_EQ(mng_store.size(), 4u);

  // Check the host-side access to the data
  EXPECT_NEAR(mng_store.at(0, ctx), 1., tol);
  EXPECT_NEAR(mng_store.at(2, ctx), 10.5, tol);
  EXPECT_NEAR(mng_store.at(1, ctx), 2., tol);
  EXPECT_NEAR(mng_store.at(3, ctx), 7.6, tol);

  // CPU sum check
  double cpu_sum{std::accumulate(mng_store.begin(), mng_store.end(), 0.)};
  EXPECT_NEAR(cpu_sum, 21.1, tol);
  EXPECT_NEAR(cpu_sum, std::accumulate(store.begin(), store.end(), 0.), tol);

  // Get the view to the managed memory and run the test
  single_store_t::view_type mng_store_view = detray::get_data(mng_store);

  vecmem::vector<double> cuda_sum(&mng_mr);
  cuda_sum.push_back(0.);
  vecmem::data::vector_view<double> sum_data = vecmem::get_data(cuda_sum);

  test_single_store(mng_store_view, sum_data);

  EXPECT_NEAR(cpu_sum, cuda_sum[0], tol);

  // Check that the store can be cleared
  mng_store.clear(ctx);

  EXPECT_TRUE(mng_store.empty());
  EXPECT_EQ(mng_store.size(), 0u);

  // Reset for next test
  cuda_sum[0] = 0.;

  // Copy the host store to device, get the view to it and run the test again
  single_store_t::buffer_type store_buffer =
      detray::get_buffer(store, dev_mr, cpy);
  single_store_t::view_type buffer_view = detray::get_data(store_buffer);

  EXPECT_NEAR(cuda_sum[0], 0.f, tol);

  test_single_store(buffer_view, sum_data);

  EXPECT_NEAR(cuda_sum[0], cpu_sum, tol);
}

/// Test the access to the tuple container by summing the contained values
TEST(container_cuda, tuple_container) {
  // Vecmem memory resources
  vecmem::host_memory_resource host_mr;
  vecmem::cuda::managed_memory_resource mng_mr;
  vecmem::cuda::device_memory_resource dev_mr;
  vecmem::cuda::copy cpy;

  // Create single store(s)
  tuple_cont_t tcont(host_mr);
  tuple_cont_t mng_tcont(mng_mr);

  // Base store function check
  EXPECT_EQ(tcont.size(), 2u);
  EXPECT_EQ(mng_tcont.size(), 2u);

  // Test the managed memory allocation
  mng_tcont.get<0>().push_back(3);
  mng_tcont.get<1>().push_back(4.);
  mng_tcont.get<1>().push_back(5.8);

  vecmem::vector<double> new_data{10.5, 7.6, 14.5};
  auto& data_coll = detail::get<1>(mng_tcont);
  std::ranges::copy(new_data, std::back_inserter(data_coll));
  std::ranges::copy(mng_tcont.get<0>(), std::back_inserter(tcont.get<0>()));
  std::ranges::copy(mng_tcont.get<1>(), std::back_inserter(tcont.get<1>()));

  EXPECT_EQ(mng_tcont.get<0>().size(), 1u);
  EXPECT_EQ(mng_tcont.get<1>().size(), 5u);
  EXPECT_EQ(tcont.get<0>().size(), 1u);
  EXPECT_EQ(tcont.get<1>().size(), 5u);

  // CPU sum check
  double cpu_sum{static_cast<double>(mng_tcont.get<0>().at(0))};
  cpu_sum = std::accumulate(data_coll.begin(), data_coll.end(), cpu_sum);
  EXPECT_NEAR(cpu_sum, 45.4, tol);

  // Get the view to the managed memory and run the test
  tuple_cont_t::view_type mng_store_view = detray::get_data(mng_tcont);

  vecmem::vector<double> cuda_sum(&mng_mr);
  cuda_sum.push_back(0.);
  vecmem::data::vector_view<double> sum_data = vecmem::get_data(cuda_sum);

  test_tuple_container(mng_store_view, sum_data);

  EXPECT_NEAR(cpu_sum, cuda_sum[0], tol);

  // Reset for next test
  cuda_sum[0] = 0.;

  // Copy the host store to device, get the view to it and run the test again
  tuple_cont_t::buffer_type store_buffer =
      detray::get_buffer(tcont, dev_mr, cpy);
  tuple_cont_t::view_type buffer_view = detray::get_data(store_buffer);

  EXPECT_NEAR(cuda_sum[0], 0.f, tol);

  test_tuple_container(buffer_view, sum_data);

  EXPECT_NEAR(cuda_sum[0], cpu_sum, tol);
}

/// Test the regular multi store functionality
TEST(container_cuda, regular_multi_store) {
  using enum reg_type_ids;

  // Vecmem memory resources
  vecmem::host_memory_resource host_mr;
  vecmem::cuda::managed_memory_resource mng_mr;
  vecmem::cuda::device_memory_resource dev_mr;
  vecmem::cuda::copy cpy;

  // Create tuple vector store
  reg_multi_store_t mng_store(mng_mr);
  reg_multi_store_t store(host_mr);

  // Base store function check
  EXPECT_EQ(mng_store.n_collections(), 3u);
  EXPECT_EQ(mng_store.empty<e_size>(), true);
  EXPECT_EQ(mng_store.empty<e_float>(), true);
  EXPECT_EQ(mng_store.empty<e_double>(), true);

  // Add elements to the store
  empty_context ctx{};
  mng_store.emplace_back<e_size>(ctx, 1u);
  mng_store.push_back<e_size>(2u, ctx);
  mng_store.emplace_back<e_float>(ctx, 3.1f);
  mng_store.push_back<e_float>(4.5f, ctx);
  mng_store.emplace_back<e_double>(ctx, 5.5);
  mng_store.push_back<e_double>(6.0, ctx);

  vecmem::vector<std::size_t> int_vec{3u, 4u, 5u};
  mng_store.insert(int_vec);

  vecmem::vector<float> float_vec{12.1f, 5.6f};
  mng_store.insert(float_vec);

  mng_store.insert(vecmem::vector<double>{10.5, 7.6});

  store.append(mng_store, ctx);

  EXPECT_EQ(mng_store.size<e_size>(), 5u);
  EXPECT_EQ(mng_store.size<e_float>(), 4u);
  EXPECT_EQ(mng_store.size<e_double>(), 4u);
  EXPECT_EQ(store.size<e_size>(), 5u);
  EXPECT_EQ(store.size<e_float>(), 4u);
  EXPECT_EQ(store.size<e_double>(), 4u);

  // CPU sum check
  double cpu_sum{std::accumulate(mng_store.get<e_size>().begin(),
                                 mng_store.get<e_size>().end(), 0.)};
  cpu_sum = std::accumulate(mng_store.get<e_float>().begin(),
                            mng_store.get<e_float>().end(), cpu_sum);
  cpu_sum = std::accumulate(mng_store.get<e_double>().begin(),
                            mng_store.get<e_double>().end(), cpu_sum);
  EXPECT_NEAR(cpu_sum, 69.9, tol);

  // CUDA sum check
  typename reg_multi_store_t::view_type store_view = get_data(mng_store);

  vecmem::vector<double> cuda_sum(&mng_mr);
  cuda_sum.push_back(0.);
  auto sum_data = vecmem::get_data(cuda_sum);

  test_reg_multi_store(store_view, sum_data);

  EXPECT_NEAR(cpu_sum, cuda_sum[0], tol);

  // Reset for next test
  cuda_sum[0] = 0.;

  // Copy the host store to device, get the view to it and run the test again
  reg_multi_store_t::buffer_type store_buffer =
      detray::get_buffer(store, dev_mr, cpy);
  reg_multi_store_t::view_type buffer_view = detray::get_data(store_buffer);

  EXPECT_NEAR(cuda_sum[0], 0.f, tol);

  test_reg_multi_store(buffer_view, sum_data);

  EXPECT_NEAR(cuda_sum[0], cpu_sum, tol);
}

/// Tets the multi store with a hierarchical memory type ( @c test<> )
TEST(container_cuda, multi_store) {
  using enum type_ids;

  // Vecmem memory resources
  vecmem::host_memory_resource host_mr;
  vecmem::cuda::managed_memory_resource mng_mr;
  vecmem::cuda::device_memory_resource dev_mr;
  vecmem::cuda::copy cpy;

  // Create tuple vector store
  multi_store_t mng_store(mng_mr);
  multi_store_t store(host_mr);

  // Base store function check
  EXPECT_EQ(mng_store.n_collections(), 2u);
  EXPECT_EQ(mng_store.empty<e_float>(), true);

  // Add elements to the store
  vecmem::vector<float> float_vec{12.1f, 5.6f};
  mng_store.insert(float_vec);

  mng_store.get<e_test_class>().first = vecmem::vector<int>{2, 3};
  mng_store.get<e_test_class>().second = vecmem::vector<double>{18., 42.6};

  std::ranges::copy(mng_store.get<e_float>(),
                    std::back_inserter(store.get<e_float>()));
  store.get<e_test_class>() = mng_store.get<e_test_class>();

  EXPECT_EQ(mng_store.size<e_float>(), 2u);
  EXPECT_EQ(mng_store.get<e_test_class>().first.size(), 2u);
  EXPECT_EQ(mng_store.get<e_test_class>().second.size(), 2u);
  EXPECT_EQ(store.size<e_float>(), 2u);
  EXPECT_EQ(store.get<e_test_class>().first.size(), 2u);
  EXPECT_EQ(store.get<e_test_class>().second.size(), 2u);

  // CPU sum check
  double cpu_sum{std::accumulate(mng_store.get<e_float>().begin(),
                                 mng_store.get<e_float>().end(), 0.)};
  cpu_sum = std::accumulate(mng_store.get<e_test_class>().first.begin(),
                            mng_store.get<e_test_class>().first.end(), cpu_sum);
  cpu_sum =
      std::accumulate(mng_store.get<e_test_class>().second.begin(),
                      mng_store.get<e_test_class>().second.end(), cpu_sum);
  EXPECT_NEAR(cpu_sum, 83.3, tol);

  // CUDA sum check
  typename multi_store_t::view_type store_view = get_data(mng_store);

  vecmem::vector<double> cuda_sum(&mng_mr);
  cuda_sum.push_back(0.);
  auto sum_data = vecmem::get_data(cuda_sum);

  test_multi_store(store_view, sum_data);

  EXPECT_NEAR(cpu_sum, cuda_sum[0], tol);

  // Reset for next test
  cuda_sum[0] = 0.;

  // Copy the host store to device, get the view to it and run the test again
  multi_store_t::buffer_type store_buffer =
      detray::get_buffer(store, dev_mr, cpy);
  multi_store_t::view_type buffer_view = detray::get_data(store_buffer);

  EXPECT_NEAR(cuda_sum[0], 0.f, tol);

  test_multi_store(buffer_view, sum_data);

  EXPECT_NEAR(cuda_sum[0], cpu_sum, tol);
}
