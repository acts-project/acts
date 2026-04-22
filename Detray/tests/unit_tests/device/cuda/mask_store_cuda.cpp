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

// Detray test include(s)
#include "mask_store_cuda_kernel.hpp"

// Vecmem include(s)
#include <vecmem/memory/cuda/managed_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <climits>
#include <cstdlib>
#include <random>

using namespace detray;

TEST(mask_store_cuda, mask_store) {
  // memory resource
  vecmem::cuda::managed_memory_resource mng_mr;

  // Types must be sorted according to their id (here: masks/mask_identifier)
  host_store_type store(mng_mr);

  ASSERT_TRUE(store.template empty<mask_id::e_annulus2D>());
  ASSERT_TRUE(store.template empty<mask_id::e_cylinder2D>());
  ASSERT_TRUE(store.template empty<mask_id::e_rectangle2D>());
  ASSERT_TRUE(store.template empty<mask_id::e_ring2D>());
  ASSERT_TRUE(store.template empty<mask_id::e_single3D>());
  ASSERT_TRUE(store.template empty<mask_id::e_trapezoid2D>());

  store.template emplace_back<mask_id::e_rectangle2D>(empty_context{}, 0u, 1.0f,
                                                      2.0f);
  store.template emplace_back<mask_id::e_trapezoid2D>(empty_context{}, 0u, 0.5f,
                                                      1.5f, 4.0f, 1.f / 8.f);
  store.template emplace_back<mask_id::e_ring2D>(empty_context{}, 0u, 1.0f,
                                                 10.0f);
  store.template emplace_back<mask_id::e_cylinder2D>(empty_context{}, 0u, 1.f,
                                                     0.5f, 2.0f);
  store.template emplace_back<mask_id::e_annulus2D>(
      empty_context{}, 0u, 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f);

  ASSERT_FALSE(store.empty<mask_id::e_annulus2D>());
  ASSERT_FALSE(store.empty<mask_id::e_cylinder2D>());
  ASSERT_FALSE(store.empty<mask_id::e_rectangle2D>());
  ASSERT_FALSE(store.empty<mask_id::e_ring2D>());
  ASSERT_TRUE(store.empty<mask_id::e_single3D>());
  ASSERT_FALSE(store.empty<mask_id::e_trapezoid2D>());

  // Generate random points for test
  vecmem::vector<point3> input_point3(n_points, &mng_mr);

  std::random_device rd;
  std::mt19937_64 gen(rd());
  std::uniform_real_distribution<dscalar<test_algebra>> dist(0.f, 9.9f);
  for (unsigned int i = 0u; i < n_points; i++) {
    point3 rand_point3 = {dist(gen), dist(gen), dist(gen)};
    input_point3[i] = rand_point3;
  }

  // host output for intersection status
  vecmem::jagged_vector<int> output_host(5, &mng_mr);

  // get mask objects
  const auto& rectangle_mask = store.get<mask_id::e_rectangle2D>()[0];
  const auto& trapezoid_mask = store.get<mask_id::e_trapezoid2D>()[0];
  const auto& ring_mask = store.get<mask_id::e_ring2D>()[0];
  const auto& cylinder_mask = store.get<mask_id::e_cylinder2D>()[0];
  const auto& annulus_mask = store.get<mask_id::e_annulus2D>()[0];

  // get host results from is_inside function
  for (unsigned int i = 0u; i < n_points; i++) {
    output_host[0].push_back(rectangle_mask.is_inside(input_point3[i]));
    output_host[1].push_back(trapezoid_mask.is_inside(input_point3[i]));
    output_host[2].push_back(ring_mask.is_inside(input_point3[i]));
    output_host[3].push_back(cylinder_mask.is_inside(input_point3[i]));
    output_host[4].push_back(annulus_mask.is_inside(input_point3[i]));
  }

  // Helper object for performing memory copies
  vecmem::cuda::copy copy;

  // device output for intersection status
  vecmem::data::jagged_vector_buffer<int> output_buffer(
      {n_points, n_points, n_points, n_points, n_points}, mng_mr, nullptr,
      vecmem::data::buffer_type::resizable);

  copy.setup(output_buffer)->wait();

  auto input_point3_data = vecmem::get_data(input_point3);
  auto store_data = get_data(store);

  // run the kernel
  mask_test(store_data, input_point3_data, output_buffer);

  vecmem::jagged_vector<int> output_device(&mng_mr);
  copy(output_buffer, output_device)->wait();

  // Compare the values
  for (unsigned int i = 0u; i < n_points; i++) {
    ASSERT_EQ(output_host[0][i], output_device[0][i]);
    ASSERT_EQ(output_host[0][i], output_device[0][i]);
    ASSERT_EQ(output_host[2][i], output_device[2][i]);
    ASSERT_EQ(output_host[3][i], output_device[3][i]);
    ASSERT_EQ(output_host[4][i], output_device[4][i]);
  }
}
