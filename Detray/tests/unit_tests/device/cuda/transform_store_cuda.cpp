// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/utils/ranges.hpp"

// Detray test include(s)
#include "transform_store_cuda_kernel.hpp"

// Vecmem include(s)
#include <vecmem/memory/cuda/managed_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <climits>

using namespace detray;

TEST(transform_store_cuda, transform_store) {
  // memory resource
  vecmem::cuda::managed_memory_resource mng_mr;

  host_transform_store_t static_store(mng_mr);
  typename host_transform_store_t::context_type ctx0{};
  typename host_transform_store_t::context_type ctx1{};

  ASSERT_TRUE(static_store.empty(ctx0));

  ASSERT_EQ(static_store.size(ctx0), 0u);

  point3 t0{1.f, 2.f, 3.f};
  point3 z0{4.f, 2.f, 3.f};
  point3 x0{1.f, 0.f, 3.f};
  transform3 tf0{t0, z0, x0};
  static_store.push_back(tf0, ctx0);
  ASSERT_EQ(static_store.size(ctx0), 1u);

  point3 t1{1.f, 1.f, 2.f};
  point3 z1{1.f, 0.f, 0.f};
  point3 x1{0.f, 0.f, 2.f};
  transform3 tf1{t1, z1, x1};
  static_store.push_back(tf1, ctx1);
  ASSERT_EQ(static_store.size(ctx1), 2u);

  point3 t2{2.f, 2.f, 5.f};
  point3 z2{0.f, 0.f, 0.f};
  point3 x2{1.f, 2.f, 3.f};
  transform3 tf2{t2, z2, x2};
  static_store.push_back(std::move(tf2), ctx0);
  ASSERT_EQ(static_store.size(ctx0), 3u);

  point3 t3{2.f, 0.f, 5.f};
  point3 z3{1.f, 2.f, 3.f};
  point3 x3{0.f, 0.f, 0.f};
  transform3 tf3{t3, z3, x3};
  static_store.emplace_back(ctx0, std::move(t3));
  ASSERT_EQ(static_store.size(ctx0), 4u);

  static_store.emplace_back(ctx0);
  ASSERT_EQ(static_store.size(ctx0), 5u);

  // Compare the transform operation results
  dindex n_transforms = static_store.size(ctx0);

  // Fill input vector
  vecmem::vector<point3> input(&mng_mr);
  for (unsigned int i = 0u; i < n_transforms; i++) {
    input.push_back({static_cast<scalar>(i), static_cast<scalar>(i + 1u),
                     static_cast<scalar>(i + 2u)});
  }

  // Get transformed vector from host side
  std::vector<point3> output_host;

  auto range = detray::ranges::subrange(static_store.get(ctx0),
                                        dindex_range{0u, n_transforms});
  std::size_t count{0u};
  for (const auto& tf : range) {
    auto output = tf.point_to_global(input[count]);
    output_host.push_back(output);
    count++;
  }

  // Get transformed vector from device side
  vecmem::vector<point3> output_device(n_transforms, &mng_mr);

  auto input_data = vecmem::get_data(input);
  auto output_data = vecmem::get_data(output_device);
  auto static_store_data = get_data(static_store);
  transform_test(input_data, static_store_data, output_data,
                 static_store.size());

  // Compare the values
  for (unsigned int i = 0u; i < n_transforms; i++) {
    ASSERT_EQ(output_host[i], output_device[i]);
  }
}
