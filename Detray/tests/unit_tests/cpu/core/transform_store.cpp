// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s)
#include "detray/core/detail/single_store.hpp"

// Detray test include(s)
#include "detray/test/framework/assert.hpp"
#include "detray/test/framework/types.hpp"

// GTest include(s)
#include <gtest/gtest.h>

// This tests the construction of a static transform store
GTEST_TEST(detray_core, static_transform_store) {
  using namespace detray;
  using transform3 = test::transform3;
  using point3 = test::point3;

  using transform_store_t = single_store<transform3, dvector, geometry_context>;
  transform_store_t static_store;
  typename transform_store_t::context_type ctx0{};
  typename transform_store_t::context_type ctx1{};

  ASSERT_TRUE(static_store.empty(ctx0));

  ASSERT_EQ(static_store.size(ctx0), 0u);

  point3 t0{0.f, 0.f, 0.f};
  transform3 tf0{t0};
  static_store.push_back(tf0, ctx0);
  ASSERT_EQ(static_store.size(ctx0), 1u);

  point3 t1{1.f, 0.f, 0.f};
  transform3 tf1{t1};
  static_store.push_back(tf1, ctx1);
  ASSERT_EQ(static_store.size(ctx1), 2u);

  point3 t2{2.f, 0.f, 0.f};
  transform3 tf2{t2};
  static_store.push_back(std::move(tf2), ctx0);
  ASSERT_EQ(static_store.size(ctx0), 3u);

  point3 t3{2.f, 0.f, 0.f};
  static_store.emplace_back(ctx0, std::move(t3));
  ASSERT_EQ(static_store.size(ctx0), 4u);

  static_store.emplace_back(ctx0);
  ASSERT_EQ(static_store.size(ctx0), 5u);
}

// This tests the construction of a static transform store
GTEST_TEST(detray_core, multicontext_transform_store) {
  using namespace detray;
  using transform3 = test::transform3;
  using point3 = test::point3;

  using transform_store_t = single_store<transform3, dvector, geometry_context>;
  transform_store_t xf_store;
  typename transform_store_t::context_type ctx0{0};
  typename transform_store_t::context_type ctx1{1};

  // Populate the default context (context 0) of the transform store
  // Add new elements one at a time
  ASSERT_TRUE(xf_store.empty(ctx0));

  point3 t0{0.f, 0.f, 0.f};
  transform3 tf0{t0};
  xf_store.push_back(tf0);

  point3 t1{1.f, 10.f, 7.f};
  transform3 tf1{t1};
  xf_store.push_back(tf1);

  point3 t2{2.f, 1.f, 4.f};
  transform3 tf2{t2};
  xf_store.push_back(tf2, ctx0);

  point3 t3{3.f, 8.f, 6.f};
  transform3 tf3{t3};
  xf_store.push_back(tf3);

  ASSERT_EQ(xf_store.size(ctx0), 4u);

  // Create a vector of 'aligned' transforms
  point3 shift{.11f, .12f, .13f};
  typename transform_store_t::base_type aligned_transforms;
  aligned_transforms.reserve(xf_store.size(ctx0));
  point3 t0_shifted = t0 + shift;
  point3 t1_shifted = t1 + shift;
  point3 t2_shifted = t2 + shift;
  point3 t3_shifted = t3 + shift;
  aligned_transforms.push_back(transform3(t0_shifted));
  aligned_transforms.push_back(transform3(t1_shifted));
  aligned_transforms.push_back(transform3(t2_shifted));
  aligned_transforms.push_back(transform3(t3_shifted));

  // Add context 1 to the transform store
  ASSERT_TRUE(xf_store.empty(ctx1));
  xf_store.add_context(aligned_transforms);
  ASSERT_EQ(xf_store.size(ctx1), 4u);

  // Check that we can access context data as expected
  point3 tr0_0 = xf_store.at(0, ctx0).translation();
  point3 tr1_0 = xf_store.at(0, ctx1).translation();
  point3 tr_diff_0 = tr1_0 - tr0_0;
  EXPECT_POINT3_NEAR(tr_diff_0, shift, 1e-6);

  point3 tr0_3 = xf_store.at(3, ctx0).translation();
  point3 tr1_3 = xf_store.at(3, ctx1).translation();
  point3 tr_diff_3 = tr1_3 - tr0_3;
  EXPECT_POINT3_NEAR(tr_diff_3, shift, 1e-6);
}
