// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s)
#include "detray/core/detail/multi_store.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes.hpp"

// Detray test include(s)
#include "detray/test/framework/types.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// Google Test include(s)
#include <gtest/gtest.h>

// This tests the construction of a static transform store
GTEST_TEST(detray_core, static_mask_store) {
  vecmem::host_memory_resource host_mr;

  using namespace detray;

  using test_algebra = test::algebra;

  enum class mask_id : unsigned int {
    e_rectangle2D = 0u,
    e_trapezoid2D = 1u,
    e_annulus2D = 2u,
    e_cylinder3D = 3u,
    e_ring2D = 4u,
    e_single3D = 5u,
  };

  using rectangle = mask<rectangle2D, test_algebra>;
  using trapezoid = mask<trapezoid2D, test_algebra>;
  using annulus = mask<annulus2D, test_algebra>;
  using cylinder = mask<cylinder2D, test_algebra>;
  using ring = mask<ring2D, test_algebra>;
  using single = mask<single3D<>, test_algebra>;

  // Types must be sorted according to their id (here: masks/mask_identifier)
  using mask_container_t =
      regular_multi_store<mask_id, empty_context, dtuple, dvector, rectangle,
                          trapezoid, annulus, cylinder, ring, single>;

  mask_container_t store(host_mr);

  ASSERT_TRUE(store.empty<mask_id::e_annulus2D>());
  ASSERT_TRUE(store.empty<mask_id::e_cylinder3D>());
  ASSERT_TRUE(store.empty<mask_id::e_rectangle2D>());
  ASSERT_TRUE(store.empty<mask_id::e_ring2D>());
  ASSERT_TRUE(store.empty<mask_id::e_single3D>());
  ASSERT_TRUE(store.empty<mask_id::e_trapezoid2D>());

  store.emplace_back<mask_id::e_cylinder3D>(empty_context{}, 0u, 1.f, 0.5f,
                                            2.0f);

  ASSERT_TRUE(store.empty<mask_id::e_annulus2D>());
  ASSERT_EQ(store.size<mask_id::e_cylinder3D>(), 1);
  ASSERT_TRUE(store.empty<mask_id::e_rectangle2D>());
  ASSERT_TRUE(store.empty<mask_id::e_ring2D>());
  ASSERT_TRUE(store.empty<mask_id::e_single3D>());
  ASSERT_TRUE(store.empty<mask_id::e_trapezoid2D>());

  store.emplace_back<mask_id::e_cylinder3D>(empty_context{}, 0u, 1.f, 1.5f,
                                            2.0f);
  store.emplace_back<mask_id::e_trapezoid2D>(empty_context{}, 0u, 0.5f, 1.5f,
                                             4.0f, 1.f / 8.f);
  store.emplace_back<mask_id::e_rectangle2D>(empty_context{}, 0u, 1.0f, 2.0f);
  store.emplace_back<mask_id::e_rectangle2D>(empty_context{}, 0u, 2.0f, 1.0f);
  store.emplace_back<mask_id::e_rectangle2D>(empty_context{}, 0u, 10.0f,
                                             100.0f);

  ASSERT_TRUE(store.empty<mask_id::e_annulus2D>());
  ASSERT_EQ(store.size<mask_id::e_cylinder3D>(), 2);
  ASSERT_EQ(store.size<mask_id::e_rectangle2D>(), 3);
  ASSERT_TRUE(store.empty<mask_id::e_ring2D>());
  ASSERT_TRUE(store.empty<mask_id::e_single3D>());
  ASSERT_EQ(store.size<mask_id::e_trapezoid2D>(), 1);

  store.emplace_back<mask_id::e_annulus2D>(empty_context{}, 0u, 1.f, 2.f, 3.f,
                                           4.f, 5.f, 6.f, 7.f);
  store.emplace_back<mask_id::e_ring2D>(empty_context{}, 0u, 10.f, 100.f);
  store.emplace_back<mask_id::e_ring2D>(empty_context{}, 0u, 10.f, 100.f);
  store.emplace_back<mask_id::e_ring2D>(empty_context{}, 0u, 10.f, 100.f);
  store.emplace_back<mask_id::e_ring2D>(empty_context{}, 0u, 10.f, 100.f);

  const auto &annulus_masks = store.get<mask_id::e_annulus2D>();
  const auto &cylinder_masks = store.get<mask_id::e_cylinder3D>();
  const auto &rectangle_masks = store.get<mask_id::e_rectangle2D>();
  const auto &ring_masks = store.get<mask_id::e_ring2D>();
  const auto &single_masks = store.get<mask_id::e_single3D>();
  const auto &trapezoid_masks = store.get<mask_id::e_trapezoid2D>();

  ASSERT_TRUE(annulus_masks.size() == 1);
  ASSERT_TRUE(cylinder_masks.size() == 2);
  ASSERT_TRUE(rectangle_masks.size() == 3);
  ASSERT_TRUE(ring_masks.size() == 4);
  ASSERT_TRUE(single_masks.empty());
  ASSERT_TRUE(trapezoid_masks.size() == 1);
}
