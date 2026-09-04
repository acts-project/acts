// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s)
#include "detray/definitions/units.hpp"
#include "detray/geometry/volume_descriptor.hpp"

// Detray test include(s)
#include "detray/test/common/build_toy_detector.hpp"
#include "detray/test/framework/types.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <algorithm>

// TODO: Move these into the test defs
namespace {

// geo object ids for testing
enum geo_objects : unsigned int {
  e_sensitive = 0u,
  e_portal = 1u,
  e_size = 2u,
  e_all = e_size,
};

// surface finder ids for testing
enum class accel_id : unsigned int {
  e_surface_default = 0u,
  e_grid = 1u,
};

}  // namespace

// This tests the detector volume class and its many links
GTEST_TEST(detray_geometry, volume_descriptor) {
  using namespace detray;

  using accel_link_t = dtyped_index<accel_id, dindex>;
  using volume_t = volume_descriptor<geo_objects, accel_link_t>;

  // Check construction, setters and getters
  volume_t v1(volume_id::e_cylinder);
  v1.set_index(12345u);
  v1.template set_accel_link<geo_objects::e_portal>(
      {accel_id::e_surface_default, 1u});
  v1.template set_accel_link<geo_objects::e_sensitive>({accel_id::e_grid, 12u});

  ASSERT_TRUE(v1.id() == volume_id::e_cylinder);
  ASSERT_TRUE(v1.index() == 12345u);
  ASSERT_TRUE(v1.template accel_link<geo_objects::e_portal>().id() ==
              accel_id::e_surface_default);
  ASSERT_TRUE(v1.template accel_link<geo_objects::e_portal>().index() == 1u);
  ASSERT_TRUE(v1.template accel_link<geo_objects::e_sensitive>().id() ==
              accel_id::e_grid);
  ASSERT_TRUE(v1.template accel_link<geo_objects::e_sensitive>().index() ==
              12u);

  // Check copy constructor
  const auto v2 = volume_t(v1);
  ASSERT_EQ(v2.id(), volume_id::e_cylinder);
  ASSERT_EQ(v2.index(), 12345u);
  ASSERT_TRUE(v2.template accel_link<geo_objects::e_portal>().id() ==
              accel_id::e_surface_default);
  ASSERT_TRUE(v2.template accel_link<geo_objects::e_portal>().index() == 1u);
  ASSERT_TRUE(v2.template accel_link<geo_objects::e_sensitive>().id() ==
              accel_id::e_grid);
  ASSERT_TRUE(v2.template accel_link<geo_objects::e_sensitive>().index() ==
              12u);
}

/// This tests the functionality of a detector volume interface
GTEST_TEST(detray_geometry, tracking_volume) {
  using namespace detray;

  using test_algebra = test::algebra;
  using scalar = test::scalar;

  constexpr scalar tol{5e-5f};

  vecmem::host_memory_resource host_mr;
  const auto [toy_det, names] = build_toy_detector<test_algebra>(host_mr);

  //
  // Volume 7 is a barrrel layer with sensitive surfaces
  //
  const auto vol7 = tracking_volume{toy_det, 7u};

  ASSERT_EQ(vol7.id(), volume_id::e_cylinder) << vol7 << std::endl;
  ASSERT_EQ(vol7.index(), 7u) << vol7 << std::endl;
  ASSERT_EQ(vol7.name(names), "barrel_7") << vol7 << std::endl;
  auto t = vol7.center();
  ASSERT_NEAR(t[0], 0.f, tol);
  ASSERT_NEAR(t[1], 0.f, tol);
  ASSERT_NEAR(t[2], 0.f, tol);
  ASSERT_EQ(vol7.surfaces().size(), 228u);

  // Access to all surfaces
  std::vector<dindex> sf_indices{};
  sf_indices.reserve(vol7.surfaces().size());
  for (const auto& sf : vol7.surfaces()) {
    sf_indices.push_back(sf.index());
  }
  auto seq = detray::views::iota(354u, 582u);
  EXPECT_TRUE(std::ranges::equal(sf_indices, seq));

  // Access to portals
  sf_indices.clear();
  for (const auto& pt : vol7.portals()) {
    sf_indices.push_back(pt.index());
  }
  seq = detray::views::iota(578u, 582u);
  EXPECT_TRUE(std::ranges::equal(sf_indices, seq));

  // Access to sensitive surfaces
  sf_indices.clear();
  for (const auto& sens : vol7.template surfaces<surface_id::e_sensitive>()) {
    sf_indices.push_back(sens.index());
  }
  seq = detray::views::iota(354u, 578u);
  EXPECT_TRUE(std::ranges::equal(sf_indices, seq));
  //
  // Volume 5 is negative endcap layer with sensitive surfaces
  //
  const auto vol5 = tracking_volume{toy_det, 5u};

  ASSERT_EQ(vol5.id(), volume_id::e_cylinder) << vol5 << std::endl;
  ASSERT_EQ(vol5.index(), 5u) << vol5 << std::endl;
  ASSERT_EQ(vol5.name(names), "endcap_5") << vol5 << std::endl;
  t = vol5.center();
  ASSERT_NEAR(t[0], 0.f, tol);
  ASSERT_NEAR(t[1], 0.f, tol);
  ASSERT_NEAR(t[2], -820.f, tol);
  ASSERT_EQ(vol5.surfaces().size(), 112u);

  // Access to all surfaces
  sf_indices.clear();
  for (const auto& sf : vol5.surfaces()) {
    sf_indices.push_back(sf.index());
  }
  seq = detray::views::iota(238u, 350u);
  EXPECT_TRUE(std::ranges::equal(sf_indices, seq));

  // Access to portals
  sf_indices.clear();
  for (const auto& pt : vol5.portals()) {
    sf_indices.push_back(pt.index());
  }
  seq = detray::views::iota(346u, 350u);
  EXPECT_TRUE(std::ranges::equal(sf_indices, seq));

  // Access to sensitive surfaces
  sf_indices.clear();
  for (const auto& sens : vol5.template surfaces<surface_id::e_sensitive>()) {
    sf_indices.push_back(sens.index());
  }
  seq = detray::views::iota(238u, 346u);
  EXPECT_TRUE(std::ranges::equal(sf_indices, seq));

  //
  // Volume 17 is the positive connector gap
  //
  const auto vol17 = tracking_volume{toy_det, 17u};

  ASSERT_EQ(vol17.id(), volume_id::e_cylinder) << vol17 << std::endl;
  ASSERT_EQ(vol17.index(), 17u) << vol17 << std::endl;
  ASSERT_EQ(vol17.name(names), "connector_gap_17") << vol17 << std::endl;
  t = vol17.center();
  ASSERT_NEAR(t[0], 0.f, tol);
  ASSERT_NEAR(t[1], 0.f, tol);
  ASSERT_NEAR(t[2], 546.25f, tol);
  ASSERT_EQ(vol17.surfaces().size(), 4u);

  // Access to all surfaces
  sf_indices.clear();
  for (const auto& sf : vol17.surfaces()) {
    sf_indices.push_back(sf.index());
  }
  seq = detray::views::iota(2994u, 2998u);
  EXPECT_TRUE(std::ranges::equal(sf_indices, seq));

  // Access to portals
  sf_indices.clear();
  for (const auto& pt : vol17.portals()) {
    sf_indices.push_back(pt.index());
  }
  EXPECT_TRUE(std::ranges::equal(sf_indices, seq));

  // Access to sensitive surfaces: None in gap volume
  sf_indices.clear();
  for (const auto& sens : vol17.template surfaces<surface_id::e_sensitive>()) {
    sf_indices.push_back(sens.index());
  }
  EXPECT_TRUE(sf_indices.empty());
}
