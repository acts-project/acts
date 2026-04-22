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

// Project include(s)
#include "detray/navigation/detail/intersection_kernel.hpp"

#include "detray/core/detail/multi_store.hpp"
#include "detray/core/detail/single_store.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes.hpp"
#include "detray/geometry/surface_descriptor.hpp"
#include "detray/navigation/intersection/helix_intersector.hpp"
#include "detray/navigation/intersection/ray_intersector.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/tracks/trajectories.hpp"
#include "detray/utils/ranges.hpp"

// Detray test include(s)
#include "detray/test/framework/types.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// Google Test include(s)
#include <gtest/gtest.h>

using namespace detray;

using test_algebra = test::algebra;
using scalar = test::scalar;
using vector3 = test::vector3;
using point3 = test::point3;
using transform_link_t = dindex;

namespace {

constexpr const scalar is_close{1e-5f};
constexpr const scalar external_tol{1.f * unit<scalar>::mm};

const intersection::config intr_cfg{
    .min_mask_tolerance = 1e-3f * unit<float>::mm,
    .max_mask_tolerance = 1e-3f * unit<float>::mm,
    .mask_tolerance_scalor = 0.f,
    .overstep_tolerance = 0.f};

enum class mask_id : unsigned int {
  e_rectangle2D = 0u,
  e_trapezoid2D = 1u,
  e_annulus2D = 2u,
  e_cylinder2D = 3u,
  e_concentric_cylinder2D = 4u,
};

enum class material_id : unsigned int {
  e_material_slab = 0u,
};

}  // anonymous namespace

/// Surface components:
using volume_link_t = std::uint_least16_t;

/// - masks, with mask identifiers 0,1,2
using rectangle_t = mask<rectangle2D, test_algebra, volume_link_t>;
using trapezoid_t = mask<trapezoid2D, test_algebra, volume_link_t>;
using annulus_t = mask<annulus2D, test_algebra, volume_link_t>;
using cylinder_t = mask<cylinder2D, test_algebra, volume_link_t>;
using cylinder_portal_t =
    mask<concentric_cylinder2D, test_algebra, volume_link_t>;

using mask_container_t =
    regular_multi_store<mask_id, empty_context, dtuple, dvector, rectangle_t,
                        trapezoid_t, annulus_t, cylinder_t, cylinder_portal_t>;
using mask_link_t = typename mask_container_t::single_link;
using material_link_t = dtyped_index<material_id, dindex>;

using transform_container_t =
    single_store<test::transform3, dvector, geometry_context>;

/// The Surface definition:
using surface_t =
    surface_descriptor<mask_link_t, material_link_t, transform_link_t>;
using surface_container_t = dvector<surface_t>;

// TODO: How about merging ray and helix tests into one to remove the code
// repetition?

// This tests the construction of a surface
GTEST_TEST(detray_intersection, intersection_kernel_ray) {
  vecmem::host_memory_resource host_mr;

  // The transforms & their store
  typename transform_container_t::context_type static_context{};
  transform_container_t transform_store;
  // Transforms of the rectangle, trapezoid, annulus and the two cylinders
  transform_store.emplace_back(static_context, vector3{0.f, 0.f, 10.f});
  transform_store.emplace_back(static_context, vector3{0.f, 0.f, 20.f});
  transform_store.emplace_back(static_context, vector3{0.f, -20.f, 30.f});
  // 90deg rotation around y-axis
  transform_store.emplace_back(static_context, vector3{0.f, 0.f, 50.f},
                               vector3{1.f, 0.f, 0.f}, vector3{0.f, 0.f, -1.f});
  // Identity transform for concentric cylinder
  transform_store.emplace_back(static_context, vector3{0.f, 0.f, 0.f});
  // Shifted rectangle
  transform_store.emplace_back(static_context, vector3{12.f, 0.f, 1000.f});

  // The masks & their store
  mask_container_t mask_store(host_mr);

  const rectangle_t rect{0u, 10.f, 10.f};
  const trapezoid_t trap{0u, 10.f, 20.f, 30.f, 1.f / 60.f};
  const annulus_t annl{0u, 15.f, 55.f, 0.75f, 1.95f, 0.f, 2.f, -2.f};
  const cylinder_t cyl{0u, 5.f, -10.f, 10.f};
  const cylinder_portal_t cyl_portal{0u, 1.f, 0.f, 1000.f};

  mask_store.template push_back<mask_id::e_rectangle2D>(rect, empty_context{});
  mask_store.template push_back<mask_id::e_rectangle2D>(rect, empty_context{});
  mask_store.template push_back<mask_id::e_trapezoid2D>(trap, empty_context{});
  mask_store.template push_back<mask_id::e_annulus2D>(annl, empty_context{});
  mask_store.template push_back<mask_id::e_cylinder2D>(cyl, empty_context{});
  mask_store.template push_back<mask_id::e_concentric_cylinder2D>(
      cyl_portal, empty_context{});

  // The surfaces and their store
  surface_t rectangle_surface(0u, {mask_id::e_rectangle2D, 0u},
                              {material_id::e_material_slab, 0u}, 0u,
                              surface_id::e_sensitive);
  surface_t trapezoid_surface(1u, {mask_id::e_trapezoid2D, 0u},
                              {material_id::e_material_slab, 1u}, 0u,
                              surface_id::e_sensitive);
  surface_t annulus_surface(2u, {mask_id::e_annulus2D, 0u},
                            {material_id::e_material_slab, 2u}, 0u,
                            surface_id::e_sensitive);
  surface_t cyl_surface(3u, {mask_id::e_cylinder2D, 0u},
                        {material_id::e_material_slab, 2u}, 0u,
                        surface_id::e_passive);
  surface_t cyl_portal_surface(4u, {mask_id::e_concentric_cylinder2D, 0u},
                               {material_id::e_material_slab, 2u}, 0u,
                               surface_id::e_portal);
  surface_t rectangle_shifted_surface(5u, {mask_id::e_rectangle2D, 1u},
                                      {material_id::e_material_slab, 0u}, 0u,
                                      surface_id::e_sensitive);
  // Easier test debugging
  rectangle_surface.set_index(0u);
  trapezoid_surface.set_index(1u);
  annulus_surface.set_index(2u);
  cyl_surface.set_index(3u);
  cyl_portal_surface.set_index(4u);
  rectangle_shifted_surface.set_index(5u);

  surface_container_t surfaces = {
      rectangle_surface, trapezoid_surface,  annulus_surface,
      cyl_surface,       cyl_portal_surface, rectangle_shifted_surface};

  const point3 pos{0.f, 0.f, 0.f};
  const vector3 mom{0.01f, 0.01f, 10.f};
  const free_track_parameters<test_algebra> track(pos, 0.f, mom, -1.f);

  // Validation data
  const point3 expected_rectangle{0.01f, 0.01f, 10.f};
  const point3 expected_trapezoid{0.02f, 0.02f, 20.f};
  const point3 expected_annulus{0.03f, 0.03f, 30.f};
  const point3 expected_cylinder1{0.045f, 0.045f, 45.0f};
  const point3 expected_cylinder2{0.055f, 0.055f, 55.0f};
  const point3 expected_cylinder_pt{0.7071068f, 0.7071068f, 707.1068f};
  const point3 expected_rectangle_shifted{1.f, 1.f, 1000.f};

  const std::vector<point3> expected_points = {
      expected_rectangle,        expected_trapezoid, expected_annulus,
      expected_cylinder1,        expected_cylinder2, expected_cylinder_pt,
      expected_rectangle_shifted};

  // Initialize kernel
  std::vector<
      intersection2D<surface_t, test_algebra, intersection::contains_pos>>
      sfi_init;
  sfi_init.reserve(expected_points.size());

  for (const auto &surface : surfaces) {
    mask_store.visit<detail::intersection_initialize<ray_intersector>>(
        surface.mask(), sfi_init, detail::ray(track), surface, transform_store,
        static_context, intr_cfg, external_tol);
  }

  EXPECT_EQ(expected_points.size(), sfi_init.size());

  // Also check intersections
  for (std::size_t i = 0u;
       i < math::min(expected_points.size(), sfi_init.size()); ++i) {
    EXPECT_TRUE(sfi_init[i].is_along());
    EXPECT_EQ(sfi_init[i].volume_link(), 0u);

    vector3 global{0.f, 0.f, 0.f};

    const surface_t sf_desc = sfi_init[i].surface();
    if (sf_desc.mask().id() == mask_id::e_rectangle2D) {
      global = rect.to_global_frame(transform_store.at(sf_desc.transform()),
                                    sfi_init[i].local());
    } else if (sf_desc.mask().id() == mask_id::e_trapezoid2D) {
      global = trap.to_global_frame(transform_store.at(sf_desc.transform()),
                                    sfi_init[i].local());
    } else if (sf_desc.mask().id() == mask_id::e_annulus2D) {
      global = annl.to_global_frame(transform_store.at(sf_desc.transform()),
                                    sfi_init[i].local());
    } else if (sf_desc.mask().id() == mask_id::e_cylinder2D) {
      global = cyl.to_global_frame(transform_store.at(sf_desc.transform()),
                                   sfi_init[i].local());
    } else if (sf_desc.mask().id() == mask_id::e_concentric_cylinder2D) {
      global = cyl_portal.to_global_frame(
          transform_store.at(sf_desc.transform()), sfi_init[i].local());
    }

    EXPECT_NEAR(global[0], expected_points[i][0], 1e-3f)
        << " at point " << i << ": surface " << sf_desc.identifier();
    EXPECT_NEAR(global[1], expected_points[i][1], 1e-3f)
        << " at point " << i << ": surface " << sf_desc.identifier();
    EXPECT_NEAR(global[2], expected_points[i][2], 1e-3f)
        << " at point " << i << ": surface " << sf_desc.identifier();
  }

  // Update kernel
  // @fixme: The intersection update kernel does not work for non-portal
  // cylinders, since it assigns the closest intersection to both
  // solutions
  /*std::vector<intersection2D_point<surface_t, test_algebra>> sfi_update;
  sfi_update.resize(5);

  for (const auto [idx, surface] : detray::views::enumerate(surfaces)) {
      sfi_update[idx].set_surface(surface);
      mask_store.visit<detail::intersection_update>(
          surface.mask(), detail::ray(track), sfi_update[idx],
          transform_store);

      if(!sfi_update[idx].is_inside()) {
  continue; } ASSERT_TRUE(sfi_update[idx].is_along()) << " at surface " <<
  sfi_update[idx]
  << ", " << sfi_init[idx]; ASSERT_EQ(sfi_update[idx].volume_link(), 0u);
      ASSERT_NEAR(sfi_update[idx].p3[0], expected_points[idx][0],
  is_close)
      << " at surface " << sfi_update[idx] << ", " << sfi_init[idx];
      ASSERT_NEAR(sfi_update[idx].p3[1], expected_points[idx][1],
  is_close) << " at surface " << sfi_update[idx] << ", " << sfi_init[idx];
      ASSERT_NEAR(sfi_update[idx].p3[2], expected_points[idx][2],
  is_close)
      << " at surface " << sfi_update[idx] << ", " << sfi_init[idx];
  }

  // Compare
  ASSERT_EQ(sfi_init.size(), 5u);
  ASSERT_EQ(sfi_update.size(), 5u);
  for (unsigned int i = 0u; i < 5u; i++) {
      ASSERT_EQ(sfi_init[i].p3, sfi_update[i].p3);
      ASSERT_EQ(sfi_init[i].local(), sfi_update[i].local());
      ASSERT_EQ(sfi_init[i].path(), sfi_update[i].path());
  }*/
}

/// Reuse the intersection kernel test for particle gun
GTEST_TEST(detray_intersection, intersection_kernel_helix) {
  vecmem::host_memory_resource host_mr;

  // The transforms & their store
  typename transform_container_t::context_type static_context{};
  transform_container_t transform_store;
  // Transforms of the rectangle, trapezoid and annulus
  transform_store.emplace_back(static_context, vector3{0.f, 0.f, 10.f});
  transform_store.emplace_back(static_context, vector3{0.f, 0.f, 20.f});
  transform_store.emplace_back(static_context, vector3{0.f, -20.f, 30.f});
  // The masks & their store
  mask_container_t mask_store(host_mr);

  const rectangle_t rect{0u, 10.f, 10.f};
  const trapezoid_t trap{0u, 10.f, 20.f, 30.f, 1.f / 60.f};
  const annulus_t annl{0u, 15.f, 55.f, 0.75f, 1.95f, 0.f, 2.f, -2.f};
  mask_store.template push_back<mask_id::e_rectangle2D>(rect, empty_context{});
  mask_store.template push_back<mask_id::e_trapezoid2D>(trap, empty_context{});
  mask_store.template push_back<mask_id::e_annulus2D>(annl, empty_context{});

  // The surfaces and their store
  const surface_t rectangle_surface(0u, {mask_id::e_rectangle2D, 0u},
                                    {material_id::e_material_slab, 0u}, 0u,
                                    surface_id::e_sensitive);
  const surface_t trapezoid_surface(1u, {mask_id::e_trapezoid2D, 0u},
                                    {material_id::e_material_slab, 1u}, 0u,
                                    surface_id::e_sensitive);
  const surface_t annulus_surface(2u, {mask_id::e_annulus2D, 0u},
                                  {material_id::e_material_slab, 2u}, 0u,
                                  surface_id::e_sensitive);
  surface_container_t surfaces = {rectangle_surface, trapezoid_surface,
                                  annulus_surface};
  const point3 pos{0.f, 0.f, 0.f};
  const vector3 mom{0.01f, 0.01f, 10.f};
  const vector3 B{0.f * unit<scalar>::T, 0.f * unit<scalar>::T,
                  is_close * unit<scalar>::T};
  const detail::helix<test_algebra> h({pos, 0.f, mom, -1.f}, B);

  // Validation data
  const point3 expected_rectangle{0.01f, 0.01f, 10.f};
  const point3 expected_trapezoid{0.02f, 0.02f, 20.f};
  const point3 expected_annulus{0.03f, 0.03f, 30.f};
  const std::vector<point3> expected_points = {
      expected_rectangle, expected_trapezoid, expected_annulus};
  std::vector<
      intersection2D<surface_t, test_algebra, intersection::contains_pos>>
      sfi_helix{};
  sfi_helix.reserve(expected_points.size());

  // Try the intersections - with automated dispatching via the kernel
  for (const auto [sf_idx, surface] : detray::views::enumerate(surfaces)) {
    mask_store.visit<detail::intersection_initialize<helix_intersector>>(
        surface.mask(), sfi_helix, h, surface, transform_store, static_context,
        intr_cfg, scalar{0.f});

    vector3 global{0.f, 0.f, 0.f};

    if (surface.mask().id() == mask_id::e_rectangle2D) {
      global =
          rect.to_global_frame(transform_store.at(0), sfi_helix[0].local());
    } else if (surface.mask().id() == mask_id::e_trapezoid2D) {
      global =
          trap.to_global_frame(transform_store.at(1), sfi_helix[0].local());
    } else if (surface.mask().id() == mask_id::e_annulus2D) {
      global =
          annl.to_global_frame(transform_store.at(2), sfi_helix[0].local());
    }

    ASSERT_NEAR(global[0], expected_points[sf_idx][0], is_close);
    ASSERT_NEAR(global[1], expected_points[sf_idx][1], is_close);
    ASSERT_NEAR(global[2], expected_points[sf_idx][2], is_close);

    sfi_helix.clear();
  }
}
