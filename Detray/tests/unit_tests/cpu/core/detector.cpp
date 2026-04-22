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
#include "detray/core/detector.hpp"

#include "detray/definitions/indexing.hpp"
#include "detray/material/predefined_materials.hpp"

// Detray test include(s)
#include "detray/test/framework/types.hpp"
#include "detray/test/utils/prefill_detector.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

/// This tests the functionality of a detector as a data store manager
GTEST_TEST(detray_core, detector) {
  using namespace detray;

  using metadata_t = test::default_metadata;
  using detector_t = detector<metadata_t>;
  using mask_id = typename detector_t::masks::id;
  using material_id = typename detector_t::material::id;
  using finder_id = typename detector_t::accel::id;

  vecmem::host_memory_resource host_mr;
  detector_t d1(host_mr);
  auto geo_ctx = typename detector_t::geometry_context{};

  // Helper lambda for checking the contents of an "empty" detector object.
  auto check_empty_detector = [](auto& d) {
    EXPECT_TRUE(d.volumes().empty());
    EXPECT_TRUE(d.portals().empty());
    EXPECT_TRUE(d.transform_store().empty());
    EXPECT_TRUE(d.mask_store().template empty<mask_id::e_rectangle2D>());
    EXPECT_TRUE(d.mask_store().template empty<mask_id::e_trapezoid2D>());
    EXPECT_TRUE(d.mask_store().template empty<mask_id::e_annulus2D>());
    EXPECT_TRUE(d.mask_store().template empty<mask_id::e_cylinder2D>());
    EXPECT_TRUE(
        d.mask_store().template empty<mask_id::e_concentric_cylinder2D>());
    EXPECT_TRUE(d.mask_store().template empty<mask_id::e_ring2D>());
    EXPECT_TRUE(d.mask_store().template empty<mask_id::e_ring2D>());
    EXPECT_TRUE(d.mask_store().template empty<mask_id::e_straw_tube>());
    EXPECT_TRUE(d.mask_store().template empty<mask_id::e_drift_cell>());
    EXPECT_TRUE(
        d.material_store().template empty<material_id::e_material_slab>());
    EXPECT_TRUE(
        d.material_store().template empty<material_id::e_material_rod>());
    EXPECT_TRUE(d.accelerator_store()
                    .template empty<finder_id::e_surface_brute_force>());
    EXPECT_TRUE(d.accelerator_store()
                    .template empty<finder_id::e_surface_ring2D_grid>());
    EXPECT_TRUE(d.accelerator_store()
                    .template empty<finder_id::e_surface_cylinder2D_grid>());
    /*EXPECT_TRUE(
        d.accelerator_store().template empty<finder_id::e_irr_disc_grid>());
    EXPECT_TRUE(d.accelerator_store()
                    .template empty<finder_id::e_irr_cylinder2_grid>());*/
    EXPECT_TRUE(
        d.accelerator_store().template empty<finder_id::e_surface_default>());
  };

  // Check the empty detector object.
  check_empty_detector(d1);

  // Add some geometrical data
  prefill_detector(d1, geo_ctx);

  // Helper lambda for checking the contents of a "filled" detector object.
  auto check_filled_detector = [](auto& d) {
    // TODO: add B-field check
    EXPECT_EQ(d.volumes().size(), 1u);
    EXPECT_EQ(d.portals().size(), 3u);
    EXPECT_EQ(d.transform_store().size(), 4u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_rectangle2D>(), 1u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_rectangle2D>(), 1u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_trapezoid2D>(), 1u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_annulus2D>(), 1u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_cylinder2D>(), 0u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_concentric_cylinder2D>(),
              0u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_ring2D>(), 0u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_ring2D>(), 0u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_straw_tube>(), 0u);
    EXPECT_EQ(d.mask_store().template size<mask_id::e_drift_cell>(), 0u);
    EXPECT_EQ(d.material_store().template size<material_id::e_material_slab>(),
              2u);
    EXPECT_EQ(d.material_store().template size<material_id::e_material_rod>(),
              1u);
    EXPECT_EQ(
        d.accelerator_store().template size<finder_id::e_surface_brute_force>(),
        1u);
    EXPECT_EQ(
        d.accelerator_store().template size<finder_id::e_surface_ring2D_grid>(),
        0u);
    EXPECT_EQ(d.accelerator_store()
                  .template size<finder_id::e_surface_cylinder2D_grid>(),
              0u);
    /*EXPECT_EQ(
        d.accelerator_store().template size<finder_id::e_irr_disc_grid>(),
        0u);
    EXPECT_EQ(d.accelerator_store()
                  .template size<finder_id::e_irr_cylinder2_grid>(),
              0u);*/
    EXPECT_EQ(
        d.accelerator_store().template size<finder_id::e_surface_default>(),
        1u);
  };

  // Check the filled detector object.
  check_filled_detector(d1);

  // Move construct a detector object.
  detector_t d2{std::move(d1)};
  check_filled_detector(d2);

  // Create a new, empty detector.
  detector_t d3{host_mr};
  check_empty_detector(d3);

  // Move assign the filled detector to the empty one.
  d3 = std::move(d2);
  check_filled_detector(d3);
}
