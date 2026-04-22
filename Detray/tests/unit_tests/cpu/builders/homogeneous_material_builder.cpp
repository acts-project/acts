// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Detray include(s)
#include "detray/builders/homogeneous_material_builder.hpp"

#include "detray/builders/cuboid_portal_generator.hpp"
#include "detray/builders/homogeneous_material_factory.hpp"
#include "detray/builders/surface_factory.hpp"
#include "detray/builders/volume_builder.hpp"
#include "detray/core/detector.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/geometry/shapes/rectangle2D.hpp"

// Detray test include(s)
#include "detray/test/framework/types.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <limits>
#include <memory>
#include <vector>

using namespace detray;

namespace {

using scalar = test::scalar;
using point3 = test::point3;

using metadata_t = test::default_metadata;
using detector_t = detector<metadata_t>;

constexpr scalar tol{std::numeric_limits<scalar>::epsilon()};

}  // anonymous namespace

/// Unittest: Test the construction of a collection of materials
TEST(detray_builders, homogeneous_material_factory) {
  using transform3 = typename detector_t::transform3_type;
  using material_id = typename detector_t::material::id;

  // Build rectangle surfaces with material slabs
  using rectangle_factory = surface_factory<detector_t, rectangle2D>;
  auto mat_factory = std::make_unique<homogeneous_material_factory<detector_t>>(
      std::make_unique<rectangle_factory>());

  vecmem::host_memory_resource host_mr;
  detector_t d(host_mr);

  EXPECT_TRUE(
      d.material_store().template empty<material_id::e_material_slab>());
  EXPECT_TRUE(d.material_store().template empty<material_id::e_material_rod>());

  EXPECT_EQ(mat_factory->size(), 0u);
  EXPECT_TRUE(mat_factory->materials().empty());
  EXPECT_TRUE(mat_factory->thickness().empty());

  // Add material for a few rectangle surfaces
  mat_factory->push_back({surface_id::e_sensitive,
                          transform3(point3{0.f, 0.f, -1.f}), 1u,
                          std::vector<scalar>{10.f, 8.f}});
  mat_factory->add_material(material_id::e_material_slab,
                            {1.f * unit<scalar>::mm, silicon<scalar>()});
  mat_factory->push_back({surface_id::e_sensitive,
                          transform3(point3{0.f, 0.f, 1.f}), 1u,
                          std::vector<scalar>{20.f, 16.f}});
  mat_factory->add_material(material_id::e_material_slab,
                            {10.f * unit<scalar>::mm, tungsten<scalar>()});
  // Pass the parameters for 'gold'
  mat_factory->push_back({surface_id::e_sensitive,
                          transform3(point3{0.f, 0.f, 1.f}), 1u,
                          std::vector<scalar>{20.f, 16.f}});
  mat_factory->add_material(
      material_id::e_material_rod,
      {0.1f * unit<scalar>::mm,
       std::vector<scalar>{
           3.344f * unit<scalar>::mm, 101.6f * unit<scalar>::mm, 196.97f, 79,
           19.32f * unit<scalar>::g / (1.f * unit<scalar>::cm3)},
       material_state::e_solid});

  EXPECT_EQ(mat_factory->size(), 3u);

  // Test the material data
  EXPECT_NEAR(mat_factory->thickness()[0], 1.f * unit<scalar>::mm, tol);
  EXPECT_NEAR(mat_factory->thickness()[1], 10.f * unit<scalar>::mm, tol);
  EXPECT_NEAR(mat_factory->thickness()[2], 0.1f * unit<scalar>::mm, tol);
  EXPECT_EQ(mat_factory->materials()[0], silicon<scalar>());
  EXPECT_EQ(mat_factory->materials()[1], tungsten<scalar>());
  EXPECT_EQ(mat_factory->materials()[2], gold<scalar>());
}
