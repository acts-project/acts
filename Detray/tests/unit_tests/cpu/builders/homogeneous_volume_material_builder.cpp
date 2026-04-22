// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Detray include(s)
#include "detray/builders/homogeneous_volume_material_builder.hpp"

#include "detray/builders/surface_factory.hpp"
#include "detray/core/detector.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/geometry/shapes/rectangle2D.hpp"
#include "detray/material/predefined_materials.hpp"

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

using scalar = detray::test::scalar;
using point3 = test::point3;

/// Unittest: Test the construction of a collection of materials
TEST(detray_builders, homogeneous_volume_material_builder) {
  using metadata_t = test::default_metadata;
  using detector_t = detector<metadata_t>;
  using transform3 = typename detector_t::transform3_type;

  constexpr auto material_id{detector_t::material::id::e_raw_material};

  vecmem::host_memory_resource host_mr;
  detector_t d(host_mr);

  EXPECT_TRUE(d.material_store().template empty<material_id>());

  // Add material to a new volume
  auto vbuilder =
      std::make_unique<volume_builder<detector_t>>(volume_id::e_cylinder);
  auto mat_builder =
      homogeneous_volume_material_builder<detector_t>{std::move(vbuilder)};

  // Add some surfaces
  using rectangle_factory = surface_factory<detector_t, rectangle2D>;
  auto sf_factory = std::make_shared<rectangle_factory>();
  sf_factory->push_back({surface_id::e_sensitive,
                         transform3(point3{0.f, 0.f, -1.f}), 1u,
                         std::vector<scalar>{10.f, 8.f}});
  mat_builder.add_surfaces(sf_factory);

  // Set volume material
  mat_builder.set_material(argon_liquid<scalar>{});

  // Add the volume to the detector
  mat_builder.build(d);

  // Test the material data
  EXPECT_EQ(d.volumes().size(), 1u);
  const auto &vol_desc = d.volumes().at(0);
  EXPECT_EQ(vol_desc.material().id(), material_id);
  EXPECT_EQ(vol_desc.material().index(), 0u);
  EXPECT_EQ(d.material_store().template size<material_id>(), 1u);
  EXPECT_EQ(d.material_store().template get<material_id>()[0],
            argon_liquid<scalar>{});
}
