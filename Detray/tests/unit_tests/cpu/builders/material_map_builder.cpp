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
#include "detray/builders/material_map_builder.hpp"

#include "detray/builders/cuboid_portal_generator.hpp"
#include "detray/builders/detector_builder.hpp"
#include "detray/builders/material_map_factory.hpp"
#include "detray/builders/surface_factory.hpp"
#include "detray/builders/volume_builder.hpp"
#include "detray/core/detector.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/utils/ranges.hpp"

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
using mat_id = typename detector_t::material::id;
using bin_index_t = axis::multi_bin<2u>;

/// Add generate input material for material maps
template <typename material_factory_t, concepts::scalar scalar_t>
auto add_material_data(const material_factory_t& mat_factory, mat_id id,
                       std::size_t sf_index, scalar_t t,
                       material<scalar_t> mat = silicon<scalar_t>()) {
  typename material_factory_t::element_type::data_type mat_data{sf_index};
  std::vector<bin_index_t> m_bins{};
  std::vector<std::size_t> n_bins{5u, 10u};
  std::vector<std::vector<scalar_t>> axis_spans = {};

  // Add material for every bin
  for (auto [i, j] : detray::views::cartesian_product{
           detray::views::iota{0u, 5u}, detray::views::iota{0u, 10u}}) {
    m_bins.push_back({i, j});
    mat_data.append(t, mat);
    t += 0.25f * unit<scalar_t>::mm;
  }

  mat_factory->add_material(id, std::move(mat_data), std::move(n_bins),
                            std::move(axis_spans), std::move(m_bins));
}

}  // anonymous namespace

/// Unittest: Test the construction of a collection of material maps
TEST(detray_builders, material_map_factory) {
  using algebra_t = typename detector_t::algebra_type;
  using transform3 = dtransform3D<algebra_t>;

  // Build rectangle surfaces with material slabs
  using rectangle_factory = surface_factory<detector_t, rectangle2D>;
  auto mat_factory =
      std::make_unique<material_map_factory<detector_t, bin_index_t>>(
          std::make_unique<rectangle_factory>());

  EXPECT_EQ(mat_factory->size(), 0u);
  EXPECT_TRUE(mat_factory->materials().empty());
  EXPECT_TRUE(mat_factory->thickness().empty());

  // Add material maps for a few rectangle surfaces
  mat_factory->push_back({surface_id::e_sensitive,
                          transform3(point3{0.f, 0.f, -1.f}), 1u,
                          std::vector<scalar>{10.f, 8.f}});

  scalar t{1.f * unit<scalar>::mm};
  add_material_data(mat_factory, mat_id::e_rectangle2D_map, 0u, t,
                    silicon<scalar>());

  mat_factory->push_back({surface_id::e_sensitive,
                          transform3(point3{0.f, 0.f, 1.f}), 1u,
                          std::vector<scalar>{20.f, 16.f}});
  t = 2.f * unit<scalar>::mm;
  add_material_data(mat_factory, mat_id::e_rectangle2D_map, 1u, t,
                    tungsten<scalar>());

  // Pass the parameters for 'gold'
  mat_factory->push_back({surface_id::e_sensitive,
                          transform3(point3{0.f, 0.f, 1.f}), 1u,
                          std::vector<scalar>{20.f, 16.f}});
  t = 3.f * unit<scalar>::mm;
  add_material_data(
      mat_factory, mat_id::e_rectangle2D_map, 2u, t,
      {3.344f * unit<scalar>::mm, 101.6f * unit<scalar>::mm, 196.97f, 79,
       19.32f * unit<scalar>::g / (1.f * unit<scalar>::cm3),
       material_state::e_solid});

  EXPECT_EQ(mat_factory->size(), 3u);

  // Test the material data
  EXPECT_EQ(mat_factory->links().size(), 3u);
  EXPECT_EQ(mat_factory->thickness().at(0).size(), 50u);
  EXPECT_EQ(mat_factory->thickness().at(1).size(), 50u);
  EXPECT_EQ(mat_factory->thickness().at(2).size(), 50u);
  EXPECT_EQ(mat_factory->materials().at(0).front(), silicon<scalar>());
  EXPECT_EQ(mat_factory->materials().at(1).front(), tungsten<scalar>());
  EXPECT_EQ(mat_factory->materials().at(2).front(), gold<scalar>());
}
