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
#include "detray/builders/volume_builder.hpp"

#include "detray/builders/cuboid_portal_generator.hpp"
#include "detray/builders/detector_builder.hpp"
#include "detray/builders/surface_factory.hpp"
#include "detray/core/detector.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes.hpp"

// Detray test include(s)
#include "detray/test/framework/types.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

namespace {

using scalar = detray::test::scalar;
using point3 = detray::test::point3;

constexpr scalar tol{1e-7f};

}  // anonymous namespace

/// This tests the functionality of a surface factory
GTEST_TEST(detray_builders, surface_factory) {
  using namespace detray;

  using metadata_t = test::default_metadata;
  using detector_t = detector<metadata_t>;
  using transform3 = dtransform3D<typename detector_t::algebra_type>;

  //
  // check portal cylinder
  //
  using portal_cylinder_factory =
      surface_factory<detector_t, concentric_cylinder2D>;

  auto pt_cyl_factory = std::make_shared<portal_cylinder_factory>();

  typename portal_cylinder_factory::sf_data_collection cyl_sf_data;
  cyl_sf_data.emplace_back(surface_id::e_portal,
                           transform3(point3{0.f, 0.f, -1000.f}), 0u,
                           std::vector<scalar>{10.f, -1000.f, 1500.f});
  cyl_sf_data.emplace_back(surface_id::e_portal,
                           transform3(point3{0.f, 0.f, 1000.f}), 2u,
                           std::vector<scalar>{20.f, -1500.f, 1000.f});

  EXPECT_EQ(pt_cyl_factory->size(), 0u);
  EXPECT_TRUE(pt_cyl_factory->types().empty());
  EXPECT_TRUE(pt_cyl_factory->bounds().empty());
  EXPECT_TRUE(pt_cyl_factory->transforms().empty());
  EXPECT_TRUE(pt_cyl_factory->volume_links().empty());

  pt_cyl_factory->push_back(std::move(cyl_sf_data));
  // data should be safely added to the factory by now
  cyl_sf_data.clear();

  EXPECT_EQ(pt_cyl_factory->size(), 2u);
  EXPECT_EQ(pt_cyl_factory->types().size(), 2u);
  EXPECT_EQ(pt_cyl_factory->bounds().size(), 2u);
  EXPECT_EQ(pt_cyl_factory->transforms().size(), 2u);
  EXPECT_EQ(pt_cyl_factory->volume_links().size(), 2u);
  const auto& portal_cyl_comps = pt_cyl_factory->bounds().front().front();
  EXPECT_NEAR(portal_cyl_comps[0], 10.f, tol);
  EXPECT_NEAR(portal_cyl_comps[1], -1000.f, tol);
  EXPECT_NEAR(portal_cyl_comps[2], 1500.f, tol);
  const auto& portal_cyl_vol_links = pt_cyl_factory->volume_links();
  EXPECT_EQ(portal_cyl_vol_links[0][0], 0u);
  EXPECT_EQ(portal_cyl_vol_links[1][0], 2u);

  //
  // check sensitive cylinder
  //
  using sensitive_cylinder_factory = surface_factory<detector_t, cylinder2D>;

  auto cyl_factory = std::make_shared<sensitive_cylinder_factory>();

  cyl_sf_data.emplace_back(surface_id::e_sensitive,
                           transform3(point3{0.f, 0.f, -50.f}), 1u,
                           std::vector<scalar>{5.f, -900.f, 900.f});
  cyl_sf_data.emplace_back(surface_id::e_sensitive,
                           transform3(point3{0.f, 0.f, 50.f}), 1u,
                           std::vector<scalar>{5.f, -900.f, 900.f});

  cyl_sf_data.emplace_back(surface_id::e_passive,
                           transform3(point3{0.f, 0.f, -20.f}), 1u,
                           std::vector<scalar>{4.9f, -900.f, 900.f});
  cyl_sf_data.emplace_back(surface_id::e_passive,
                           transform3(point3{0.f, 0.f, 20.f}), 1u,
                           std::vector<scalar>{4.9f, -900.f, 900.f});

  EXPECT_EQ(cyl_factory->size(), 0u);
  EXPECT_TRUE(cyl_factory->types().empty());
  EXPECT_TRUE(cyl_factory->bounds().empty());
  EXPECT_TRUE(cyl_factory->transforms().empty());
  EXPECT_TRUE(cyl_factory->volume_links().empty());

  cyl_factory->push_back(std::move(cyl_sf_data));
  cyl_sf_data.clear();

  EXPECT_EQ(cyl_factory->size(), 4u);
  EXPECT_EQ(cyl_factory->types().size(), 4u);
  EXPECT_EQ(cyl_factory->bounds().size(), 4u);
  EXPECT_EQ(cyl_factory->transforms().size(), 4u);
  EXPECT_EQ(cyl_factory->volume_links().size(), 4u);
  const auto& sens_cyl_comps = cyl_factory->bounds().front().front();
  EXPECT_NEAR(sens_cyl_comps[0], 5.f, tol);
  EXPECT_NEAR(sens_cyl_comps[1], -900.f, tol);
  EXPECT_NEAR(sens_cyl_comps[2], 900.f, tol);
  const auto& sens_cyl_vol_links = cyl_factory->volume_links().front();
  EXPECT_EQ(sens_cyl_vol_links[0], 1u);

  //
  // check the other mask types for this detector
  //

  // annulus
  using annulus_factory = surface_factory<detector_t, annulus2D>;

  auto ann_factory = std::make_shared<annulus_factory>();

  typename annulus_factory::sf_data_collection ann_sf_data;
  ann_sf_data.emplace_back(
      surface_id::e_sensitive, transform3(point3{0.f, 0.f, 0.f}), 1u,
      std::vector<scalar>{300.f, 350.f, -0.1f, 0.1f, 0.5f, 0.6f, 1.4f});
  ann_factory->push_back(std::move(ann_sf_data));
  ann_sf_data.clear();

  const auto& ann_comps = ann_factory->bounds().front().front();
  EXPECT_NEAR(ann_comps[0], 300.f, tol);
  EXPECT_NEAR(ann_comps[1], 350.f, tol);
  EXPECT_NEAR(ann_comps[2], -0.1f, tol);
  EXPECT_NEAR(ann_comps[3], 0.1f, tol);
  EXPECT_NEAR(ann_comps[4], 0.5f, tol);
  EXPECT_NEAR(ann_comps[5], 0.6f, tol);
  EXPECT_NEAR(ann_comps[6], 1.4f, tol);

  // rectangles
  using rectangle_factory = surface_factory<detector_t, rectangle2D>;

  auto rect_factory = std::make_shared<rectangle_factory>();

  typename rectangle_factory::sf_data_collection rect_sf_data;
  rect_sf_data.emplace_back(surface_id::e_sensitive,
                            transform3(point3{0.f, 0.f, 0.f}), 1u,
                            std::vector<scalar>{10.f, 8.f});
  rect_factory->push_back(std::move(rect_sf_data));
  rect_sf_data.clear();

  const auto& rectgl_comps = rect_factory->bounds().front().front();
  EXPECT_NEAR(rectgl_comps[0], 10.f, tol);
  EXPECT_NEAR(rectgl_comps[1], 8.f, tol);

  // ring
  using disc_factory = surface_factory<detector_t, ring2D>;

  auto sf_disc_factory = std::make_shared<disc_factory>();

  typename disc_factory::sf_data_collection disc_sf_data;
  disc_sf_data.emplace_back(surface_id::e_passive,
                            transform3(point3{0.f, 0.f, 0.f}), 1u,
                            std::vector<scalar>{0.f, 5.f});
  sf_disc_factory->push_back(std::move(disc_sf_data));
  disc_sf_data.clear();

  const auto& ring_comps = sf_disc_factory->bounds().front().front();
  EXPECT_NEAR(ring_comps[0], 0.f, tol);
  EXPECT_NEAR(ring_comps[1], 5.f, tol);

  // trapezoid
  using trapezoid_factory = surface_factory<detector_t, trapezoid2D>;

  auto trpz_factory = std::make_shared<trapezoid_factory>();

  typename trapezoid_factory::sf_data_collection trpz_sf_data;
  trpz_sf_data.emplace_back(surface_id::e_sensitive,
                            transform3(point3{0.f, 0.f, 0.f}), 1u,
                            std::vector<scalar>{1.f, 3.f, 2.f, 0.25f});
  trpz_factory->push_back(std::move(trpz_sf_data));
  trpz_sf_data.clear();

  const auto& trpz_comps = trpz_factory->bounds().front().front();
  EXPECT_NEAR(trpz_comps[0], 1.f, tol);
  EXPECT_NEAR(trpz_comps[1], 3.f, tol);
  EXPECT_NEAR(trpz_comps[2], 2.f, tol);
  EXPECT_NEAR(trpz_comps[3], 0.25f, tol);
}

/// This tests the initialization of a detector volume using a volume builder
GTEST_TEST(detray_builders, volume_builder) {
  using namespace detray;

  vecmem::host_memory_resource host_mr;

  using metadata_t = test::default_metadata;
  using detector_t = detector<metadata_t>;
  using transform3 = typename detector_t::transform3_type;

  detector_t d(host_mr);

  EXPECT_TRUE(d.volumes().size() == 0u);

  using rectangle_factory = surface_factory<detector_t, rectangle2D>;
  auto sf_factory = std::make_shared<rectangle_factory>();
  sf_factory->push_back({surface_id::e_sensitive,
                         transform3(point3{0.f, 0.f, -1.f}), 1u,
                         std::vector<scalar>{10.f, 8.f}});

  volume_builder<detector_t> vbuilder{volume_id::e_cylinder};
  vbuilder.add_surfaces(sf_factory);
  vbuilder.build(d);

  const auto& vol = d.volumes().back();

  EXPECT_TRUE(d.volumes().size() == 1u);
  EXPECT_EQ(vol.index(), 0u);
  EXPECT_EQ(vol.index(), vbuilder.vol_index());
  EXPECT_EQ(vol.id(), volume_id::e_cylinder);
}
