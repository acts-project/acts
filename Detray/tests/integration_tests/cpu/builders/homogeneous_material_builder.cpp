// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Detray include(s)
#include "detray/builders/homogeneous_material_builder.hpp"

#include "detray/builders/cuboid_portal_generator.hpp"
#include "detray/builders/detector_builder.hpp"
#include "detray/builders/homogeneous_material_factory.hpp"
#include "detray/builders/surface_factory.hpp"
#include "detray/builders/volume_builder.hpp"
#include "detray/core/detector.hpp"
#include "detray/definitions/indexing.hpp"

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

/// Integration test: material builder as volume builder decorator
GTEST_TEST(detray_builders, decorator_homogeneous_material_builder) {
  using transform3 = typename detector_t::transform3_type;
  using mask_id = typename detector_t::masks::id;
  using material_id = typename detector_t::material::id;

  using pt_cylinder_t = concentric_cylinder2D;
  using pt_cylinder_factory_t = surface_factory<detector_t, pt_cylinder_t>;
  using rectangle_factory = surface_factory<detector_t, rectangle2D>;
  using trapezoid_factory = surface_factory<detector_t, trapezoid2D>;
  using cylinder_factory = surface_factory<detector_t, cylinder2D>;

  using mat_factory_t = homogeneous_material_factory<detector_t>;

  vecmem::host_memory_resource host_mr;
  detector_t d(host_mr);
  auto geo_ctx = typename detector_t::geometry_context{};

  auto vbuilder =
      std::make_unique<volume_builder<detector_t>>(volume_id::e_cylinder);
  auto mat_builder =
      homogeneous_material_builder<detector_t>{std::move(vbuilder)};

  EXPECT_TRUE(d.volumes().empty());

  // Add some portals first
  auto pt_cyl_factory = std::make_unique<pt_cylinder_factory_t>();

  pt_cyl_factory->push_back({surface_id::e_portal,
                             transform3(point3{0.f, 0.f, 0.f}), 1u,
                             std::vector<scalar>{10.f, -1500.f, 1500.f}});
  pt_cyl_factory->push_back({surface_id::e_portal,
                             transform3(point3{0.f, 0.f, 0.f}), 2u,
                             std::vector<scalar>{20.f, -1500.f, 1500.f}});

  // Then some passive and sensitive surfaces
  auto rect_factory = std::make_unique<rectangle_factory>();

  typename rectangle_factory::sf_data_collection rect_sf_data;
  rect_sf_data.emplace_back(surface_id::e_sensitive,
                            transform3(point3{0.f, 0.f, -10.f}), 0u,
                            std::vector<scalar>{10.f, 8.f});
  rect_sf_data.emplace_back(surface_id::e_sensitive,
                            transform3(point3{0.f, 0.f, -20.f}), 0u,
                            std::vector<scalar>{10.f, 8.f});
  rect_sf_data.emplace_back(surface_id::e_sensitive,
                            transform3(point3{0.f, 0.f, -30.f}), 0u,
                            std::vector<scalar>{10.f, 8.f});
  rect_factory->push_back(std::move(rect_sf_data));

  auto trpz_factory = std::make_unique<trapezoid_factory>();

  trpz_factory->push_back({surface_id::e_sensitive,
                           transform3(point3{0.f, 0.f, 1000.f}), 0u,
                           std::vector<scalar>{1.f, 3.f, 2.f, 0.25f}});

  auto cyl_factory = std::make_unique<cylinder_factory>();

  cyl_factory->push_back({surface_id::e_passive,
                          transform3(point3{0.f, 0.f, 0.f}), 0u,
                          std::vector<scalar>{5.f, -1300.f, 1300.f}});

  // Now add the material for each surface
  auto mat_pt_cyl_factory =
      std::make_shared<mat_factory_t>(std::move(pt_cyl_factory));
  mat_pt_cyl_factory->add_material(material_id::e_material_slab,
                                   {1.f * unit<scalar>::mm, silicon<scalar>()});
  mat_pt_cyl_factory->add_material(
      material_id::e_material_slab,
      {1.5f * unit<scalar>::mm, silicon<scalar>()});

  auto mat_rect_factory =
      std::make_shared<mat_factory_t>(std::move(rect_factory));
  mat_rect_factory->add_material(material_id::e_material_slab,
                                 {1.f * unit<scalar>::mm, silicon<scalar>()});
  mat_rect_factory->add_material(material_id::e_material_slab,
                                 {2.f * unit<scalar>::mm, silicon<scalar>()});
  mat_rect_factory->add_material(material_id::e_material_slab,
                                 {3.f * unit<scalar>::mm, silicon<scalar>()});

  auto mat_trpz_factory =
      std::make_shared<mat_factory_t>(std::move(trpz_factory));
  mat_trpz_factory->add_material(material_id::e_material_slab,
                                 {1.f * unit<scalar>::mm, silicon<scalar>()});

  auto mat_cyl_factory =
      std::make_shared<mat_factory_t>(std::move(cyl_factory));
  mat_cyl_factory->add_material(material_id::e_material_slab,
                                {1.5f * unit<scalar>::mm, tungsten<scalar>()});

  // Add surfaces and material to detector
  mat_builder.add_surfaces(mat_pt_cyl_factory, geo_ctx);
  mat_builder.add_surfaces(mat_rect_factory, geo_ctx);
  mat_builder.add_surfaces(mat_trpz_factory, geo_ctx);
  mat_builder.add_surfaces(mat_cyl_factory, geo_ctx);

  // Add the volume to the detector
  mat_builder.build(d);

  //
  // check results
  //
  const auto &vol = d.volumes().back();
  EXPECT_TRUE(d.volumes().size() == 1u);
  EXPECT_EQ(vol.index(), 0u);
  EXPECT_EQ(vol.id(), volume_id::e_cylinder);

  EXPECT_EQ(d.surfaces().size(), 7u);
  EXPECT_EQ(d.transform_store().size(), 8u);
  EXPECT_EQ(d.mask_store().template size<mask_id::e_concentric_cylinder2D>(),
            2u);
  EXPECT_EQ(d.mask_store().template size<mask_id::e_ring2D>(), 0u);
  EXPECT_EQ(d.mask_store().template size<mask_id::e_cylinder2D>(), 1u);
  EXPECT_EQ(d.mask_store().template size<mask_id::e_rectangle2D>(), 3u);
  EXPECT_EQ(d.mask_store().template size<mask_id::e_trapezoid2D>(), 1u);

  EXPECT_EQ(d.material_store().template size<material_id::e_material_slab>(),
            7u);
  EXPECT_EQ(d.material_store().template size<material_id::e_material_rod>(),
            0u);

  for (auto [idx, sf_desc] : detray::views::enumerate(d.surfaces())) {
    const auto &mat_link = sf_desc.material();
    EXPECT_EQ(mat_link.id(), material_id::e_material_slab);
    EXPECT_EQ(mat_link.index(), idx);
  }

  for (const auto &mat_slab :
       d.material_store().template get<material_id::e_material_slab>()) {
    EXPECT_TRUE(mat_slab.get_material() == silicon<scalar>() ||
                mat_slab.get_material() == tungsten<scalar>());
  }
}

/// Integration test to build an empty cuboid volume with material
GTEST_TEST(detray_builders, detector_builder_with_material) {
  using namespace detray;

  using transform3 = typename detector_t::transform3_type;
  using mask_id = typename detector_t::masks::id;
  using material_id = typename detector_t::material::id;

  // Surface factories
  using trapezoid_factory = surface_factory<detector_t, trapezoid2D>;

  // detector builder
  detector_builder<typename detector_t::metadata> det_builder{};
  auto geo_ctx = typename detector_t::geometry_context{};

  // Vanilla volume builder
  auto vbuilder = det_builder.new_volume(volume_id::e_cuboid);
  const auto vol_idx{
      static_cast<typename detector_t::surface_type::navigation_link>(
          vbuilder->vol_index())};

  // Add material
  auto mv_builder =
      det_builder.template decorate<homogeneous_material_builder<detector_t>>(
          vbuilder->vol_index());

  assert(mv_builder != nullptr);

  typename detector_t::point3_type t{0.f, 0.f, 20.f};
  mv_builder->add_volume_placement(t);

  // Add a sensitive surface
  auto trpz_factory = std::make_unique<trapezoid_factory>();
  // Add material to the surface
  auto mat_sf_factory =
      std::make_shared<homogeneous_material_factory<detector_t>>(
          std::move(trpz_factory));

  mat_sf_factory->push_back({surface_id::e_sensitive,
                             transform3(point3{0.f, 0.f, 1000.f}), vol_idx,
                             std::vector<scalar>{1.f, 3.f, 2.f, 0.25f}});
  mat_sf_factory->add_material(material_id::e_material_slab,
                               {1.f * unit<scalar>::mm, silicon<scalar>()});

  // Add a portal box around the cuboid volume with a min distance of 'env'
  constexpr auto env{0.1f * unit<scalar>::mm};
  auto portal_generator =
      std::make_unique<cuboid_portal_generator<detector_t>>(env);

  // Add homogeneous material to every portal
  auto mat_portal_factory =
      std::make_shared<homogeneous_material_factory<detector_t>>(
          std::move(portal_generator));
  mat_portal_factory->add_material(material_id::e_material_slab,
                                   {2.f * unit<scalar>::mm, silicon<scalar>()});
  mat_portal_factory->add_material(material_id::e_material_slab,
                                   {3.f * unit<scalar>::mm, silicon<scalar>()});
  mat_portal_factory->add_material(material_id::e_material_slab,
                                   {4.f * unit<scalar>::mm, silicon<scalar>()});
  mat_portal_factory->add_material(material_id::e_material_slab,
                                   {5.f * unit<scalar>::mm, silicon<scalar>()});
  mat_portal_factory->add_material(material_id::e_material_slab,
                                   {6.f * unit<scalar>::mm, silicon<scalar>()});
  mat_portal_factory->add_material(material_id::e_material_slab,
                                   {7.f * unit<scalar>::mm, silicon<scalar>()});

  mv_builder->add_surfaces(mat_sf_factory, geo_ctx);
  mv_builder->add_surfaces(mat_portal_factory);

  //
  // build the detector
  //
  vecmem::host_memory_resource host_mr;
  const detector_t d = det_builder.build(host_mr);
  const auto vol = tracking_volume{d, 0u};

  // check the results
  EXPECT_EQ(d.volumes().size(), 1u);
  EXPECT_EQ(vol.id(), volume_id::e_cuboid);
  EXPECT_EQ(vol.index(), 0u);

  // Check the volume placement
  typename detector_t::transform3_type trf{t};
  EXPECT_TRUE(vol.transform() == trf);
  EXPECT_TRUE(d.transform_store().at(0u) == trf);

  EXPECT_EQ(d.surfaces().size(), 7u);
  EXPECT_EQ(d.mask_store().template size<mask_id::e_rectangle2D>(), 3u);
  EXPECT_EQ(d.mask_store().template size<mask_id::e_trapezoid2D>(), 1u);
  EXPECT_EQ(d.material_store().template size<material_id::e_material_slab>(),
            7u);

  // Check the material links
  for (const auto [idx, sf_desc] : detray::views::enumerate(d.surfaces())) {
    EXPECT_EQ(sf_desc.material().id(), material_id::e_material_slab);
    EXPECT_EQ(sf_desc.material().index(), idx);
  }

  // Check the material
  scalar thickness{1.f * unit<scalar>::mm};
  for (const auto &slab :
       d.material_store().template get<material_id::e_material_slab>()) {
    EXPECT_EQ(slab.get_material(), silicon<scalar>());
    EXPECT_NEAR(slab.thickness(), thickness, tol);
    thickness += 1.f * unit<scalar>::mm;
  }
}
