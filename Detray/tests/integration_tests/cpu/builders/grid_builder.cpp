// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Detray include(s)
#include "detray/builders/grid_builder.hpp"

#include "detray/builders/surface_factory.hpp"
#include "detray/builders/volume_builder.hpp"
#include "detray/core/detector.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/utils/type_list.hpp"

// Detray test include(s)
#include "detray/test/framework/types.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// Gtest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <limits>

using namespace detray;
using namespace detray::axis;

namespace {

using scalar = test::scalar;
using point3 = test::point3;
using vector3 = test::vector3;

using detector_t = detector<test::toy_metadata>;

}  // anonymous namespace

/// Integration test: grid builder as volume builder decorator
GTEST_TEST(detray_builders, decorator_grid_builder) {
  using algebra_t = typename detector_t::algebra_type;
  using transform3 = dtransform3D<algebra_t>;
  using geo_obj_id = typename detector_t::geo_obj_ids;
  using acc_ids = typename detector_t::accel::id;
  using mask_id = typename detector_t::masks::id;

  // cylinder grid type of the toy detector
  using cyl_grid_t = grid<algebra_t, axes<concentric_cylinder2D>,
                          bins::static_array<detector_t::surface_type, 1>,
                          simple_serializer, host_container_types, false>;

  using pt_cylinder_t = concentric_cylinder2D;
  using pt_cylinder_factory_t = surface_factory<detector_t, pt_cylinder_t>;
  using rectangle_factory = surface_factory<detector_t, rectangle2D>;
  using trapezoid_factory = surface_factory<detector_t, trapezoid2D>;
  using cylinder_factory = surface_factory<detector_t, pt_cylinder_t>;

  vecmem::host_memory_resource host_mr;
  detector_t d(host_mr);
  auto geo_ctx = typename detector_t::geometry_context{};
  const auto vol_idx{
      static_cast<typename detector_t::surface_type::navigation_link>(
          d.volumes().size())};

  auto vbuilder = std::make_unique<volume_builder<detector_t>>(
      volume_id::e_cylinder, vol_idx);
  auto gbuilder = grid_builder<detector_t, cyl_grid_t>{std::move(vbuilder)};
  // passive surfaces are added to the grid
  // gbuilder.set_add_surfaces();

  // The cylinder portals are at the end of the surface range by construction
  const auto cyl_mask =
      mask<concentric_cylinder2D, algebra_t>{0u, 10.f, -500.f, 500.f};
  std::size_t n_phi_bins{5u};
  std::size_t n_z_bins{4u};

  // Build empty grid
  gbuilder.init_grid(cyl_mask, {n_phi_bins, n_z_bins});

  EXPECT_TRUE(d.volumes().size() == 0);

  // Add some portals first
  auto pt_cyl_factory = std::make_shared<pt_cylinder_factory_t>();

  typename pt_cylinder_factory_t::sf_data_collection cyl_sf_data;
  cyl_sf_data.emplace_back(surface_id::e_portal,
                           transform3(point3{0.f, 0.f, 0.f}), 0u,
                           std::vector<scalar>{10.f, -1500.f, 1500.f});
  cyl_sf_data.emplace_back(surface_id::e_portal,
                           transform3(point3{0.f, 0.f, 0.f}), 2u,
                           std::vector<scalar>{20.f, -1500.f, 1500.f});
  pt_cyl_factory->push_back(std::move(cyl_sf_data));

  // Then some passive and sensitive surfaces
  auto rect_factory = std::make_shared<rectangle_factory>();

  typename rectangle_factory::sf_data_collection rect_sf_data;
  rect_sf_data.emplace_back(surface_id::e_sensitive,
                            transform3(point3{7.07f, 7.07f, -500.f}), vol_idx,
                            std::vector<scalar>{10.f, 8.f});
  rect_sf_data.emplace_back(surface_id::e_sensitive,
                            transform3(point3{7.07f, 7.07f, -250.f}), vol_idx,
                            std::vector<scalar>{10.f, 8.f});
  rect_sf_data.emplace_back(surface_id::e_sensitive,
                            transform3(point3{7.07f, 7.07f, 100.f}), vol_idx,
                            std::vector<scalar>{10.f, 8.f});
  rect_factory->push_back(std::move(rect_sf_data));

  auto trpz_factory = std::make_shared<trapezoid_factory>();

  typename trapezoid_factory::sf_data_collection trpz_sf_data;
  trpz_sf_data.emplace_back(surface_id::e_sensitive,
                            transform3(point3{7.07f, 7.07f, 600.f}), vol_idx,
                            std::vector<scalar>{1.f, 3.f, 2.f, 0.25f});
  trpz_factory->push_back(std::move(trpz_sf_data));

  auto cyl_factory = std::make_shared<cylinder_factory>();

  cyl_sf_data.clear();
  cyl_sf_data.emplace_back(surface_id::e_passive,
                           transform3(point3{0.f, 0.f, 0.f}), vol_idx,
                           std::vector<scalar>{5.f, -1300.f, 1300.f});
  cyl_factory->push_back(std::move(cyl_sf_data));

  // senstivies should end up in the grid, portals in the volume
  gbuilder.add_surfaces(pt_cyl_factory, geo_ctx);
  gbuilder.add_surfaces(rect_factory, geo_ctx);
  gbuilder.add_surfaces(trpz_factory, geo_ctx);
  gbuilder.add_surfaces(cyl_factory, geo_ctx);

  const auto& cyl_axis_z = gbuilder.get().template get_axis<label::e_cyl_z>();
  EXPECT_EQ(cyl_axis_z.label(), label::e_cyl_z);
  EXPECT_EQ(cyl_axis_z.nbins(), 4u);

  gbuilder.build(d);

  //
  // check results
  //
  const auto& vol = d.volumes().back();
  EXPECT_TRUE(d.volumes().size() == 1u);
  EXPECT_EQ(vol.index(), 0u);
  EXPECT_EQ(vol.id(), volume_id::e_cylinder);

  // only the portals are referenced through the volume
  test::toy_metadata::object_link_type sf_range{};
  sf_range[0] = {acc_ids::e_surface_default, 0u};
  sf_range[1] = {acc_ids::e_surface_concentric_cylinder2D_grid, 0u};

  // toy detector makes no distinction between the surface types
  EXPECT_EQ(vol.template accel_link<geo_obj_id::e_portal>(),
            sf_range[geo_obj_id::e_portal]);
  EXPECT_EQ(vol.template accel_link<geo_obj_id::e_sensitive>(),
            sf_range[geo_obj_id::e_sensitive]);
  EXPECT_EQ(vol.template accel_link<geo_obj_id::e_passive>(),
            sf_range[geo_obj_id::e_passive]);

  // Only the portals should be in the detector's surface container now
  EXPECT_EQ(d.surfaces().size(), 7u);
  EXPECT_EQ(d.mask_store().template size<mask_id::e_concentric_cylinder2D>(),
            3u);
  EXPECT_EQ(d.mask_store().template size<mask_id::e_ring2D>(), 0u);
  EXPECT_EQ(d.mask_store().template size<mask_id::e_rectangle2D>(), 3u);
  EXPECT_EQ(d.mask_store().template size<mask_id::e_trapezoid2D>(), 1u);

  // check the portals in the detector
  const auto& bf_finder =
      d.accelerator_store()
          .template get<detector_t::accel::id::e_surface_brute_force>()[0];
  for (const auto& sf : bf_finder.all()) {
    EXPECT_TRUE((sf.id() == surface_id::e_portal) ||
                (sf.id() == surface_id::e_passive));
    EXPECT_EQ(sf.volume(), 0u);
  }

  // check the sensitive surfaces in the grid
  const auto& cyl_grid =
      d.accelerator_store()
          .template get<
              detector_t::accel::id::e_surface_concentric_cylinder2D_grid>()[0];
  dindex trf_idx{3u};
  for (const auto& sf : cyl_grid.all()) {
    EXPECT_TRUE(sf.is_sensitive());
    EXPECT_EQ(sf.volume(), 0u);
    EXPECT_EQ(sf.transform(), trf_idx++);
  }
}
