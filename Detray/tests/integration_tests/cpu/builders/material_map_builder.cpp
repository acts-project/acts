// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Detray include(s)
#include "detray/builders/material_map_builder.hpp"

#include "detray/builders/cuboid_portal_generator.hpp"
#include "detray/builders/detector_builder.hpp"
#include "detray/builders/material_map_factory.hpp"
#include "detray/builders/surface_factory.hpp"
#include "detray/builders/volume_builder.hpp"
#include "detray/core/detector.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/utils/consistency_checker.hpp"
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
using point2 = test::point2;
using point3 = test::point3;
using transform3 = test::transform3;

using metadata_t = test::default_metadata;
using detector_t = host::detector<metadata_t>;
using mat_id = typename detector_t::material::id;
using bin_index_t = axis::multi_bin<2u>;

constexpr dindex n_volumes{3u};

/// Volume local surface indices
/// @{
// Material map: identical in all volumes
constexpr dindex shared_map_sf{0u};
// Material map: identical in all volumes, but on a smaller surface. The grid
// extent is given by the axis spans, so this surface shares the map as well
constexpr dindex shared_map_small_sf{1u};
// Material map: different in every volume
constexpr dindex unique_map_sf{2u};
/// @}

/// Add generate input material for material maps
template <typename material_factory_t, concepts::scalar scalar_t>
auto add_material_data(const material_factory_t &mat_factory, mat_id id,
                       std::size_t sf_index, scalar_t t,
                       material<scalar_t> mat = silicon<scalar_t>()) {
  typename material_factory_t::element_type::data_type mat_data{sf_index};
  std::vector<bin_index_t> bins{};
  std::vector<std::size_t> n_bins{5u, 10u};
  std::vector<std::vector<scalar_t>> axis_spans = {};

  // Add material for every bin
  for (auto [i, j] : detray::views::cartesian_product{
           detray::views::iota{0u, 5u}, detray::views::iota{0u, 10u}}) {
    bins.push_back({i, j});
    mat_data.append(t, mat);
    t += 0.25f * unit<scalar_t>::mm;
  }

  mat_factory->add_material(id, std::move(mat_data), std::move(n_bins),
                            std::move(axis_spans), std::move(bins));
}

/// Build a detector with @c n_volumes cuboid volumes, in which some surfaces
/// carry identical material
detector_t build_detector(vecmem::memory_resource &mr,
                          const volume_builder_options &builder_opts) {
  detector_builder<metadata_t> det_builder{};
  det_builder.set_name("test_detector");

  // Add material maps on rectangle surfaces to a volume builder
  using rectangle_factory = surface_factory<detector_t, rectangle2D>;
  using map_factory_t = material_map_factory<detector_t, bin_index_t>;

  auto map_factory =
      std::make_shared<map_factory_t>(std::make_unique<rectangle_factory>());

  // Add a portal box around the volumes
  auto portal_generator = std::make_shared<cuboid_portal_generator<detector_t>>(
      0.1f * unit<scalar>::mm);

  for (dindex v_idx = 0u; v_idx < n_volumes; ++v_idx) {
    auto vbuilder = det_builder.new_volume(volume_id::e_cuboid);
    vbuilder->add_volume_placement(point3{
        0.f, 0.f, static_cast<scalar>(v_idx) * 100.f * unit<scalar>::mm});

    auto mat_builder =
        det_builder.template decorate<material_map_builder<detector_t>>(v_idx);

    // Surfaces with material maps
    const auto vol_link{
        static_cast<typename detector_t::surface_type::navigation_link>(v_idx)};

    map_factory->push_back({surface_id::e_sensitive,
                            transform3(point3{0.f, 0.f, 0.f}), vol_link,
                            std::vector<scalar>{10.f, 10.f}});
    add_material_data(map_factory, mat_id::e_rectangle2D_map, shared_map_sf,
                      1.f * unit<scalar>::mm, silicon<scalar>());

    // Different surface extent than the first surface
    map_factory->push_back({surface_id::e_sensitive,
                            transform3(point3{0.f, 0.f, 10.f}), vol_link,
                            std::vector<scalar>{6.f, 4.f}});
    add_material_data(map_factory, mat_id::e_rectangle2D_map,
                      shared_map_small_sf, 1.f * unit<scalar>::mm,
                      silicon<scalar>());

    map_factory->push_back({surface_id::e_sensitive,
                            transform3(point3{0.f, 0.f, 20.f}), vol_link,
                            std::vector<scalar>{10.f, 10.f}});
    // Scale material thickness with volume index to make it unique
    add_material_data(map_factory, mat_id::e_rectangle2D_map, unique_map_sf,
                      (1.f + static_cast<scalar>(v_idx)) * unit<scalar>::mm,
                      gold<scalar>());

    mat_builder->add_surfaces(map_factory);
    mat_builder->add_surfaces(portal_generator);

    // Reset for next volume
    map_factory->clear();
    portal_generator->clear();
  }

  return det_builder.build(mr, builder_opts);
}

/// @returns the material slab of the surface material at a local point
struct get_material_slab {
  template <typename mat_coll_t, concepts::index index_t>
  auto operator()(const mat_coll_t &mat_coll, const index_t idx,
                  const point2 &loc_p) const -> material_slab<scalar> {
    using material_t = typename mat_coll_t::value_type;

    // The test surfaces only carry 2D material maps
    if constexpr (concepts::material_map<material_t>) {
      if constexpr (material_t::dim == 2) {
        return detail::material_accessor::get(mat_coll, idx, loc_p);
      }
    }
    return {};
  }
};

}  // anonymous namespace

/// Integration test: material builder as volume builder decorator
GTEST_TEST(detray_builders, decorator_material_map_builder) {
  using test_algebra = typename detector_t::algebra_type;
  using scalar = dscalar<test_algebra>;
  using transform3 = dtransform3D<test_algebra>;
  using mask_id = typename detector_t::masks::id;

  using pt_cylinder_factory_t =
      surface_factory<detector_t, concentric_cylinder2D>;
  using rectangle_factory = surface_factory<detector_t, rectangle2D>;
  using trapezoid_factory = surface_factory<detector_t, trapezoid2D>;
  using cylinder_factory = surface_factory<detector_t, cylinder2D>;

  using mat_factory_t = material_map_factory<detector_t, bin_index_t>;

  vecmem::host_memory_resource host_mr;
  detector_t d(host_mr);
  auto geo_ctx = typename detector_t::geometry_context{};
  volume_builder_options builder_opts{};
  builder_opts.deduplicate(false);

  auto vbuilder =
      std::make_unique<volume_builder<detector_t>>(volume_id::e_cylinder);
  auto mat_builder = material_map_builder<detector_t>{std::move(vbuilder)};

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
  scalar t{1.f * unit<scalar>::mm};
  add_material_data(mat_pt_cyl_factory, mat_id::e_concentric_cylinder2D_map, 0u,
                    t, silicon<scalar>());
  add_material_data(mat_pt_cyl_factory, mat_id::e_concentric_cylinder2D_map, 1u,
                    t, silicon<scalar>());

  auto mat_rect_factory =
      std::make_shared<mat_factory_t>(std::move(rect_factory));
  t = 1.f * unit<scalar>::mm;
  add_material_data(mat_rect_factory, mat_id::e_rectangle2D_map, 2u, t,
                    tungsten<scalar>());
  // No material for surface with index 3
  t = 3.f * unit<scalar>::mm;
  add_material_data(mat_rect_factory, mat_id::e_rectangle2D_map, 4u, t,
                    tungsten<scalar>());

  auto mat_trpz_factory =
      std::make_shared<mat_factory_t>(std::move(trpz_factory));
  t = 1.f * unit<scalar>::mm;
  add_material_data(mat_trpz_factory, mat_id::e_trapezoid2D_map, 5u, t,
                    tungsten<scalar>());

  auto mat_cyl_factory =
      std::make_shared<mat_factory_t>(std::move(cyl_factory));
  t = 1.5f * unit<scalar>::mm;
  add_material_data(mat_cyl_factory, mat_id::e_cylinder2D_map, 6u, t,
                    gold<scalar>());

  // Add surfaces and material to detector
  mat_builder.add_surfaces(mat_pt_cyl_factory, geo_ctx);
  mat_builder.add_surfaces(mat_rect_factory, geo_ctx);
  mat_builder.add_surfaces(mat_trpz_factory, geo_ctx);
  mat_builder.add_surfaces(mat_cyl_factory, geo_ctx);

  // Add the volume to the detector
  mat_builder.build(d, builder_opts);

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

  EXPECT_EQ(d.material_store().template size<mat_id::e_material_slab>(), 0u);
  EXPECT_EQ(d.material_store().template size<mat_id::e_material_rod>(), 0u);
  EXPECT_EQ(d.material_store().template size<mat_id::e_ring2D_map>(), 0u);
  EXPECT_EQ(d.material_store().template size<mat_id::e_annulus2D_map>(), 0u);
  EXPECT_EQ(d.material_store().template size<mat_id::e_cylinder2D_map>(), 1u);
  EXPECT_EQ(
      d.material_store().template size<mat_id::e_concentric_cylinder2D_map>(),
      2u);
  // Rectangle and trapezoid surfaces have the same grid geometry
  EXPECT_EQ(d.material_store().template size<mat_id::e_rectangle2D_map>(), 3u);
  EXPECT_EQ(d.material_store().template size<mat_id::e_trapezoid2D_map>(), 3u);

  // Check the material links
  std::size_t pt_cyl_idx{0u};
  std::size_t cyl_idx{0u};
  std::size_t cart_idx{0u};
  for (auto &sf_desc : d.surfaces()) {
    const auto &mat_link = sf_desc.material();
    switch (mat_link.id()) {
      case mat_id::e_cylinder2D_map: {
        EXPECT_EQ(mat_link.index(), cyl_idx++) << sf_desc;
        break;
      }
      case mat_id::e_concentric_cylinder2D_map: {
        EXPECT_EQ(mat_link.index(), pt_cyl_idx++) << sf_desc;
        break;
      }
      case mat_id::e_rectangle2D_map: {
        EXPECT_EQ(mat_link.index(), cart_idx++) << sf_desc;
        break;
      }
      case mat_id::e_none: {
        // No material on surface 3
        EXPECT_TRUE(sf_desc.index() == 3u) << "No material on: " << sf_desc;
        EXPECT_TRUE(mat_link.id() == mat_id::e_none) << sf_desc;
        EXPECT_TRUE(mat_link.is_invalid()) << sf_desc;
        break;
      }
      default: {
        EXPECT_TRUE(false) << sf_desc;
      }
    }
  }

  // Check the material map content
  for (auto cyl_mat_grid :
       d.material_store().template get<mat_id::e_concentric_cylinder2D_map>()) {
    EXPECT_EQ(cyl_mat_grid.nbins(), 50u);
    EXPECT_EQ(cyl_mat_grid.size(), 50u);

    auto r_axis = cyl_mat_grid.template get_axis<axis::label::e_rphi>();
    EXPECT_EQ(r_axis.nbins(), 5u);
    auto z_axis = cyl_mat_grid.template get_axis<axis::label::e_cyl_z>();
    EXPECT_EQ(z_axis.nbins(), 10u);

    for (const auto &mat_slab : cyl_mat_grid.all()) {
      EXPECT_TRUE(mat_slab.get_material() == silicon<scalar>() ||
                  mat_slab.get_material() == gold<scalar>());
    }
  }

  for (auto cart_mat_grid :
       d.material_store().template get<mat_id::e_rectangle2D_map>()) {
    EXPECT_EQ(cart_mat_grid.nbins(), 50u);
    EXPECT_EQ(cart_mat_grid.size(), 50u);

    auto x_axis = cart_mat_grid.template get_axis<axis::label::e_x>();
    EXPECT_EQ(x_axis.nbins(), 5u);
    auto y_axis = cart_mat_grid.template get_axis<axis::label::e_y>();
    EXPECT_EQ(y_axis.nbins(), 10u);

    for (const auto &mat_slab : cart_mat_grid.all()) {
      EXPECT_TRUE(mat_slab.get_material() == tungsten<scalar>());
    }
  }
}

/// Build a detector with identical material on surfaces in different volumes
/// with and without material deduplication
GTEST_TEST(detray_builders, material_map_deduplication) {
  vecmem::host_memory_resource host_mr;
  volume_builder_options builder_opts{};

  builder_opts.deduplicate(false);
  const detector_t ref_det = build_detector(host_mr, builder_opts);
  builder_opts.deduplicate(true);
  const detector_t det = build_detector(host_mr, builder_opts);

  EXPECT_TRUE(detail::check_consistency(ref_det));
  EXPECT_TRUE(detail::check_consistency(det));

  ASSERT_EQ(ref_det.volumes().size(), n_volumes);
  ASSERT_EQ(det.volumes().size(), n_volumes);
  ASSERT_EQ(ref_det.surfaces().size(), det.surfaces().size());

  // Test the material map content

  // Without deduplication: One material entry per surface
  const auto &ref_store = ref_det.material_store();
  const auto &maps = ref_store.template get<mat_id::e_rectangle2D_map>();
  ASSERT_EQ(maps.size(), 3u * n_volumes);

  // The map content is identical for the first two map surfaces
  EXPECT_TRUE(detail::is_identical_grid(maps[0], maps[0]));
  EXPECT_TRUE(detail::is_identical_grid(maps[0], maps[1]));
  EXPECT_TRUE(detail::is_identical_grid(maps[0], maps[3]));
  EXPECT_FALSE(detail::is_identical_grid(maps[0], maps[2]));
  // Unique map in the first and second volume
  EXPECT_FALSE(detail::is_identical_grid(maps[2], maps[5]));

  // In contrast, the grid operator== compares the storage position of the
  // non-owning grids
  EXPECT_FALSE(maps[0] == maps[1]);

  // With deduplication: The shared material is only present once
  const auto &store = det.material_store();
  EXPECT_EQ(store.template size<mat_id::e_rectangle2D_map>(), 1u + n_volumes);

  // Check the material links
  std::set<dindex> ref_map_links;
  std::set<dindex> unique_map_links;
  for (dindex v = 0u; v < n_volumes; ++v) {
    const auto vol_desc = det.volume(v);
    ASSERT_EQ(vol_desc, ref_det.volume(v));

    // Every surface has its own material entry in the ref detector
    for (dindex sf_idx : {shared_map_sf, shared_map_small_sf, unique_map_sf}) {
      const auto sf_desc = ref_det.surface(vol_desc.to_global_sf_index(sf_idx));
      const auto link = sf_desc.material();
      EXPECT_EQ(link.id(), mat_id::e_rectangle2D_map);
      // Insertion is successful if the index is uniuqe
      EXPECT_TRUE(ref_map_links.insert(link.index()).second);
    }

    // Shared material links to the first map
    for (dindex sf_idx : {shared_map_sf, shared_map_small_sf}) {
      const auto sf_desc = det.surface(vol_desc.to_global_sf_index(sf_idx));
      const auto link = sf_desc.material();
      EXPECT_EQ(link.id(), mat_id::e_rectangle2D_map);
      EXPECT_EQ(link.index(), 0u);
    }

    const auto sf_desc =
        det.surface(vol_desc.to_global_sf_index(unique_map_sf));
    const auto link = sf_desc.material();
    EXPECT_EQ(link.id(), mat_id::e_rectangle2D_map);
    EXPECT_TRUE(unique_map_links.insert(link.index()).second);
  }
  EXPECT_FALSE(unique_map_links.contains(0u));

  // The material lookups are unchanged on every surface
  const std::array<point2, 5> loc_points{point2{0.f, 0.f}, point2{-9.f, -9.f},
                                         point2{3.5f, -2.f}, point2{9.f, 9.f},
                                         point2{-5.5f, 7.f}};

  std::size_t n_checked{0u};
  for (const auto ref_sf_desc : ref_det.surfaces()) {
    const geometry::surface ref_sf{ref_det, ref_sf_desc};
    const auto sf_desc = det.surface(ref_sf.index());
    const geometry::surface sf{det, sf_desc};

    // Check material link
    ASSERT_EQ(ref_sf.has_material(), sf.has_material());
    if (!sf.has_material()) {
      continue;
    }
    ASSERT_EQ(ref_sf_desc.material().id(), sf_desc.material().id());

    // Check lookup
    for (const auto &p : loc_points) {
      const auto ref_slab =
          ref_sf.template visit_material<get_material_slab>(p);
      const auto slab = sf.template visit_material<get_material_slab>(p);

      EXPECT_EQ(ref_slab, slab) << "surface " << ref_sf.index() << ", point ("
                                << p[0] << ", " << p[1] << ")";
      ++n_checked;
    }
  }
  // Have 3 material surfaces per volume
  EXPECT_EQ(n_checked, 3u * loc_points.size() * n_volumes);
}
