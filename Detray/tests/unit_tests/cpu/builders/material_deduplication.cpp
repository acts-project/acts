// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s)
#include "detray/builders/detail/material_deduplication.hpp"

#include "detray/builders/cuboid_portal_generator.hpp"
#include "detray/builders/detector_builder.hpp"
#include "detray/builders/homogeneous_material_builder.hpp"
#include "detray/builders/homogeneous_material_factory.hpp"
#include "detray/builders/material_map_builder.hpp"
#include "detray/builders/material_map_factory.hpp"
#include "detray/builders/surface_factory.hpp"
#include "detray/core/detector.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/material/detail/material_accessor.hpp"
#include "detray/material/predefined_materials.hpp"
#include "detray/utils/consistency_checker.hpp"
#include "detray/utils/ranges.hpp"

// Detray test include(s)
#include "detray/test/framework/types.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <array>
#include <memory>
#include <set>
#include <vector>

using namespace detray;

namespace {

using scalar = test::scalar;
using point2 = test::point2;
using point3 = test::point3;

using metadata_t = test::default_metadata;
using detector_builder_t = detector_builder<metadata_t>;
using detector_t = typename detector_builder_t::detector_type;
using transform3 = dtransform3D<typename detector_t::algebra_type>;
using mat_id = typename detector_t::material::id;
using bin_index_t = axis::multi_bin<2u>;

using rectangle_factory = surface_factory<detector_t, rectangle2D>;
using hom_mat_factory_t = homogeneous_material_factory<detector_t>;
using map_factory_t = material_map_factory<detector_t, bin_index_t>;

constexpr dindex n_volumes{3u};

/// Volume local surface indices
/// @{
// Homogeneous material: identical in all volumes
constexpr dindex shared_slab_sf{0u};
// Homogeneous material: different in every volume
constexpr dindex unique_slab_sf{1u};
// Material map: identical in all volumes
constexpr dindex shared_map_sf{2u};
// Material map: identical in all volumes, but on a smaller surface. The grid
// extent is given by the axis spans, so this surface shares the map as well
constexpr dindex shared_map_small_sf{3u};
// Material map: different in every volume
constexpr dindex unique_map_sf{4u};
/// @}

/// Half length of the material maps in x and y
constexpr scalar map_half_length{10.f * unit<scalar>::mm};

/// Add material map data to the factory @param mat_factory for the surface
/// with volume local index @param sf_index
void add_map_data(map_factory_t &mat_factory, std::size_t sf_index, scalar t,
                  const material<scalar> &mat) {
  typename map_factory_t::data_type mat_data{sf_index};
  std::vector<bin_index_t> bins{};
  std::vector<std::size_t> n_bins{5u, 10u};
  std::vector<std::vector<scalar>> axis_spans = {
      {-map_half_length, map_half_length}, {-map_half_length, map_half_length}};

  // Add material for every bin
  for (auto [i, j] : detray::views::cartesian_product{
           detray::views::iota{0u, 5u}, detray::views::iota{0u, 10u}}) {
    bins.push_back({i, j});
    mat_data.append(t, mat);
    t += 0.25f * unit<scalar>::mm;
  }

  mat_factory.add_material(mat_id::e_rectangle2D_map, std::move(mat_data),
                           std::move(n_bins), std::move(axis_spans),
                           std::move(bins));
}

/// Build a detector with @c n_volumes cuboid volumes, in which some surfaces
/// carry identical material
detector_t build_detector(vecmem::memory_resource &mr, const bool dedup) {
  detector_builder_t det_builder{};
  det_builder.set_name("material_deduplication");

  // Off by default
  EXPECT_FALSE(det_builder.deduplicate_material());
  det_builder.deduplicate_material(dedup);
  EXPECT_EQ(det_builder.deduplicate_material(), dedup);

  for (dindex v = 0u; v < n_volumes; ++v) {
    auto vbuilder = det_builder.new_volume(volume_id::e_cuboid);
    vbuilder->add_volume_placement(
        point3{0.f, 0.f, static_cast<scalar>(v) * 100.f * unit<scalar>::mm});
    const auto vol_idx{vbuilder->vol_index()};

    det_builder.template decorate<homogeneous_material_builder<detector_t>>(
        vol_idx);
    auto mat_builder =
        det_builder.template decorate<material_map_builder<detector_t>>(
            vol_idx);

    const auto vol_link{
        static_cast<typename detector_t::surface_type::navigation_link>(
            vol_idx)};
    const scalar vs{static_cast<scalar>(v)};

    // Surfaces with homogeneous material
    auto hom_factory = std::make_shared<hom_mat_factory_t>(
        std::make_unique<rectangle_factory>());

    hom_factory->push_back({surface_id::e_sensitive,
                            transform3(point3{0.f, 0.f, -20.f}), vol_link,
                            std::vector<scalar>{10.f, 10.f}});
    hom_factory->add_material(mat_id::e_material_slab,
                              {1.f * unit<scalar>::mm, silicon<scalar>()});

    hom_factory->push_back({surface_id::e_sensitive,
                            transform3(point3{0.f, 0.f, -10.f}), vol_link,
                            std::vector<scalar>{10.f, 10.f}});
    hom_factory->add_material(
        mat_id::e_material_slab,
        {(1.f + vs) * unit<scalar>::mm, tungsten<scalar>()});

    // Surfaces with material maps
    auto map_factory =
        std::make_shared<map_factory_t>(std::make_unique<rectangle_factory>());

    map_factory->push_back({surface_id::e_sensitive,
                            transform3(point3{0.f, 0.f, 0.f}), vol_link,
                            std::vector<scalar>{10.f, 10.f}});
    add_map_data(*map_factory, shared_map_sf, 1.f * unit<scalar>::mm,
                 silicon<scalar>());

    map_factory->push_back({surface_id::e_sensitive,
                            transform3(point3{0.f, 0.f, 10.f}), vol_link,
                            std::vector<scalar>{6.f, 4.f}});
    add_map_data(*map_factory, shared_map_small_sf, 1.f * unit<scalar>::mm,
                 silicon<scalar>());

    map_factory->push_back({surface_id::e_sensitive,
                            transform3(point3{0.f, 0.f, 20.f}), vol_link,
                            std::vector<scalar>{10.f, 10.f}});
    add_map_data(*map_factory, unique_map_sf, (1.f + vs) * unit<scalar>::mm,
                 gold<scalar>());

    // Add a portal box around the volume
    auto portal_generator =
        std::make_shared<cuboid_portal_generator<detector_t>>(0.1f *
                                                              unit<scalar>::mm);

    mat_builder->add_surfaces(hom_factory);
    mat_builder->add_surfaces(map_factory);
    mat_builder->add_surfaces(portal_generator);
  }

  return det_builder.build(mr);
}

/// @returns the material slab of the surface material at a local point
struct get_material_slab {
  template <typename mat_coll_t, concepts::index index_t>
  auto operator()(const mat_coll_t &mat_coll, const index_t idx,
                  const point2 &loc_p) const -> material_slab<scalar> {
    using material_t = typename mat_coll_t::value_type;

    if constexpr (concepts::material_map<material_t>) {
      return detail::material_accessor::get(
          mat_coll, idx, typename material_t::point_type{loc_p[0], loc_p[1]});
    } else if constexpr (std::same_as<material_t, material_slab<scalar>>) {
      return detail::material_accessor::get(mat_coll, idx, loc_p);
    } else {
      return {};
    }
  }
};

/// @returns the material link of the surface with volume local index
/// @param sf_idx in the volume @param vol_idx
auto material_link(const detector_t &det, dindex vol_idx, dindex sf_idx) {
  // The surfaces with material are the sensitive surfaces, which are added
  // first to every volume
  const auto &vol_desc = det.volumes()[vol_idx];
  const dindex first_sf{
      vol_desc.template sf_link<surface_id::e_sensitive>()[0]};

  const auto &sf_desc = det.surfaces()[first_sf + sf_idx];
  EXPECT_EQ(sf_desc.volume(), vol_idx);

  return sf_desc.material();
}

}  // anonymous namespace

/// Test the content comparison of material grids
GTEST_TEST(detray_builders, material_deduplication_grid_comparison) {
  vecmem::host_memory_resource host_mr;
  const detector_t det = build_detector(host_mr, false);

  const auto &maps =
      det.material_store().template get<mat_id::e_rectangle2D_map>();
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

  // Material slabs
  const material_slab<scalar> slab{silicon<scalar>(), 1.f * unit<scalar>::mm};
  EXPECT_TRUE(detail::is_identical(slab, slab));
  EXPECT_FALSE(detail::is_identical(
      slab, material_slab<scalar>{silicon<scalar>(), 2.f * unit<scalar>::mm}));
  EXPECT_FALSE(detail::is_identical(
      slab, material_slab<scalar>{tungsten<scalar>(), 1.f * unit<scalar>::mm}));
  EXPECT_EQ(detail::material_hash(slab),
            detail::material_hash(material_slab<scalar>{
                silicon<scalar>(), 1.f * unit<scalar>::mm}));
}

/// Build a detector with identical material on surfaces in different volumes
/// with and without material deduplication
GTEST_TEST(detray_builders, material_deduplication) {
  vecmem::host_memory_resource host_mr;

  const detector_t ref_det = build_detector(host_mr, false);
  const detector_t det = build_detector(host_mr, true);

  EXPECT_TRUE(detail::check_consistency(ref_det));
  EXPECT_TRUE(detail::check_consistency(det));

  ASSERT_EQ(ref_det.volumes().size(), n_volumes);
  ASSERT_EQ(det.volumes().size(), n_volumes);
  ASSERT_EQ(ref_det.surfaces().size(), det.surfaces().size());

  // Without deduplication: One material entry per surface
  const auto &ref_store = ref_det.material_store();
  EXPECT_EQ(ref_store.template size<mat_id::e_material_slab>(), 2u * n_volumes);
  EXPECT_EQ(ref_store.template size<mat_id::e_rectangle2D_map>(),
            3u * n_volumes);

  // With deduplication: The shared material is only present once
  const auto &store = det.material_store();
  EXPECT_EQ(store.template size<mat_id::e_material_slab>(), 1u + n_volumes);
  EXPECT_EQ(store.template size<mat_id::e_rectangle2D_map>(), 1u + n_volumes);

  // Check the material links
  std::set<dindex> ref_slab_links;
  std::set<dindex> ref_map_links;
  std::set<dindex> unique_slab_links;
  std::set<dindex> unique_map_links;
  for (dindex v = 0u; v < n_volumes; ++v) {
    // Every surface has its own material entry
    for (dindex sf_idx : {shared_slab_sf, unique_slab_sf}) {
      const auto link = material_link(ref_det, v, sf_idx);
      EXPECT_EQ(link.id(), mat_id::e_material_slab);
      EXPECT_TRUE(ref_slab_links.insert(link.index()).second);
    }
    for (dindex sf_idx : {shared_map_sf, shared_map_small_sf, unique_map_sf}) {
      const auto link = material_link(ref_det, v, sf_idx);
      EXPECT_EQ(link.id(), mat_id::e_rectangle2D_map);
      EXPECT_TRUE(ref_map_links.insert(link.index()).second);
    }

    // The shared material is the first entry that was added
    auto link = material_link(det, v, shared_slab_sf);
    EXPECT_EQ(link.id(), mat_id::e_material_slab);
    EXPECT_EQ(link.index(), 0u);

    link = material_link(det, v, unique_slab_sf);
    EXPECT_EQ(link.id(), mat_id::e_material_slab);
    EXPECT_TRUE(unique_slab_links.insert(link.index()).second);

    for (dindex sf_idx : {shared_map_sf, shared_map_small_sf}) {
      link = material_link(det, v, sf_idx);
      EXPECT_EQ(link.id(), mat_id::e_rectangle2D_map);
      EXPECT_EQ(link.index(), 0u);
    }

    link = material_link(det, v, unique_map_sf);
    EXPECT_EQ(link.id(), mat_id::e_rectangle2D_map);
    EXPECT_TRUE(unique_map_links.insert(link.index()).second);
  }
  EXPECT_FALSE(unique_slab_links.contains(0u));
  EXPECT_FALSE(unique_map_links.contains(0u));

  // The material lookups are unchanged on every surface
  const std::array<point2, 5> loc_points{point2{0.f, 0.f}, point2{-9.f, -9.f},
                                         point2{3.5f, -2.f}, point2{9.f, 9.f},
                                         point2{-5.5f, 7.f}};

  std::size_t n_checked{0u};
  for (const auto &[i, ref_sf_desc] :
       detray::views::enumerate(ref_det.surfaces())) {
    const auto &sf_desc = det.surfaces()[static_cast<dindex>(i)];

    ASSERT_EQ(ref_sf_desc.has_material(), sf_desc.has_material());
    if (!sf_desc.has_material()) {
      continue;
    }
    ASSERT_EQ(ref_sf_desc.material().id(), sf_desc.material().id());

    for (const auto &p : loc_points) {
      const auto ref_slab = ref_store.template visit<get_material_slab>(
          ref_sf_desc.material(), p);
      const auto slab =
          store.template visit<get_material_slab>(sf_desc.material(), p);

      EXPECT_TRUE(detail::is_identical(ref_slab, slab))
          << "surface " << i << ", point (" << p[0] << ", " << p[1] << ")";
      ++n_checked;
    }
  }
  EXPECT_EQ(n_checked, 5u * loc_points.size() * n_volumes);
}
