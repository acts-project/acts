// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s)
#include "detray/builders/volume_builder.hpp"

#include "detray/builders/cuboid_portal_generator.hpp"
#include "detray/builders/surface_factory.hpp"
#include "detray/core/concepts.hpp"
#include "detray/core/detector.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes.hpp"

// Detray test include(s)
#include "detray/test/framework/types.hpp"
#include "detray/test/utils/prefill_detector.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

namespace {

using scalar = detray::test::scalar;
using point3 = detray::test::point3;

/// Check volume links for a collection of masks in a given detector
template <detray::concepts::detector detector_t,
          typename detector_t::masks::id mask_id>
inline void check_mask(const detector_t& d,
                       const std::vector<detray::dindex>& vol_links) {
  for (const auto [idx, mask] :
       detray::views::enumerate(d.mask_store().template get<mask_id>())) {
    EXPECT_EQ(mask.volume_link(), vol_links.at(idx))
        << "mask no. " << idx << ": " << mask.to_string();
  }
}

/// Build a cylinder volume with contained surfaces into the detector @p d
///
/// @param d the detector (prefilled)
/// @param builder_opts the volume builder options (e.g. deduplication)
/// @param t the volume placement
template <detray::concepts::detector detector_t>
void build_test_volume(detector_t& d,
                       const detray::volume_builder_options& builder_opts,
                       const typename detector_t::point3_type& t) {
  using namespace detray;

  using transform3 = typename detector_t::transform3_type;

  // Surface factories
  using portal_cylinder_factory =
      surface_factory<detector_t, concentric_cylinder2D>;
  using annulus_factory = surface_factory<detector_t, annulus2D>;
  using cylinder_factory = surface_factory<detector_t, cylinder2D>;
  using rectangle_factory = surface_factory<detector_t, rectangle2D>;
  using disc_factory = surface_factory<detector_t, ring2D>;
  using trapezoid_factory = surface_factory<detector_t, trapezoid2D>;

  auto geo_ctx = typename detector_t::geometry_context{};
  const auto vol_idx{
      static_cast<typename detector_t::surface_type::navigation_link>(
          d.volumes().size())};

  // volume builder
  volume_builder<detector_t> vbuilder{volume_id::e_cylinder};
  vbuilder.add_volume_placement(t);

  //
  // Fill the surface factories with data
  //

  // portal surfaces
  auto pt_cyl_factory = std::make_shared<portal_cylinder_factory>();
  typename portal_cylinder_factory::sf_data_collection cyl_sf_data;
  // Creates two cylinders at radius 0mm and 10mm with an extent in z
  // of -5mm to 5mm and linking to volumes 0 and 2, respectively.
  cyl_sf_data.emplace_back(surface_id::e_portal,
                           transform3(point3{0.f, 0.f, 0.f}), 0u,
                           std::vector<scalar>{10.f, -1500.f, 1500.f});
  cyl_sf_data.emplace_back(surface_id::e_portal,
                           transform3(point3{0.f, 0.f, 0.f}), 2u,
                           std::vector<scalar>{20.f, -1500.f, 1500.f});
  pt_cyl_factory->push_back(std::move(cyl_sf_data));

  auto pt_disc_factory = std::make_shared<disc_factory>();
  typename disc_factory::sf_data_collection disc_sf_data;
  // Creates two discs with a radius of 10mm, linking to volumes 3 and 4
  disc_sf_data.emplace_back(surface_id::e_portal,
                            transform3(point3{0.f, 0.f, -1500.f}), 3u,
                            std::vector<scalar>{0.f, 10.f});
  disc_sf_data.emplace_back(surface_id::e_portal,
                            transform3(point3{0.f, 0.f, 1500.f}), 4u,
                            std::vector<scalar>{0.f, 10.f});
  pt_disc_factory->push_back(std::move(disc_sf_data));

  // sensitive surfaces
  auto ann_factory = std::make_shared<annulus_factory>();
  typename annulus_factory::sf_data_collection ann_sf_data;
  ann_sf_data.emplace_back(
      surface_id::e_sensitive, transform3(point3{0.f, 0.f, 1000.f}), vol_idx,
      std::vector<scalar>{300.f, 350.f, -0.1f, 0.1f, 0.5f, 0.6f, 1.4f});
  ann_sf_data.emplace_back(
      surface_id::e_sensitive, transform3(point3{0.f, 0.f, 1000.f}), vol_idx,
      std::vector<scalar>{350.f, 400.f, -0.1f, 0.1f, 0.5f, 0.6f, 1.4f});
  ann_factory->push_back(std::move(ann_sf_data));

  auto rect_factory = std::make_shared<rectangle_factory>();
  typename rectangle_factory::sf_data_collection rect_sf_data;
  rect_sf_data.emplace_back(surface_id::e_sensitive,
                            transform3(point3{0.f, 0.f, -10.f}), vol_idx,
                            std::vector<scalar>{-10.f, -8.f, 10.f, 8.f});
  rect_sf_data.emplace_back(surface_id::e_sensitive,
                            transform3(point3{0.f, 0.f, -20.f}), vol_idx,
                            std::vector<scalar>{-10.f, -8.f, 10.f, 8.f});
  rect_sf_data.emplace_back(surface_id::e_sensitive,
                            transform3(point3{0.f, 0.f, -30.f}), vol_idx,
                            std::vector<scalar>{-10.f, -8.f, 10.f, 8.f});
  rect_factory->push_back(std::move(rect_sf_data));

  auto trpz_factory = std::make_shared<trapezoid_factory>();
  typename trapezoid_factory::sf_data_collection trpz_sf_data;
  trpz_sf_data.emplace_back(surface_id::e_sensitive,
                            transform3(point3{0.f, 0.f, 1000.f}), vol_idx,
                            std::vector<scalar>{1.f, 3.f, 2.f, 0.25f});
  trpz_factory->push_back(std::move(trpz_sf_data));

  // passive surfaces
  auto cyl_factory = std::make_shared<cylinder_factory>();
  cyl_sf_data.clear();
  cyl_sf_data.emplace_back(surface_id::e_passive,
                           transform3(point3{0.f, 0.f, 0.f}), vol_idx,
                           std::vector<scalar>{5.f, -1300.f, 1300.f});
  cyl_factory->push_back(std::move(cyl_sf_data));

  auto sf_disc_factory = std::make_shared<disc_factory>();
  disc_sf_data.clear();
  disc_sf_data.emplace_back(surface_id::e_passive,
                            transform3(point3{0.f, 0.f, -1300.f}), vol_idx,
                            std::vector<scalar>{0.f, 5.f});
  disc_sf_data.emplace_back(surface_id::e_passive,
                            transform3(point3{0.f, 0.f, 1300.f}), vol_idx,
                            std::vector<scalar>{0.f, 5.f});
  sf_disc_factory->push_back(std::move(disc_sf_data));

  //
  // Fill everything into volume
  //
  vbuilder.add_surfaces(pt_cyl_factory, geo_ctx);
  vbuilder.add_surfaces(pt_disc_factory, geo_ctx);

  vbuilder.add_surfaces(ann_factory, geo_ctx);
  vbuilder.add_surfaces(rect_factory, geo_ctx);
  vbuilder.add_surfaces(trpz_factory, geo_ctx);

  vbuilder.add_surfaces(cyl_factory, geo_ctx);
  vbuilder.add_surfaces(sf_disc_factory, geo_ctx);

  // try adding something extra later...
  rect_factory->clear();
  rect_sf_data.clear();
  rect_sf_data.emplace_back(surface_id::e_sensitive,
                            transform3(point3{0.f, 0.f, 10.f}), vol_idx,
                            std::vector<scalar>{-10.f, -8.f, 10.f, 8.f});
  rect_sf_data.emplace_back(surface_id::e_sensitive,
                            transform3(point3{0.f, 0.f, 20.f}), vol_idx,
                            std::vector<scalar>{-10.f, -8.f, 10.f, 8.f});
  rect_factory->push_back(std::move(rect_sf_data));
  rect_sf_data.clear();
  rect_sf_data.emplace_back(surface_id::e_sensitive,
                            transform3(point3{0.f, 0.f, 30.f}), vol_idx,
                            std::vector<scalar>{-10.f, -8.f, 10.f, 8.f});
  rect_factory->push_back(std::move(rect_sf_data));

  vbuilder.add_surfaces(rect_factory, geo_ctx);

  //
  // Adds all surfaces to the detector
  //
  vbuilder.build(d, builder_opts);
}

}  // anonymous namespace

/// Integration test to build a cylinder volume with contained surfaces
GTEST_TEST(detray_builders, tracking_volume_construction) {
  using namespace detray;

  using metadata_t = test::default_metadata;
  using detector_t = host::detector<metadata_t>;
  using geo_obj_id = typename detector_t::geo_obj_ids;
  using mask_id = typename detector_t::masks::id;
  using accel_id = typename detector_t::accel::id;

  // detector
  vecmem::host_memory_resource host_mr;
  detector_t d(host_mr);
  auto geo_ctx = typename detector_t::geometry_context{};
  // Check the data layout without deduplication
  const auto builder_opts = volume_builder_options{}.deduplicate(false);

  // ensure there is a data offset that needs to be handled correctly
  prefill_detector(d, geo_ctx);
  const dindex first_trf{d.transform_store().size()};

  // initial checks
  EXPECT_EQ(d.volumes().size(), 1u);
  EXPECT_EQ(d.portals().size(), 3u);

  //
  // Adds all surfaces to the detector
  //
  typename detector_t::point3_type t{0.f, 0.f, 20.f};
  build_test_volume(d, builder_opts, t);

  //
  // check results
  //
  const auto& vol = d.volumes().back();

  EXPECT_EQ(d.volumes().size(), 2u);
  EXPECT_EQ(vol.index(), 1u);
  EXPECT_EQ(vol.id(), volume_id::e_cylinder);

  // Check the volume placement
  typename detector_t::transform3_type trf{t};
  EXPECT_TRUE(d.transform_store().at(first_trf) == trf);

  // Check the acceleration data structure link
  dtyped_index<accel_id, dindex> acc_link{accel_id::e_surface_default, 1u};
  ASSERT_TRUE(vol.accel_link().size() == geo_obj_id::e_size);
  EXPECT_EQ(vol.accel_link<geo_obj_id::e_portal>(), acc_link);
  EXPECT_EQ(vol.accel_link<geo_obj_id::e_passive>(), acc_link);
  // Not set by the vanilla volume builder
  EXPECT_TRUE(
      detail::is_invalid_value(vol.accel_link<geo_obj_id::e_sensitive>()));

  EXPECT_EQ(d.portals().size(), 19u);
  EXPECT_EQ(d.mask_store().template size<mask_id::e_concentric_cylinder2D>(),
            2u);
  EXPECT_EQ(d.mask_store().template size<mask_id::e_ring2D>(), 4u);
  EXPECT_EQ(d.mask_store().template size<mask_id::e_annulus2D>(), 3u);
  EXPECT_EQ(d.mask_store().template size<mask_id::e_cylinder2D>(), 1u);
  EXPECT_EQ(d.mask_store().template size<mask_id::e_rectangle2D>(), 7u);
  EXPECT_EQ(d.mask_store().template size<mask_id::e_ring2D>(), 4u);
  EXPECT_EQ(d.mask_store().template size<mask_id::e_trapezoid2D>(), 2u);

  // check surface type and volume link
  std::vector<surface_id> sf_ids{};
  sf_ids.reserve(d.surfaces().size());
  sf_ids.insert(sf_ids.end(), 3u, surface_id::e_sensitive);
  sf_ids.insert(sf_ids.end(), 4u, surface_id::e_portal);
  sf_ids.insert(sf_ids.end(), 6u, surface_id::e_sensitive);
  sf_ids.insert(sf_ids.end(), 3u, surface_id::e_passive);
  sf_ids.insert(sf_ids.end(), 3u, surface_id::e_sensitive);

  std::vector<dindex> volume_links{};
  volume_links.reserve(d.surfaces().size());
  volume_links.insert(volume_links.end(), 3u, 0u);
  volume_links.insert(volume_links.end(), 16u, 1u);

  // Check surface id and volume links
  for (const auto [idx, sf_id] : detray::views::enumerate(sf_ids)) {
    geometry::identifier geo_id{};
    geo_id.set_index(idx);
    const auto& sf = d.surface(geo_id);
    EXPECT_EQ(sf.id(), sf_id) << "error at index: " << idx;
    EXPECT_EQ(sf.volume(), volume_links.at(idx)) << "error at index: " << idx;
  }

  // check that the transform indices are continuous for the newly added
  // surfaces. The first new transform belongs to the volume itself
  for (std::size_t idx :
       detray::views::iota(dindex_range{3, d.surfaces().size()})) {
    geometry::identifier geo_id{};
    geo_id.set_index(idx);
    // Add a shift to the index for the volume placement transforms
    EXPECT_EQ(d.surface(geo_id).transform(), idx + 2)
        << "error at index: " << idx;
  }

  // check surface mask links
  std::vector<typename detector_t::surface_type::mask_link> mask_links{
      {mask_id::e_rectangle2D, {0u, 1u}},
      {mask_id::e_annulus2D, {0u, 1u}},
      {mask_id::e_trapezoid2D, {0u, 1u}},
      {mask_id::e_concentric_cylinder2D, {0u, 1u}},
      {mask_id::e_concentric_cylinder2D, {1u, 1u}},
      {mask_id::e_ring2D, {0u, 1u}},
      {mask_id::e_ring2D, {1u, 1u}},
      {mask_id::e_annulus2D, {1u, 1u}},
      {mask_id::e_annulus2D, {2u, 1u}},
      {mask_id::e_rectangle2D, {1u, 1u}},
      {mask_id::e_rectangle2D, {2u, 1u}},
      {mask_id::e_rectangle2D, {3u, 1u}},
      {mask_id::e_trapezoid2D, {1u, 1u}},
      {mask_id::e_cylinder2D, {0u, 1u}},
      {mask_id::e_ring2D, {2u, 1u}},
      {mask_id::e_ring2D, {3u, 1u}},
      {mask_id::e_rectangle2D, {4u, 1u}},
      {mask_id::e_rectangle2D, {5u, 1u}},
      {mask_id::e_rectangle2D, {6u, 1u}}};
  for (const auto [idx, m_link] : detray::views::enumerate(mask_links)) {
    geometry::identifier geo_id{};
    geo_id.set_index(idx);
    EXPECT_EQ(d.surface(geo_id).mask(), m_link) << "error at index: " << idx;
  }

  // check mask volume links
  volume_links.clear();
  volume_links = {0u, 2u};
  check_mask<detector_t, mask_id::e_concentric_cylinder2D>(d, volume_links);

  volume_links.clear();
  volume_links = {1u};
  check_mask<detector_t, mask_id::e_cylinder2D>(d, volume_links);

  volume_links.clear();
  volume_links = {3u, 4u, 1u, 1u};
  check_mask<detector_t, mask_id::e_ring2D>(d, volume_links);
  check_mask<detector_t, mask_id::e_ring2D>(d, volume_links);

  volume_links.clear();
  volume_links = {0u, 1u, 1u};
  check_mask<detector_t, mask_id::e_annulus2D>(d, volume_links);

  volume_links.clear();
  volume_links.reserve(7u);
  volume_links.push_back(0u);
  volume_links.insert(volume_links.end(), 6u, 1u);
  check_mask<detector_t, mask_id::e_rectangle2D>(d, volume_links);

  volume_links.clear();
  volume_links = {0u, 1u};
  check_mask<detector_t, mask_id::e_trapezoid2D>(d, volume_links);
}

/// Integration test to build a cylinder volume with contained surfaces and
/// deduplicated transforms and masks
GTEST_TEST(detray_builders, tracking_volume_construction_deduplication) {
  using namespace detray;

  using metadata_t = test::default_metadata;
  using detector_t = host::detector<metadata_t>;
  using mask_id = typename detector_t::masks::id;

  // detector
  vecmem::host_memory_resource host_mr;
  detector_t d(host_mr);
  auto geo_ctx = typename detector_t::geometry_context{};
  // Deduplication is the default
  const volume_builder_options builder_opts{};
  ASSERT_TRUE(builder_opts.deduplicate());

  // ensure there is a data offset that needs to be handled correctly
  prefill_detector(d, geo_ctx);
  const dindex first_trf{d.transform_store().size()};
  EXPECT_EQ(first_trf, 4u);

  typename detector_t::point3_type t{0.f, 0.f, 20.f};
  build_test_volume(d, builder_opts, t);

  //
  // check results
  //
  EXPECT_EQ(d.volumes().size(), 2u);
  EXPECT_EQ(d.surfaces().size(), 19u);
  EXPECT_EQ(d.portals().size(), 19u);

  // The volume placement is new
  const auto& vol = d.volumes().back();
  EXPECT_EQ(vol.transform(), first_trf);
  EXPECT_TRUE(d.transform_store().at(vol.transform()) ==
              typename detector_t::transform3_type{t});

  // Without deduplication: 4 (prefill) + 1 (volume) + 16 (surfaces) = 21
  EXPECT_EQ(d.transform_store().size(), 15u);

  // Identical transforms are shared between surfaces:
  // - identity: prefilled rectangle, portal cylinders and passive cylinder
  // - z = 1000mm: annuli and trapezoid
  // - z = 20mm: volume placement and the rectangle that is added last
  const std::vector<dindex> trf_links{1u, 2u,  3u, 1u,  1u,  5u,  6u,
                                      7u, 7u,  8u, 9u,  10u, 7u,  1u,
                                      11u, 12u, 13u, 4u,  14u};
  ASSERT_EQ(trf_links.size(), d.surfaces().size());
  for (const auto [idx, trf_link] : detray::views::enumerate(trf_links)) {
    geometry::identifier geo_id{};
    geo_id.set_index(idx);
    EXPECT_EQ(d.surface(geo_id).transform(), trf_link)
        << "error at index: " << idx;
  }

  // No transform appears twice in the store
  for (dindex i = 0u; i < d.transform_store().size(); ++i) {
    for (dindex j = i + 1u; j < d.transform_store().size(); ++j) {
      EXPECT_FALSE(d.transform_store().at(i) == d.transform_store().at(j))
          << "duplicate transforms at " << i << " and " << j;
    }
  }

  // Masks with identical values and volume links are shared
  EXPECT_EQ(d.mask_store().template size<mask_id::e_concentric_cylinder2D>(),
            2u);
  EXPECT_EQ(d.mask_store().template size<mask_id::e_ring2D>(), 3u);
  EXPECT_EQ(d.mask_store().template size<mask_id::e_annulus2D>(), 3u);
  EXPECT_EQ(d.mask_store().template size<mask_id::e_cylinder2D>(), 1u);
  EXPECT_EQ(d.mask_store().template size<mask_id::e_rectangle2D>(), 2u);
  EXPECT_EQ(d.mask_store().template size<mask_id::e_trapezoid2D>(), 2u);

  // check surface mask links
  std::vector<typename detector_t::surface_type::mask_link> mask_links{
      {mask_id::e_rectangle2D, {0u, 1u}},
      {mask_id::e_annulus2D, {0u, 1u}},
      {mask_id::e_trapezoid2D, {0u, 1u}},
      {mask_id::e_concentric_cylinder2D, {0u, 1u}},
      {mask_id::e_concentric_cylinder2D, {1u, 1u}},
      {mask_id::e_ring2D, {0u, 1u}},
      {mask_id::e_ring2D, {1u, 1u}},
      {mask_id::e_annulus2D, {1u, 1u}},
      {mask_id::e_annulus2D, {2u, 1u}},
      {mask_id::e_rectangle2D, {1u, 1u}},
      {mask_id::e_rectangle2D, {1u, 1u}},
      {mask_id::e_rectangle2D, {1u, 1u}},
      {mask_id::e_trapezoid2D, {1u, 1u}},
      {mask_id::e_cylinder2D, {0u, 1u}},
      {mask_id::e_ring2D, {2u, 1u}},
      {mask_id::e_ring2D, {2u, 1u}},
      {mask_id::e_rectangle2D, {1u, 1u}},
      {mask_id::e_rectangle2D, {1u, 1u}},
      {mask_id::e_rectangle2D, {1u, 1u}}};
  for (const auto [idx, m_link] : detray::views::enumerate(mask_links)) {
    geometry::identifier geo_id{};
    geo_id.set_index(idx);
    EXPECT_EQ(d.surface(geo_id).mask(), m_link) << "error at index: " << idx;
  }

  // check mask volume links
  check_mask<detector_t, mask_id::e_concentric_cylinder2D>(d, {0u, 2u});
  check_mask<detector_t, mask_id::e_cylinder2D>(d, {1u});
  check_mask<detector_t, mask_id::e_ring2D>(d, {3u, 4u, 1u});
  check_mask<detector_t, mask_id::e_annulus2D>(d, {0u, 1u, 1u});
  check_mask<detector_t, mask_id::e_rectangle2D>(d, {0u, 1u});
  check_mask<detector_t, mask_id::e_trapezoid2D>(d, {0u, 1u});
}

/// Integration test for the deduplication of surfaces with multiple masks
GTEST_TEST(detray_builders, tracking_volume_multi_mask_deduplication) {
  using namespace detray;

  using metadata_t = test::default_metadata;
  using detector_t = host::detector<metadata_t>;
  using transform3 = typename detector_t::transform3_type;
  using mask_id = typename detector_t::masks::id;
  using nav_link_t = typename detector_t::surface_type::navigation_link;
  using rectangle_factory = surface_factory<detector_t, rectangle2D>;
  using sf_data_t = typename rectangle_factory::surface_data_t;

  vecmem::host_memory_resource host_mr;
  detector_t d(host_mr);
  auto geo_ctx = typename detector_t::geometry_context{};
  const volume_builder_options builder_opts{};

  const std::vector<scalar> r1{1.f, 2.f};
  const std::vector<scalar> r2{3.f, 4.f};
  const std::vector<scalar> r3{5.f, 6.f};

  // Build a volume with mask ranges that have to be compared as a block
  auto build_volume = [&](std::vector<std::vector<std::vector<scalar>>> bounds) {
    volume_builder<detector_t> vbuilder{volume_id::e_cuboid};
    vbuilder.add_volume_placement(point3{0.f, 0.f, 0.f});

    auto rect_factory = std::make_shared<rectangle_factory>();
    for (const auto [i, sf_bounds] : detray::views::enumerate(bounds)) {
      // Unique surface placements (different from the volume placement)
      const scalar z{static_cast<scalar>(10u * (d.surfaces().size() + i) + 1u)};
      rect_factory->push_back(sf_data_t{
          surface_id::e_sensitive, transform3(point3{0.f, 0.f, z}),
          std::vector<nav_link_t>(sf_bounds.size(), 0u), sf_bounds});
    }
    vbuilder.add_surfaces(rect_factory, geo_ctx);
    vbuilder.build(d, builder_opts);
  };

  // Surfaces: {r1, r2}, {r1, r2}, {r2, r3}, {r2}, {r2, r3}
  build_volume({{r1, r2}, {r1, r2}, {r2, r3}, {r2}, {r2, r3}});
  // The second volume reuses the masks of the first volume
  build_volume({{r1, r2}, {r2, r3}});

  ASSERT_EQ(d.surfaces().size(), 7u);

  // {r2, r3} is not found in {r1, r2} as a block, so r2 is added again
  EXPECT_EQ(d.mask_store().template size<mask_id::e_rectangle2D>(), 4u);

  std::vector<typename detector_t::surface_type::mask_link> mask_links{
      {mask_id::e_rectangle2D, {0u, 2u}}, {mask_id::e_rectangle2D, {0u, 2u}},
      {mask_id::e_rectangle2D, {2u, 2u}}, {mask_id::e_rectangle2D, {1u, 1u}},
      {mask_id::e_rectangle2D, {2u, 2u}}, {mask_id::e_rectangle2D, {0u, 2u}},
      {mask_id::e_rectangle2D, {2u, 2u}}};
  for (const auto [idx, m_link] : detray::views::enumerate(mask_links)) {
    geometry::identifier geo_id{};
    geo_id.set_index(idx);
    EXPECT_EQ(d.surface(geo_id).mask(), m_link) << "error at index: " << idx;
  }

  // Check the mask values that the surfaces point to
  const auto& rect_masks = d.mask_store().template get<mask_id::e_rectangle2D>();
  const std::vector<std::vector<scalar>> expected_vals{r1, r2, r2, r3};
  for (const auto [idx, vals] : detray::views::enumerate(expected_vals)) {
    EXPECT_EQ(rect_masks.at(idx)[0], vals[0]) << "mask no. " << idx;
    EXPECT_EQ(rect_masks.at(idx)[1], vals[1]) << "mask no. " << idx;
  }

  // Both volume placements are the identity
  EXPECT_EQ(d.volumes().at(0).transform(), d.volumes().at(1).transform());
  // One transform per surface and one for the volumes
  EXPECT_EQ(d.transform_store().size(), 8u);
}
