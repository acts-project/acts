// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Detray include(s)
#include "detray/builders/grid_builder.hpp"

#include "detray/builders/surface_factory.hpp"
#include "detray/builders/volume_builder.hpp"
#include "detray/core/detector.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes/annulus2D.hpp"
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

test::transform3 Identity{};

using detector_t = detector<test::toy_metadata>;
using algebra_t = typename detector_t::algebra_type;
using scalar = dscalar<algebra_t>;
using point3 = dpoint3D<algebra_t>;
using vector3 = dvector3D<algebra_t>;

/// Mock volume builder for unit testing
struct mock_volume_builder : public volume_builder_interface<detector_t> {
  auto vol_index() const -> dindex override { return 0; }
  void has_accel(bool) override { /*Do nothing*/ }
  bool has_accel() const override { return false; }
  void set_name(std::string) override { /*Do nothing*/ }
  std::string_view name() override { return "test_volume_0"; }
  auto operator()() const -> const typename detector_t::volume_type& override {
    return m_vol;
  }
  auto operator()() -> typename detector_t::volume_type& override {
    return m_vol;
  }
  auto build(detector_t& /*unused*/,
             typename detector_t::geometry_context /*mask*/ = {}) ->
      typename detector_t::volume_type* override {
    return &m_vol;
  }
  void add_volume_placement(
      const typename detector_t::transform3_type& /*mask*/ = {})
      override { /*Do nothing*/ }
  void add_volume_placement(const typename detector_t::point3_type& /*unused*/)
      override { /*Do nothing*/ }
  void add_volume_placement(const typename detector_t::point3_type& /*unused*/,
                            const typename detector_t::vector3_type& /*unused*/,
                            const typename detector_t::vector3_type& /*unused*/)
      override { /*Do nothing*/ }
  void add_surfaces(
      std::shared_ptr<surface_factory_interface<detector_t>>,
      typename detector_t::geometry_context = {}) override { /*Do nothing*/ }

 protected:
  typename detector_t::surface_lookup_container& surfaces() override {
    return m_sf_lookup;
  }
  typename detector_t::transform_container& transforms() override {
    return m_transforms;
  }
  typename detector_t::mask_container& masks() override { return m_masks; }

 private:
  detector_t::volume_type m_vol{};
  detector_t::surface_lookup_container m_sf_lookup{};
  detector_t::transform_container m_transforms{};
  detector_t::mask_container m_masks{};
};

}  // anonymous namespace

/// Unittest: Test the construction of a collection of grids
GTEST_TEST(detray_builders, grid_factory_static) {
  // Data-owning grid collection
  vecmem::host_memory_resource host_mr;
  auto gr_factory =
      grid_factory<bins::static_array<dindex, 3>, simple_serializer, algebra_t>{
          host_mr};

  // Build from existing mask of a surface
  const scalar minR{0.f};
  const scalar maxR{10.f};
  const scalar minPhi{0.f};
  const scalar maxPhi{constant<scalar>::pi};
  mask<annulus2D, algebra_t> ann2{0u,     minR, maxR, minPhi,
                                  maxPhi, 0.f,  0.f,  0.f};

  // Grid with correctly initialized axes, but empty bin content
  auto ann_gr = gr_factory.new_grid(ann2, {5, 10});

  // Test axis
  const auto& ann_axis_r = ann_gr.template get_axis<label::e_r>();
  EXPECT_EQ(ann_axis_r.label(), label::e_r);
  EXPECT_EQ(ann_axis_r.bounds(), bounds::e_closed);
  EXPECT_EQ(ann_axis_r.binning(), binning::e_regular);
  EXPECT_EQ(ann_axis_r.nbins(), 5u);
  EXPECT_NEAR(ann_axis_r.span()[0], 0.f,
              std::numeric_limits<scalar>::epsilon());
  EXPECT_NEAR(ann_axis_r.span()[1], 10.f,
              std::numeric_limits<scalar>::epsilon());

  // Test fill a bin to see, if bin content was correctly initialized
  point3 p = {0.5f, 2.f, 0.f};
  vector3 d{};
  auto loc_p = ann_gr.project(Identity, p, d);
  ann_gr.template populate<attach<>>(loc_p, 3u);
  ann_gr.template populate<attach<>>(loc_p, 5u);
  auto bin2 = ann_gr.bin(loc_p);

  EXPECT_TRUE(bin2.size() == 2u);
  EXPECT_FALSE(bin2.empty());
  EXPECT_EQ(bin2[0], 3u);
  EXPECT_EQ(bin2[1], 5u);

  // Build from parameters
  const std::vector<scalar> bin_edges_z{-10.f, -8.f, -6.5f, -1.f,
                                        4.f,   5.f,  6.f,   9.f};
  const std::vector<scalar> bin_edges_phi{};

  auto cyl_gr = gr_factory.template new_grid<cylinder2D>(
      std::vector<scalar>{0.f, 2.f * constant<scalar>::pi, bin_edges_z.front(),
                          bin_edges_z.back()},
      {10u, bin_edges_z.size() - 1}, {}, {bin_edges_phi, bin_edges_z},
      types::list<circular<label::e_rphi>, closed<label::e_cyl_z>>{},
      types::list<regular<scalar>, irregular<scalar>>{});

  // Test axis
  const auto& cyl_axis_rphi = cyl_gr.template get_axis<label::e_rphi>();
  EXPECT_EQ(cyl_axis_rphi.label(), label::e_rphi);
  EXPECT_EQ(cyl_axis_rphi.bounds(), bounds::e_circular);
  EXPECT_EQ(cyl_axis_rphi.binning(), binning::e_regular);
  EXPECT_EQ(cyl_axis_rphi.nbins(), 10u);
  EXPECT_NEAR(cyl_axis_rphi.span()[0], 0.f,
              std::numeric_limits<scalar>::epsilon());
  EXPECT_NEAR(cyl_axis_rphi.span()[1], 2.f * constant<scalar>::pi,
              std::numeric_limits<scalar>::epsilon());

  const auto& cyl_axis_z = cyl_gr.template get_axis<label::e_cyl_z>();
  EXPECT_EQ(cyl_axis_z.label(), label::e_cyl_z);
  EXPECT_EQ(cyl_axis_z.bounds(), bounds::e_closed);
  EXPECT_EQ(cyl_axis_z.binning(), binning::e_irregular);
  EXPECT_EQ(cyl_axis_z.nbins(), 7u);
  EXPECT_NEAR(cyl_axis_z.span()[0], -10.f,
              std::numeric_limits<scalar>::epsilon());
  EXPECT_NEAR(cyl_axis_z.span()[1], 9.f,
              std::numeric_limits<scalar>::epsilon());

  // Test fill a bin to see, if bin content was correctly initialized
  loc_p = cyl_gr.project(Identity, p, d);
  cyl_gr.template populate<attach<>>(loc_p, 33u);
  cyl_gr.template populate<attach<>>(loc_p, 55u);

  EXPECT_EQ(cyl_gr.bin(loc_p)[0], 33u);
  EXPECT_EQ(cyl_gr.bin(loc_p)[1], 55u);

  // Build the same cylinder grid from a mask
  const scalar r{5.f};
  const scalar n_half_z{-10.f};
  const scalar p_half_z{9.f};
  mask<cylinder2D, algebra_t> cyl2{0u, r, n_half_z, p_half_z};

  auto cyl_gr2 =
      gr_factory
          .template new_grid<circular<label::e_rphi>, closed<label::e_cyl_z>,
                             regular<scalar>, irregular<scalar>>(
              cyl2, {bin_edges_z.size() - 1, 10u}, {},
              {bin_edges_phi, bin_edges_z});
}

/// Unittest: Test the construction of a collection of grids
GTEST_TEST(detray_builders, grid_factory_dynamic) {
  // Data-owning grid collection
  vecmem::host_memory_resource host_mr;
  auto gr_factory =
      grid_factory<bins::dynamic_array<dindex>, simple_serializer, algebra_t>{
          host_mr};

  // Build from existing mask of a surface
  const scalar minR{0.f};
  const scalar maxR{10.f};
  const scalar minPhi{0.f};
  const scalar maxPhi{constant<scalar>::pi};
  mask<annulus2D, algebra_t> ann2{0u,     minR, maxR, minPhi,
                                  maxPhi, 0.f,  0.f,  0.f};

  // Grid with correctly initialized axes and bins, but empty bin content
  using ann_grid_t =
      typename decltype(gr_factory)::template grid_type<annulus2D>;
  std::vector<std::pair<typename ann_grid_t::loc_bin_index, dindex>> capacities;

  auto bin_indexer2D = detray::views::cartesian_product{
      detray::views::iota{0u, 2u}, detray::views::iota{0u, 3u}};
  dindex capacity{0u};
  for (const auto [bin_idx0, bin_idx1] : bin_indexer2D) {
    typename ann_grid_t::loc_bin_index mbin{bin_idx0, bin_idx1};
    capacities.emplace_back(mbin, ++capacity);
  }

  ann_grid_t ann_gr = gr_factory.new_grid(ann2, {2, 3}, capacities);

  // Test axis
  const auto& ann_axis_r = ann_gr.template get_axis<label::e_r>();
  EXPECT_EQ(ann_axis_r.label(), label::e_r);
  EXPECT_EQ(ann_axis_r.bounds(), bounds::e_closed);
  EXPECT_EQ(ann_axis_r.binning(), binning::e_regular);
  EXPECT_EQ(ann_axis_r.nbins(), 2u);
  EXPECT_NEAR(ann_axis_r.span()[0], 0.f,
              std::numeric_limits<scalar>::epsilon());
  EXPECT_NEAR(ann_axis_r.span()[1], 10.f,
              std::numeric_limits<scalar>::epsilon());

  // Test bin storage
  EXPECT_EQ(ann_gr.bins().size(), 6u);
  for (const auto& bin : ann_gr.bins()) {
    EXPECT_EQ(bin.size(), 0u);
    EXPECT_TRUE(bin.capacity() >= 1u);
  }
  EXPECT_EQ(ann_gr.size(), 0u);

  EXPECT_EQ(ann_gr.bins().entry_data().size(), 21u);
  EXPECT_EQ(ann_gr.bins().entry_data().capacity(), 21u);

  // Test fill a bin to see, if bin content was correctly initialized
  point3 p = {0.5f, 2.f, 0.f};
  vector3 d{};
  auto loc_p = ann_gr.project(Identity, p, d);
  auto bin2 = ann_gr.bin(loc_p);
  EXPECT_TRUE(bin2.capacity() == 2u);
  EXPECT_TRUE(bin2.empty());
  EXPECT_TRUE(bin2.size() == 0u);
  ann_gr.template populate<attach<>>(loc_p, 3u);
  ann_gr.template populate<attach<>>(loc_p, 5u);

  EXPECT_TRUE(bin2.capacity() == 2u);
  EXPECT_TRUE(bin2.size() == 2u);
  EXPECT_FALSE(bin2.empty());
  EXPECT_EQ(bin2[0], 3u);
  EXPECT_EQ(bin2[1], 5u);

  EXPECT_EQ(ann_gr.size(), 2u);

  // Build from parameters
  const std::vector<scalar> bin_edges_z{-10.f, -8.f, -6.5f, -1.f,
                                        4.f,   5.f,  6.f,   9.f};
  const std::vector<scalar> bin_edges_phi{};

  bin_indexer2D = detray::views::cartesian_product{detray::views::iota{0u, 10u},
                                                   detray::views::iota{0u, 7u}};
  capacity = 0u;
  capacities.clear();
  for (const auto [bin_idx0, bin_idx1] : bin_indexer2D) {
    typename ann_grid_t::loc_bin_index mbin{bin_idx0, bin_idx1};
    capacities.emplace_back(mbin, ++capacity);
  }

  auto cyl_gr = gr_factory.template new_grid<concentric_cylinder2D>(
      std::vector<scalar>{0.f, 2.f * constant<scalar>::pi, bin_edges_z.front(),
                          bin_edges_z.back()},
      {10u, bin_edges_z.size() - 1}, capacities, {bin_edges_phi, bin_edges_z},
      types::list<circular<label::e_rphi>, closed<label::e_cyl_z>>{},
      types::list<regular<scalar>, irregular<scalar>>{});

  // Test axis
  const auto& cyl_axis_rphi = cyl_gr.template get_axis<label::e_rphi>();
  EXPECT_EQ(cyl_axis_rphi.label(), label::e_rphi);
  EXPECT_EQ(cyl_axis_rphi.bounds(), bounds::e_circular);
  EXPECT_EQ(cyl_axis_rphi.binning(), binning::e_regular);
  EXPECT_EQ(cyl_axis_rphi.nbins(), 10u);
  EXPECT_NEAR(cyl_axis_rphi.span()[0], 0.f,
              std::numeric_limits<scalar>::epsilon());
  EXPECT_NEAR(cyl_axis_rphi.span()[1], 2.f * constant<scalar>::pi,
              std::numeric_limits<scalar>::epsilon());

  const auto& cyl_axis_z = cyl_gr.template get_axis<label::e_cyl_z>();
  EXPECT_EQ(cyl_axis_z.label(), label::e_cyl_z);
  EXPECT_EQ(cyl_axis_z.bounds(), bounds::e_closed);
  EXPECT_EQ(cyl_axis_z.binning(), binning::e_irregular);
  EXPECT_EQ(cyl_axis_z.nbins(), 7u);
  EXPECT_NEAR(cyl_axis_z.span()[0], -10.f,
              std::numeric_limits<scalar>::epsilon());
  EXPECT_NEAR(cyl_axis_z.span()[1], 9.f,
              std::numeric_limits<scalar>::epsilon());

  // Test bin storage
  EXPECT_EQ(cyl_gr.bins().size(), 70u);
  for (const auto& bin : cyl_gr.bins()) {
    EXPECT_EQ(bin.size(), 0u);
    EXPECT_TRUE(bin.capacity() >= 1u);
  }
  EXPECT_EQ(cyl_gr.size(), 0u);

  EXPECT_EQ(cyl_gr.bins().entry_data().size(), 2485u);
  EXPECT_EQ(cyl_gr.bins().entry_data().capacity(), 2485u);

  // Test fill a bin to see, if bin content was correctly initialized
  loc_p = cyl_gr.project(Identity, p, d);
  auto bin3 = cyl_gr.bin(loc_p);

  EXPECT_TRUE(bin3.capacity() == 18u);
  EXPECT_TRUE(bin3.empty());
  EXPECT_TRUE(bin3.size() == 0u);

  cyl_gr.template populate<attach<>>(loc_p, 33u);
  cyl_gr.template populate<attach<>>(loc_p, 55u);
  cyl_gr.template populate<attach<>>(loc_p, 77u);

  EXPECT_TRUE(bin3.capacity() == 18u);
  EXPECT_FALSE(bin3.empty());
  EXPECT_TRUE(bin3.size() == 3u);
  EXPECT_EQ(bin3[0], 33u);
  EXPECT_EQ(bin3[1], 55u);
  EXPECT_EQ(bin3[2], 77u);

  EXPECT_EQ(cyl_gr.size(), 3u);

  // Build the same cylinder grid from a mask
  const scalar r{5.f};
  const scalar n_half_z{-10.f};
  const scalar p_half_z{9.f};
  mask<cylinder2D, algebra_t> cyl2{0u, r, n_half_z, p_half_z};

  auto cyl_gr2 =
      gr_factory
          .template new_grid<circular<label::e_rphi>, closed<label::e_cyl_z>,
                             regular<scalar>, irregular<scalar>>(
              cyl2, {10u, bin_edges_z.size() - 1}, capacities,
              {bin_edges_phi, bin_edges_z});
}

/// Unittest: Test the grid builder
GTEST_TEST(detray_builders, grid_builder) {
  // cylinder grid type of the toy detector
  using cyl_grid_t = grid<algebra_t, axes<concentric_cylinder2D>,
                          bins::static_array<detector_t::surface_type, 1>,
                          simple_serializer, host_container_types, false>;

  auto gbuilder = grid_builder<detector_t, cyl_grid_t, detray::bin_associator>{
      std::make_unique<mock_volume_builder>()};

  // The cylinder portals are at the end of the surface range by construction
  const auto cyl_mask =
      mask<concentric_cylinder2D, algebra_t>{0u, 10.f, -500.f, 500.f};
  std::size_t n_phi_bins{5u};
  std::size_t n_z_bins{4u};

  // Build empty grid
  gbuilder.init_grid(cyl_mask, {n_phi_bins, n_z_bins});
  auto cyl_grid = gbuilder.get();

  const auto& cyl_axis_z = cyl_grid.template get_axis<label::e_cyl_z>();
  EXPECT_EQ(cyl_axis_z.label(), label::e_cyl_z);
  EXPECT_EQ(cyl_axis_z.bounds(), bounds::e_closed);
  EXPECT_EQ(cyl_axis_z.binning(), binning::e_regular);
  EXPECT_EQ(cyl_axis_z.nbins(), 4u);
  EXPECT_NEAR(cyl_axis_z.span()[0], -500.f,
              std::numeric_limits<scalar>::epsilon());
  EXPECT_NEAR(cyl_axis_z.span()[1], 500.f,
              std::numeric_limits<scalar>::epsilon());
}
