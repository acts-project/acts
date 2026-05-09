// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Detray include(s)
#include "detray/utils/grid/grid_collection.hpp"

#include "detray/definitions/indexing.hpp"
#include "detray/geometry/shapes/cylinder3D.hpp"
#include "detray/utils/grid/grid.hpp"
#include "detray/utils/grid/populators.hpp"
#include "detray/utils/grid/serializers.hpp"

// Detray test include(s)
#include "detray/test/framework/types.hpp"

// System include(s)
#include <algorithm>
#include <limits>

// GTest include(s)
#include <gtest/gtest.h>

using namespace detray;
using namespace detray::axis;

using test_algebra = test::algebra;
using scalar = test::scalar;

namespace {

// non-owning multi-axis: Takes external containers
bool constexpr is_n_owning = false;

constexpr dindex inf{std::numeric_limits<dindex>::max()};

// Create some bin data for non-owning grid
template <typename populator_t, typename bin_t>
struct bin_content_sequence {
  using entry_t = typename bin_t::entry_type;
  entry_t entry{0};

  auto operator()() {
    entry += entry_t{1};
    bin_t bin{};
    populator_t{}(bin, entry);
    return bin;
  }
};

}  // anonymous namespace

/// Unittest: Test the construction of a collection of grids
GTEST_TEST(detray_grid, grid_collection) {
  // Non-owning grid type with array<dindex, 3> as bin content
  using grid_t =
      grid<test_algebra, axes<cylinder3D>, bins::static_array<dindex, 3>,
           simple_serializer, host_container_types, is_n_owning>;

  // Build test data
  grid_t::bin_container_type bin_data{};
  bin_data.resize(197u);
  std::ranges::generate_n(
      bin_data.begin(), 197u,
      bin_content_sequence<attach<>, typename grid_t::bin_type>());
  dvector<dindex> grid_offsets = {0u, 48u, 72u};

  // Offsets into edges container and #bins for all axes
  dvector<dsized_index_range> edge_ranges = {{0u, 2u},  {2u, 4u},  {4u, 6u},
                                             {6u, 1u},  {8u, 3u},  {10u, 8u},
                                             {12u, 5u}, {14u, 5u}, {16u, 5u}};

  // Bin edges for all axes
  dvector<scalar> bin_edges = {-10.f, 10.f, -20.f, 20.f, 0.f, 120.f,
                               -5.f,  5.f,  -15.f, 15.f, 0.f, 50.f,
                               -15.f, 15.f, -35.f, 35.f, 0.f, 550.f};

  // Data-owning grid collection
  auto grid_coll =
      grid_collection<grid_t>(std::move(grid_offsets), std::move(bin_data),
                              std::move(edge_ranges), std::move(bin_edges));

  // Tests

  // Basics
  EXPECT_EQ(grid_coll.size(), 3u);
  EXPECT_EQ(grid_coll.bin_storage().size(), 197u);
  EXPECT_EQ(grid_coll.axes_storage().size(), 9u);
  EXPECT_EQ(grid_coll.bin_edges_storage().size(), 18u);

  // Get a grid instance
  auto single_grid = grid_coll[1];

  static_assert(std::is_same_v<decltype(single_grid), grid_t>,
                "Grid from collection has wrong type");

  EXPECT_EQ(single_grid.dim, 3);
  EXPECT_EQ(single_grid.nbins(), 24u);
  auto r_axis = single_grid.get_axis<label::e_r>();
  EXPECT_EQ(r_axis.nbins(), 1u);
  auto phi_axis = single_grid.get_axis<label::e_phi>();
  EXPECT_EQ(phi_axis.nbins(), 3u);
  using z_axis_t = single_axis<closed<label::e_z>, regular<scalar>>;
  auto z_axis = single_grid.get_axis<z_axis_t>();
  EXPECT_EQ(z_axis.nbins(), 8u);

  // The generator starts counting at one instead of zero
  EXPECT_EQ(single_grid.bin(0u, 0u, 0u)[0u], 49u);
  EXPECT_EQ(single_grid.bin(0u, 0u, 0u)[1u], inf);
  EXPECT_EQ(single_grid.bin(0u, 0u, 0u)[2u], inf);

  // Test the bin view
  auto& bin_view = grid_coll[2].bin(101u);
  grid_coll[2].template populate<attach<>>(101u, 42u);
  EXPECT_EQ(bin_view[0u], 102u + 72u);
  EXPECT_EQ(bin_view[1u], 42u);
  EXPECT_EQ(bin_view[2u], inf);

  // Test the global bin iteration. Take the middle grid!
  auto seq = detray::views::iota(49, 73);
  auto flat_bin_view = single_grid.all();
  EXPECT_EQ(seq.size(), 24u);
  EXPECT_EQ(flat_bin_view.size(), 24u);
  EXPECT_TRUE(
      std::equal(flat_bin_view.begin(), flat_bin_view.end(), seq.begin()));

  auto grid_coll_view = get_data(grid_coll);
  static_assert(std::is_same_v<decltype(grid_coll_view),
                               typename grid_collection<grid_t>::view_type>,
                "Grid collection view incorrectly assembled");

  const grid_collection<grid_t>& const_coll = grid_coll;
  auto const_coll_view = get_data(const_coll);
  static_assert(
      std::is_same_v<decltype(const_coll_view),
                     typename grid_collection<grid_t>::const_view_type>,
      "Grid collection const view incorrectly assembled");
}

/// Unittest: Test the construction of a collection of grids with dynamic grid
/// capacities
GTEST_TEST(detray_grid, grid_collection_dynamic_bin) {
  // Non-owning grid type with vector<dindex> as bin content
  using grid_t =
      grid<test_algebra, axes<cylinder3D>, bins::dynamic_array<dindex>,
           simple_serializer, host_container_types, is_n_owning>;

  // Build test data
  grid_t::bin_container_type bin_data{};
  bin_data.bins.resize(197u);
  bin_data.entries.resize(4u * 197u);

  int i{0};
  dindex entry{0u};
  dindex offset{0u};
  attach<> attacher{};
  for (auto& data : bin_data.bins) {
    data.offset = offset;
    // Every second bin holds one element, otherwise three
    data.capacity = (i % 2) != 0u ? 1u : 3u;

    detray::bins::dynamic_array bin{bin_data.entries.data(), data};

    ASSERT_TRUE(bin.capacity() == (i % 2 != 0u ? 1u : 3u));
    ASSERT_TRUE(bin.size() == 0);

    offset += bin.capacity();

    // Populate the bin
    attacher(bin, entry);

    ASSERT_TRUE(bin.size() == 1);

    std::size_t idx{0u};
    for (auto e : bin) {
      if (idx == 0u) {
        ASSERT_TRUE(e == entry);
      } else {
        ASSERT_TRUE(e == inf);
      }
      ++idx;
    }
    ++entry;
    ++i;
  }

  dvector<dindex> grid_offsets = {0u, 48u, 72u};

  // Offsets into edges container and #bins for all axes
  dvector<dsized_index_range> edge_ranges = {{0u, 2u},  {2u, 4u},  {4u, 6u},
                                             {6u, 1u},  {8u, 3u},  {10u, 8u},
                                             {12u, 5u}, {14u, 5u}, {16u, 5u}};

  // Bin edges for all axes
  dvector<scalar> bin_edges = {-10.f, 10.f, -20.f, 20.f, 0.f, 120.f,
                               -5.f,  5.f,  -15.f, 15.f, 0.f, 50.f,
                               -15.f, 15.f, -35.f, 35.f, 0.f, 550.f};

  // Data-owning grid collection
  auto grid_coll =
      grid_collection<grid_t>(std::move(grid_offsets), std::move(bin_data),
                              std::move(edge_ranges), std::move(bin_edges));

  // Tests

  // Basics
  EXPECT_EQ(grid_coll.size(), 3u);
  EXPECT_EQ(grid_coll.bin_storage().bins.size(), 197u);
  EXPECT_EQ(grid_coll.bin_storage().entries.size(), 4u * 197u);
  EXPECT_EQ(grid_coll.axes_storage().size(), 9u);
  EXPECT_EQ(grid_coll.bin_edges_storage().size(), 18u);

  // Get a grid instance
  auto single_grid = grid_coll[1];

  static_assert(std::is_same_v<decltype(single_grid), grid_t>,
                "Grid from collection has wrong type");

  EXPECT_EQ(single_grid.dim, 3);
  EXPECT_EQ(single_grid.nbins(), 24u);
  auto r_axis = single_grid.get_axis<label::e_r>();
  EXPECT_EQ(r_axis.nbins(), 1u);
  auto phi_axis = single_grid.get_axis<label::e_phi>();
  EXPECT_EQ(phi_axis.nbins(), 3u);
  using z_axis_t = single_axis<closed<label::e_z>, regular<scalar>>;
  auto z_axis = single_grid.get_axis<z_axis_t>();
  EXPECT_EQ(z_axis.nbins(), 8u);

  auto bin = single_grid.bin(0u, 0u, 0u);
  EXPECT_EQ(bin.capacity(), 3u);
  EXPECT_EQ(bin.size(), 1u);
  EXPECT_EQ(bin[0u], 48u);

  // Test the bin view
  EXPECT_EQ(grid_coll[2].nbins(), 125u);
  grid_coll[2].template populate<attach<>>(102u, 42u);
  auto bin_view = grid_coll[2].bin(102u);
  EXPECT_EQ(bin_view.capacity(), 3u);
  EXPECT_EQ(bin_view.size(), 2u);
  EXPECT_EQ(bin_view[0u], 102u + 72u);
  EXPECT_EQ(bin_view[1u], 42u);

  // Test the global bin iteration. Take the middle grid!
  auto seq = detray::views::iota(48, 72);
  auto flat_bin_view = single_grid.all();
  EXPECT_EQ(seq.size(), 24u);
  EXPECT_EQ(flat_bin_view.size(), 24u);
  EXPECT_TRUE(
      std::equal(flat_bin_view.begin(), flat_bin_view.end(), seq.begin()));

  auto grid_coll_view = get_data(grid_coll);
  static_assert(std::is_same_v<decltype(grid_coll_view),
                               typename grid_collection<grid_t>::view_type>,
                "Grid collection view incorrectly assembled");

  const grid_collection<grid_t>& const_coll = grid_coll;
  auto const_coll_view = get_data(const_coll);
  static_assert(
      std::is_same_v<decltype(const_coll_view),
                     typename grid_collection<grid_t>::const_view_type>,
      "Grid collection const view incorrectly assembled");
}
