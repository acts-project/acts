// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Detray include(s)
#include "detray/definitions/indexing.hpp"
#include "detray/geometry/shapes/cuboid3D.hpp"
#include "detray/navigation/accelerators/concepts.hpp"
#include "detray/utils/grid/concepts.hpp"
#include "detray/utils/grid/grid.hpp"

// Detray test include(s)
#include "detray/test/framework/types.hpp"

// Vecmem include(s)
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/containers/vector.hpp>
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <algorithm>
#include <limits>
#include <random>

using namespace detray;
using namespace detray::axis;

namespace {

// Algebra definitions
using test_algebra = test::algebra;
using scalar = test::scalar;
using point3 = test::point3;

// Either a data owning or non-owning 3D cartesian multi-axis
template <bool ownership = true, typename containers = host_container_types>
using cartesian_3D =
    coordinate_axes<axes<cuboid3D>, test_algebra, ownership, containers>;

// non-owning multi-axis: Takes external containers
bool constexpr is_owning = true;
bool constexpr is_n_owning = false;

// Bin edges for all axes
dvector<scalar> bin_edges = {-10.f, 10.f, -20.f, 20.f, 0.f, 100.f};
// Offsets into edges container and #bins for all axes
dvector<dsized_index_range> edge_ranges = {{0u, 20u}, {2u, 40u}, {4u, 50u}};

// non-owning multi-axis for the non-owning grid
cartesian_3D<is_n_owning, host_container_types> ax_n_own(edge_ranges,
                                                         bin_edges);

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

/// Test bin entry retrieval
GTEST_TEST(detray_acceleration_structures, spatial_grid) {
  // Non-owning, 3D cartesian, replacing grid
  using grid_t = grid<test_algebra, axes<cuboid3D>, bins::single<scalar>>;
  using spatial_grid_t = spatial_grid_impl<grid_t>;

  static_assert(concepts::grid<spatial_grid_t>);
  static_assert(concepts::accelerator<spatial_grid_t>);

  // Fill the bin data for every test
  // bin test entries
  spatial_grid_t::bin_container_type bin_data{};
  bin_data.resize(40'000u);
  std::ranges::generate_n(
      bin_data.begin(), 40'000u,
      bin_content_sequence<replace<>, bins::single<scalar>>());

  // Copy data that will be moved into the data owning types
  dvector<scalar> bin_edges_cp(bin_edges);
  dvector<dsized_index_range> edge_ranges_cp(edge_ranges);
  spatial_grid_t::bin_container_type bin_data_cp(bin_data);

  // Data-owning axes and grid
  cartesian_3D<is_owning, host_container_types> axes_own(
      std::move(edge_ranges_cp), std::move(bin_edges_cp));
  spatial_grid_t grid_3D(std::move(bin_data_cp), std::move(axes_own));

  // Test the bin view
  point3 p = {-10.f, -20.f, 0.f};

  std::array<dindex, 2> search_window_size{0, 0};
  axis::multi_bin_range<3> search_window{
      axis::bin_range{0, 1}, axis::bin_range{0, 1}, axis::bin_range{0, 1}};

  const auto bview1 = axis::detail::bin_view(grid_3D, search_window);
  const auto joined_view1 = detray::views::join(bview1);
  const auto grid_search1 = grid_3D.search(p, search_window_size);

  static_assert(detray::ranges::bidirectional_range<decltype(bview1)>);
  static_assert(detray::ranges::bidirectional_range<decltype(joined_view1)>);
  static_assert(detray::ranges::bidirectional_range<decltype(grid_search1)>);

  ASSERT_EQ(bview1.size(), 1u);
  ASSERT_EQ(joined_view1.size(), 1u);
  ASSERT_EQ(grid_search1.size(), 1u);

  for (auto bin : bview1) {
    for (auto entry : bin) {
      EXPECT_EQ(entry, grid_3D.bin(0).value()) << "bin entry: " << entry;
    }
  }

  for (scalar entry : joined_view1) {
    EXPECT_EQ(entry, grid_3D.bin(0).value()) << "bin entry: " << entry;
  }

  for (scalar entry : grid_search1) {
    EXPECT_EQ(entry, grid_3D.bin(0).value()) << "bin entry: " << entry;
  }

  //
  // In the corner of the grid, include nearest neighbor
  //
  search_window_size[0] = 1;
  search_window_size[1] = 1;
  search_window[0] = axis::bin_range{0, 2};
  search_window[1] = axis::bin_range{0, 2};
  search_window[2] = axis::bin_range{0, 2};

  std::vector<scalar> expected{1, 801, 21, 821, 2, 802, 22, 822};

  const auto bview2 = axis::detail::bin_view(grid_3D, search_window);
  const auto joined_view2 = detray::views::join(bview2);
  const auto grid_search2 = grid_3D.search(p, search_window_size);

  static_assert(detray::ranges::bidirectional_range<decltype(bview2)>);
  static_assert(detray::ranges::bidirectional_range<decltype(joined_view2)>);
  static_assert(detray::ranges::bidirectional_range<decltype(grid_search2)>);

  ASSERT_EQ(bview2.size(), 8u);
  ASSERT_EQ(joined_view2.size(), 8u);
  ASSERT_EQ(grid_search2.size(), 8u);

  for (auto [i, bin] : detray::views::enumerate(bview2)) {
    for (scalar entry : bin) {
      EXPECT_EQ(entry, expected[i]) << "bin entry: " << entry;
    }
  }

  for (auto [i, entry] : detray::views::enumerate(joined_view2)) {
    EXPECT_EQ(entry, expected[i]) << "bin entry: " << entry;
  }

  for (auto [i, entry] : detray::views::enumerate(grid_search2)) {
    EXPECT_EQ(entry, expected[i]) << "bin entry: " << entry;
  }

  //
  // Bin 1 and nearest neighbors
  //
  p = {-9.f, -19.f, 2.f};
  search_window[0] = axis::bin_range{0, 3};
  search_window[1] = axis::bin_range{0, 3};
  search_window[2] = axis::bin_range{0, 3};

  expected = {1, 801, 1601, 21, 821, 1621, 41, 841, 1641,
              2, 802, 1602, 22, 822, 1622, 42, 842, 1642,
              3, 803, 1603, 23, 823, 1623, 43, 843, 1643};

  const auto bview3 = axis::detail::bin_view(grid_3D, search_window);
  const auto joined_view3 = detray::views::join(bview3);
  const auto grid_search3 = grid_3D.search(p, search_window_size);

  static_assert(detray::ranges::bidirectional_range<decltype(bview3)>);
  static_assert(detray::ranges::bidirectional_range<decltype(joined_view3)>);
  static_assert(detray::ranges::bidirectional_range<decltype(grid_search3)>);

  ASSERT_EQ(bview3.size(), 27u);
  ASSERT_EQ(joined_view3.size(), 27u);
  ASSERT_EQ(grid_search3.size(), 27u);

  for (auto [i, bin] : detray::views::enumerate(bview3)) {
    for (scalar entry : bin) {
      EXPECT_EQ(entry, expected[i]) << "bin entry: " << entry;
    }
  }

  for (auto [i, entry] : detray::views::enumerate(joined_view3)) {
    EXPECT_EQ(entry, expected[i]) << "bin entry: " << entry;
  }

  for (auto [i, entry] : detray::views::enumerate(grid_search3)) {
    EXPECT_EQ(entry, expected[i]) << "bin entry: " << entry;
  }
}

/// Test bin entry retrieval
GTEST_TEST(detray_acceleration_structures, spatial_grid_complete_population) {
  // Non-owning, 3D cartesian, completing grid (4 dims and sort)
  using grid_t =
      grid<test_algebra, decltype(ax_n_own), bins::static_array<scalar, 4>,
           simple_serializer, host_container_types, false>;
  using spatial_grid_t = spatial_grid_impl<grid_t>;

  static_assert(concepts::grid<spatial_grid_t>);
  static_assert(concepts::accelerator<spatial_grid_t>);

  // init
  spatial_grid_t::bin_container_type bin_data{};
  bin_data.resize(40'000u, spatial_grid_t::bin_type{});

  // Create non-owning grid
  spatial_grid_t g3c(&bin_data, ax_n_own);
  // Fill and read
  point3 p = {-4.5f, -4.5f, 4.5f};
  g3c.template populate<complete<>>(p, 4.f);

  // Test search without search window
  for (scalar entry : g3c.search(p)) {
    EXPECT_EQ(entry, 4.f);
  }

  std::array<dindex, 2> search_window_size{0, 0};

  const auto grid_search1 = g3c.search(p, search_window_size);

  static_assert(detray::ranges::bidirectional_range<decltype(grid_search1)>);

  ASSERT_EQ(grid_search1.size(), 4u);

  for (scalar entry : grid_search1) {
    EXPECT_EQ(entry, 4.f);
  }

  // No neighbors were filled, expect the same candidates
  search_window_size[0] = 1;
  search_window_size[1] = 1;

  const auto grid_search2 = g3c.search(p, search_window_size);

  static_assert(detray::ranges::bidirectional_range<decltype(grid_search2)>);

  ASSERT_EQ(grid_search2.size(), 4u);

  for (scalar entry : grid_search2) {
    EXPECT_EQ(entry, 4.f);
  }

  // Populate some neighboring bins
  point3 p2 = {-5.5f, -4.5f, 4.5f};
  g3c.template populate<complete<>>(p2, 4.f);
  p2 = {-3.5f, -4.5f, 4.5f};
  g3c.template populate<complete<>>(p2, 4.f);
  p2 = {-4.5f, -5.5f, 4.5f};
  g3c.template populate<complete<>>(p2, 4.f);
  p2 = {-4.5f, -3.5f, 4.5f};
  g3c.template populate<complete<>>(p2, 4.f);

  const auto grid_search3 = g3c.search(p, search_window_size);
  ASSERT_EQ(grid_search3.size(), 20u);

  for (scalar entry : grid_search3) {
    EXPECT_EQ(entry, 4.f);
  }
}

/// Test bin entry retrieval
GTEST_TEST(detray_acceleration_structures,
           spatial_grid_regular_attach_population) {
  // Non-owning, 3D cartesian, completing grid (4 dims and sort)
  using grid_t =
      grid<test_algebra, decltype(ax_n_own), bins::static_array<scalar, 4>,
           simple_serializer, host_container_types, false>;
  using spatial_grid_t = spatial_grid_impl<grid_t>;

  static_assert(concepts::grid<spatial_grid_t>);
  static_assert(concepts::accelerator<spatial_grid_t>);

  // init
  spatial_grid_t::bin_container_type bin_data{};
  bin_data.resize(40'000u, spatial_grid_t::bin_type{});

  // Create non-owning grid
  spatial_grid_t g3ra(&bin_data, ax_n_own);
  // Fill and read
  point3 p = {-4.5f, -4.5f, 4.5f};
  g3ra.template populate<attach<>>(p, 5.f);

  // Test search without search window
  for (scalar entry : g3ra.search(p)) {
    EXPECT_EQ(entry, 5.f);
  }

  std::array<dindex, 2> search_window_size{0, 0};

  const auto grid_search1 = g3ra.search(p, search_window_size);

  static_assert(detray::ranges::bidirectional_range<decltype(grid_search1)>);

  ASSERT_EQ(grid_search1.size(), 1u);

  for (scalar entry : grid_search1) {
    EXPECT_EQ(entry, 5.f);
  }

  // No neighbors were filled, expect the same candidates
  search_window_size[0] = 1;
  search_window_size[1] = 1;

  const auto grid_search2 = g3ra.search(p, search_window_size);

  static_assert(detray::ranges::bidirectional_range<decltype(grid_search2)>);

  ASSERT_EQ(grid_search2.size(), 1u);

  for (scalar entry : grid_search2) {
    EXPECT_EQ(entry, 5.f);
  }

  // Put more candidates into bin
  g3ra.template populate<attach<>>(p, 6.f);
  g3ra.template populate<attach<>>(p, 7.f);
  std::vector<scalar> entries{5.f, 6.f, 7.f};

  const auto grid_search3 = g3ra.search(p, search_window_size);
  ASSERT_EQ(grid_search3.size(), 3u);

  for (auto [i, entry] : detray::views::enumerate(grid_search3)) {
    EXPECT_EQ(entry, entries[i]);
  }

  // Populate some neighboring bins
  point3 p2 = {-5.5f, -4.5f, 4.5f};
  g3ra.template populate<attach<>>(p2, 5.f);
  p2 = {-3.5f, -4.5f, 4.5f};
  g3ra.template populate<attach<>>(p2, 5.f);
  p2 = {-4.5f, -5.5f, 4.5f};
  g3ra.template populate<attach<>>(p2, 5.f);
  p2 = {-4.5f, -3.5f, 4.5f};
  g3ra.template populate<attach<>>(p2, 5.f);

  const auto grid_search4 = g3ra.search(p, search_window_size);
  ASSERT_EQ(grid_search4.size(), 7u);

  for (scalar entry : grid_search4) {
    EXPECT_TRUE(entry == 5.f || entry == 6.f || entry == 7.f);
  }
}
