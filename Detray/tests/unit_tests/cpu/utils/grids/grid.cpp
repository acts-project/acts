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
#include "detray/utils/grid/grid.hpp"

#include "detray/definitions/indexing.hpp"
#include "detray/geometry/shapes/cuboid3D.hpp"
#include "detray/utils/grid/concepts.hpp"

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

constexpr scalar inf{std::numeric_limits<scalar>::max()};
constexpr scalar tol{1e-7f};

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

/// Test bin content element by element
template <concepts::grid grid_t, typename content_t>
void test_content(const grid_t& g, const point3& p, const content_t& expected) {
  dindex i = 0u;
  for (const auto& entry : g.bin(p)) {
    ASSERT_NEAR(entry, expected[i++], tol) << " at index " << i - 1u;
  }
}

}  // anonymous namespace

/// Unittest: Test single-bin grid construction
GTEST_TEST(detray_grid, single_grid) {
  // Owning and non-owning, cartesian, 3-dimensional grids
  using grid_owning_t =
      grid<test_algebra, axes<cuboid3D>, bins::single<scalar>>;

  using grid_n_owning_t =
      grid<test_algebra, axes<cuboid3D>, bins::single<scalar>,
           simple_serializer, host_container_types, false>;

  using grid_device_t = grid<test_algebra, axes<cuboid3D>, bins::single<scalar>,
                             simple_serializer, device_container_types>;

  static_assert(concepts::grid<grid_owning_t>);
  static_assert(concepts::grid<grid_n_owning_t>);
  static_assert(concepts::grid<grid_device_t>);

  // Fill the bin data for every test
  // bin test entries
  grid_owning_t::bin_container_type bin_data{};
  bin_data.resize(40'000u);
  std::ranges::generate_n(
      bin_data.begin(), 40'000u,
      bin_content_sequence<replace<>, bins::single<scalar>>());

  // Copy data that will be moved into the data owning types
  dvector<scalar> bin_edges_cp(bin_edges);
  dvector<dsized_index_range> edge_ranges_cp(edge_ranges);
  grid_owning_t::bin_container_type bin_data_cp(bin_data);

  // Data-owning axes and grid
  cartesian_3D<is_owning, host_container_types> axes_own(
      std::move(edge_ranges_cp), std::move(bin_edges_cp));
  grid_owning_t grid_own(std::move(bin_data_cp), std::move(axes_own));

  // Copy a second time for the comparison
  dvector<scalar> bin_edges_cp2(bin_edges);
  dvector<dsized_index_range> edge_ranges_cp2(edge_ranges);
  grid_owning_t::bin_container_type bin_data_cp2(bin_data);

  // Make a second grid
  cartesian_3D<is_owning, host_container_types> axes_own2(
      std::move(edge_ranges_cp2), std::move(bin_edges_cp2));
  grid_owning_t grid_own2(std::move(bin_data_cp2), std::move(axes_own2));

  // CHECK equality
  EXPECT_TRUE(grid_own == grid_own2);

  // Check a few basics
  EXPECT_EQ(grid_own.dim, 3u);
  EXPECT_EQ(grid_own.nbins(), 40'000u);
  auto y_axis = grid_own.get_axis<label::e_y>();
  EXPECT_EQ(y_axis.nbins(), 40u);
  auto z_axis =
      grid_own.get_axis<single_axis<closed<label::e_z>, regular<scalar>>>();
  EXPECT_EQ(z_axis.nbins(), 50u);

  // Create non-owning grid
  grid_n_owning_t grid_n_own(&bin_data, ax_n_own);

  // Test for consistency with owning grid
  EXPECT_EQ(grid_n_own.dim, grid_own.dim);
  y_axis = grid_n_own.get_axis<label::e_y>();
  EXPECT_EQ(y_axis.nbins(), grid_own.get_axis<label::e_y>().nbins());
  z_axis = grid_n_own.get_axis<label::e_z>();
  EXPECT_EQ(z_axis.nbins(), grid_own.get_axis<label::e_z>().nbins());

  // Construct a grid from a view
  grid_owning_t::view_type grid_view = get_data(grid_own);
  grid_device_t device_grid(grid_view);

  // Test for consistency with non-owning grid
  EXPECT_EQ(device_grid.dim, grid_n_own.dim);
  auto y_axis_dev = device_grid.get_axis<label::e_y>();
  EXPECT_EQ(y_axis_dev.nbins(), grid_n_own.get_axis<label::e_y>().nbins());
  auto z_axis_dev = device_grid.get_axis<label::e_z>();
  EXPECT_EQ(z_axis_dev.nbins(), grid_n_own.get_axis<label::e_z>().nbins());

  // Test the global bin iteration: owning grid
  auto seq = detray::views::iota(1, 40'001);
  auto flat_bin_view = grid_own.all();

  static_assert(detray::ranges::random_access_range<decltype(flat_bin_view)>);

  EXPECT_EQ(seq.size(), 40'000u);
  EXPECT_EQ(flat_bin_view.size(), 40'000u);
  EXPECT_EQ(flat_bin_view[42], 43u);
  EXPECT_TRUE(
      std::equal(flat_bin_view.begin(), flat_bin_view.end(), seq.begin()));

  // Test the global bin iteration: non-owning grid
  auto flat_bin_view2 = grid_n_own.all();

  static_assert(detray::ranges::random_access_range<decltype(flat_bin_view2)>);

  EXPECT_EQ(seq.size(), 40'000u);
  EXPECT_EQ(flat_bin_view2.size(), 40'000u);
  EXPECT_EQ(flat_bin_view2[42], 43u);
  EXPECT_TRUE(
      std::equal(flat_bin_view2.begin(), flat_bin_view2.end(), seq.begin()));

  // Test const grid view
  /*auto const_grid_view = get_data(const_cast<const
  grid_owning_t&>(grid_own));

  static_assert(
      std::is_same_v<decltype(const_grid_view),
                     typename grid<test_algebra,cartesian_3D<>,
  bins::single<const scalar>>::view_type>, "Const grid view was not correctly
  constructed!");

  grid<test_algebra,cartesian_3D<is_owning, device_container_types>,
  bins::single<const scalar>> const_device_grid(const_grid_view);

  static_assert(
      std::is_same_v<typename decltype(const_device_grid)::bin_type,
                     typename replace::template bin_type<const scalar>>,
      "Const grid was not correctly constructed from view!");*/
}

/// Unittest: Test dynamic array grid
GTEST_TEST(detray_grid, dynamic_array) {
  // Owning and non-owning, cartesian, 3-dimensional grids
  using grid_owning_t =
      grid<test_algebra, axes<cuboid3D>, bins::dynamic_array<scalar>>;

  using grid_n_owning_t =
      grid<test_algebra, axes<cuboid3D>, bins::dynamic_array<scalar>,
           simple_serializer, host_container_types, false>;

  using grid_device_t =
      grid<test_algebra, axes<cuboid3D>, bins::dynamic_array<scalar>,
           simple_serializer, device_container_types>;

  static_assert(concepts::grid<grid_owning_t>);
  static_assert(concepts::grid<grid_n_owning_t>);
  static_assert(concepts::grid<grid_device_t>);

  // Fill the bin data for every test
  // bin test entries
  grid_owning_t::bin_container_type bin_data{};
  // 40 000 entries
  bin_data.entries.resize(80'000u);
  // 20 000 bins
  bin_data.bins.resize(40'000u);

  int i{0};
  dindex offset{0u};
  scalar entry{0.f};
  complete<> completer{};

  // Test data to compare bin content against
  std::vector<scalar> seq;
  seq.reserve(80'000);

  for (auto& data : bin_data.bins) {
    data.offset = offset;
    // Every second bin holds one element, otherwise three
    data.capacity = (i % 2) ? 1u : 3u;

    detray::bins::dynamic_array bin{bin_data.entries.data(), data};

    ASSERT_TRUE(bin.capacity() == (i % 2 ? 1u : 3u));
    ASSERT_TRUE(bin.size() == 0);

    offset += bin.capacity();

    // Populate the bin
    completer(bin, entry);

    for (auto e : bin) {
      ASSERT_TRUE(e == entry);
      seq.push_back(e);
    }
    entry += 1.f;
    ++i;
  }

  // Copy data that will be moved into the data owning types
  dvector<scalar> bin_edges_cp(bin_edges);
  dvector<dsized_index_range> edge_ranges_cp(edge_ranges);
  grid_owning_t::bin_container_type bin_data_cp(bin_data);

  // Data-owning axes and grid
  cartesian_3D<is_owning, host_container_types> axes_own(
      std::move(edge_ranges_cp), std::move(bin_edges_cp));
  grid_owning_t grid_own(std::move(bin_data_cp), std::move(axes_own));

  // Check a few basics
  EXPECT_EQ(grid_own.dim, 3u);
  EXPECT_EQ(grid_own.nbins(), 40'000u);
  auto y_axis = grid_own.get_axis<label::e_y>();
  EXPECT_EQ(y_axis.nbins(), 40u);
  auto z_axis =
      grid_own.get_axis<single_axis<closed<label::e_z>, regular<scalar>>>();
  EXPECT_EQ(z_axis.nbins(), 50u);

  // Check equality operator:
  // - Copy a second time for the comparison
  dvector<scalar> bin_edges_cp2(bin_edges);
  dvector<dsized_index_range> edge_ranges_cp2(edge_ranges);
  grid_owning_t::bin_container_type bin_data_cp2(bin_data);

  // Make a second grid
  cartesian_3D<is_owning, host_container_types> axes_own2(
      std::move(edge_ranges_cp2), std::move(bin_edges_cp2));

  grid_owning_t grid_own2(std::move(bin_data_cp2), std::move(axes_own2));

  // CHECK equality
  EXPECT_TRUE(grid_own == grid_own2);

  // Create non-owning grid
  grid_n_owning_t grid_n_own(&bin_data, ax_n_own);

  // Test for consistency with owning grid
  EXPECT_EQ(grid_n_own.dim, grid_own.dim);
  y_axis = grid_n_own.get_axis<label::e_y>();
  EXPECT_EQ(y_axis.nbins(), grid_own.get_axis<label::e_y>().nbins());
  z_axis = grid_n_own.get_axis<label::e_z>();
  EXPECT_EQ(z_axis.nbins(), grid_own.get_axis<label::e_z>().nbins());

  // Construct a grid from a view
  grid_owning_t::view_type grid_view = get_data(grid_own);
  grid_device_t device_grid(grid_view);

  // Test for consistency with non-owning grid
  EXPECT_EQ(device_grid.dim, grid_n_own.dim);
  auto y_axis_dev = device_grid.get_axis<label::e_y>();
  EXPECT_EQ(y_axis_dev.nbins(), grid_n_own.get_axis<label::e_y>().nbins());
  auto z_axis_dev = device_grid.get_axis<label::e_z>();
  EXPECT_EQ(z_axis_dev.nbins(), grid_n_own.get_axis<label::e_z>().nbins());

  // Test the global bin iteration
  auto flat_bin_view = grid_own.all();

  static_assert(detray::ranges::bidirectional_range<decltype(flat_bin_view)>);
  // TODO: Const-correctness issue
  // static_assert(detray::ranges::random_access_range<typename decltype(
  //                  flat_bin_view)>);

  EXPECT_EQ(seq.size(), 80'000u);
  EXPECT_EQ(flat_bin_view.size(), 80'000u);
  EXPECT_TRUE(
      std::equal(flat_bin_view.begin(), flat_bin_view.end(), seq.begin()));

  // Test the global bin iteration: non-owning grid
  auto flat_bin_view2 = grid_n_own.all();

  static_assert(detray::ranges::bidirectional_range<decltype(flat_bin_view2)>);

  EXPECT_EQ(seq.size(), 80'000u);
  EXPECT_EQ(flat_bin_view2.size(), 80'000u);
  EXPECT_TRUE(
      std::equal(flat_bin_view2.begin(), flat_bin_view2.end(), seq.begin()));
}

/// Integration test: Test replace population
GTEST_TEST(detray_grid, replace_population) {
  // Non-owning, 3D cartesian  grid
  using grid_t = grid<test_algebra, decltype(ax_n_own), bins::single<scalar>,
                      simple_serializer, host_container_types, false>;

  static_assert(concepts::grid<grid_t>);

  // init
  using bin_t = grid_t::bin_type;
  grid_t::bin_container_type bin_data{};
  bin_data.resize(40'000u, bin_t{});

  // Create non-owning grid
  grid_t g3r(&bin_data, ax_n_own);

  // Test the initialization
  point3 p = {-10.f, -20.f, 0.f};
  for (int ib0 = 0; ib0 < 20; ++ib0) {
    for (int ib1 = 0; ib1 < 40; ++ib1) {
      for (int ib2 = 0; ib2 < 100; ib2 += 2) {
        p = {static_cast<scalar>(-10 + ib0), static_cast<scalar>(-20 + ib1),
             static_cast<scalar>(0 + ib2)};
        EXPECT_NEAR(g3r.bin(p)[0], std::numeric_limits<scalar>::max(), tol);
      }
    }
  }

  p = {-4.5f, -4.5f, 4.5f};
  // Fill and read
  g3r.template populate<replace<>>(p, 3.f);
  EXPECT_NEAR(g3r.bin(p)[0], static_cast<scalar>(3u), tol);

  // Fill and read two times, fill first 0-99, then 100-199
  for (unsigned int il = 0u; il < 2u; ++il) {
    scalar counter{static_cast<scalar>(il * 100u)};
    for (int ib0 = 0; ib0 < 20; ++ib0) {
      for (int ib1 = 0; ib1 < 40; ++ib1) {
        for (int ib2 = 0; ib2 < 100; ib2 += 2) {
          p = {static_cast<scalar>(-10 + ib0), static_cast<scalar>(-20 + ib1),
               static_cast<scalar>(0 + ib2)};
          g3r.template populate<replace<>>(p, counter);
          EXPECT_NEAR(g3r.bin(p)[0], counter, tol);
          counter += 1.f;
        }
      }
    }
  }
}

/// Test bin entry retrieval
GTEST_TEST(detray_grid, complete_population) {
  // Non-owning, 3D cartesian, completing grid (4 dims and sort)
  using grid_t =
      grid<test_algebra, decltype(ax_n_own), bins::static_array<scalar, 4>,
           simple_serializer, host_container_types, false>;
  using bin_t = grid_t::bin_type;
  using bin_content_t = std::array<scalar, 4>;

  static_assert(concepts::grid<grid_t>);

  // init
  grid_t::bin_container_type bin_data{};
  bin_data.resize(40'000u, bin_t{});
  // Create non-owning grid
  grid_t g3c(&bin_data, ax_n_own);

  // Test the initialization
  point3 p = {-10.f, -20.f, 0.f};
  bin_t invalid{};
  for (int ib0 = 0; ib0 < 20; ++ib0) {
    for (int ib1 = 0; ib1 < 40; ++ib1) {
      for (int ib2 = 0; ib2 < 100; ib2 += 2) {
        p = {static_cast<scalar>(-10 + ib0), static_cast<scalar>(-20 + ib1),
             static_cast<scalar>(0 + ib2)};
        test_content(g3c, p, invalid);
      }
    }
  }

  p = {-4.5f, -4.5f, 4.5f};
  bin_content_t expected{4.f, 4.f, 4.f, 4.f};
  // Fill and read
  g3c.template populate<complete<>>(p, 4.f);
  test_content(g3c, p, expected);

  // Test search without search window
  for (scalar entry : g3c.bin(p)) {
    EXPECT_EQ(entry, 4.f);
  }
}

/// Test bin entry retrieval
GTEST_TEST(detray_grid, regular_attach_population) {
  // Non-owning, 3D cartesian, completing grid (4 dims and sort)
  using grid_t =
      grid<test_algebra, decltype(ax_n_own), bins::static_array<scalar, 4>,
           simple_serializer, host_container_types, false>;
  using bin_t = grid_t::bin_type;
  using bin_content_t = std::array<scalar, 4>;

  static_assert(concepts::grid<grid_t>);

  // init
  grid_t::bin_container_type bin_data{};
  bin_data.resize(40'000u, bin_t{});

  // Create non-owning grid
  grid_t g3ra(&bin_data, ax_n_own);

  // Test the initialization
  point3 p = {-10.f, -20.f, 0.f};
  bin_t invalid{};
  for (int ib0 = 0; ib0 < 20; ++ib0) {
    for (int ib1 = 0; ib1 < 40; ++ib1) {
      for (int ib2 = 0; ib2 < 100; ib2 += 2) {
        p = {static_cast<scalar>(-10 + ib0), static_cast<scalar>(-20 + ib1),
             static_cast<scalar>(0 + ib2)};
        test_content(g3ra, p, invalid);
      }
    }
  }

  p = {-4.5f, -4.5f, 4.5f};
  bin_content_t expected{5.f, inf, inf, inf};
  // Fill and read
  g3ra.template populate<attach<>>(p, 5.f);
  test_content(g3ra, p, expected);

  // Test search without search window
  for (scalar entry : g3ra.bin(p)) {
    EXPECT_EQ(entry, 5.f);
  }
}
