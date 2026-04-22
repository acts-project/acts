// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Detray test include(s)
#include <detray/test/framework/assert.hpp>

#include "sf_finders_grid_cuda_kernel.hpp"

// VecMem include(s).
#include <vecmem/memory/cuda/managed_memory_resource.hpp>

// GTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <limits>

using namespace detray;

namespace {

/// Test single entry bin content element by element
template <typename bin_content_t, typename content_t>
void test_content(const bin_content_t& bin_content, const content_t& expected) {
  unsigned int i{0u};
  for (const auto& elem : bin_content) {
    ASSERT_NEAR(elem, expected[i++], tol);
  }
}

/// Test collection entries in a bin element by element
template <typename content_t, typename expected_content_t>
void test_entry_collection(const content_t& bin_content,
                           const expected_content_t& expected) {
  // Running indices for the points in the collection and the elements of a
  // single point3
  unsigned int i = 0u;
  unsigned int j = 0u;
  for (const auto& entry : bin_content) {
    const auto& expt_entry = expected[i++];
    j = 0u;
    for (const auto& elem : entry) {
      ASSERT_NEAR(elem, expt_entry[j++], tol);
    }
  }
}

}  // anonymous namespace

TEST(grids_cuda, grid3_replace_populator) {
  // memory resource
  vecmem::cuda::managed_memory_resource mng_mr;

  // Build multi-axis
  using axes_t = host_grid3_single::axes_type;
  using bin_t = host_grid3_single::bin_type;

  typename axes_t::edge_offset_container_type axis_data(&mng_mr);
  typename axes_t::edges_container_type bin_edges(&mng_mr);

  axis_data.reserve(3);
  axis_data.insert(axis_data.begin(),
                   {dsized_index_range{0u, 3u}, dsized_index_range{2u, 6u},
                    dsized_index_range{4u, 10u}});
  bin_edges.reserve(6);
  bin_edges.insert(bin_edges.begin(), {-1.f, 2.f, 0.f, 6.f, -5.f, 5.f});

  axes_t axes(std::move(axis_data), std::move(bin_edges));

  // build host grid
  host_grid3_single::bin_container_type bin_data(&mng_mr);
  bin_data.resize(3u * 6u * 10u, bin_t{});

  host_grid3_single g3(std::move(bin_data), std::move(axes));

  const auto& axis_x = g3.template get_axis<axis::label::e_x>();
  const auto& axis_y = g3.template get_axis<axis::label::e_y>();
  const auto& axis_z = g3.template get_axis<axis::label::e_z>();

  // pre-check
  for (unsigned int i_x = 0u; i_x < axis_x.nbins(); i_x++) {
    for (unsigned int i_y = 0u; i_y < axis_y.nbins(); i_y++) {
      for (unsigned int i_z = 0u; i_z < axis_z.nbins(); i_z++) {
        const auto& bin = g3.bin(i_x, i_y, i_z);
        auto invalid_bin = bin_t{};
        test_content(bin.value(), invalid_bin.value());
      }
    }
  }
  // run device side test, which populates the grid
  grid_replace_test(get_data(g3), axis_x.nbins(), axis_y.nbins(),
                    axis_z.nbins());
  // post-check
  for (unsigned int i_x = 0u; i_x < axis_x.nbins(); i_x++) {
    for (unsigned int i_y = 0u; i_y < axis_y.nbins(); i_y++) {
      for (unsigned int i_z = 0u; i_z < axis_z.nbins(); i_z++) {
        const dindex gbin = g3.serialize({i_x, i_y, i_z});

        const auto& bin = g3.bin(gbin);

        const scalar gbin_f{static_cast<scalar>(gbin)};
        const point3 tp{axis_x.min() + gbin_f * axis_x.bin_width(),
                        axis_y.min() + gbin_f * axis_y.bin_width(),
                        axis_z.min() + gbin_f * axis_z.bin_width()};

        test_content(bin.value(), tp);
      }
    }
  }

  // Print the grid on device
  // print_grid<device_grid3_single>(get_data(g3), axis_x.nbins(),
  // axis_y.nbins(), axis_z.nbins());
}

TEST(grids_cuda, grid2_replace_populator_ci) {
  vecmem::cuda::managed_memory_resource mng_mr;

  // Build multi-axis
  using axes_t = host_grid2_single_ci::axes_type;
  using bin_t = host_grid2_single_ci::bin_type;

  typename axes_t::edge_offset_container_type axis_data(&mng_mr);
  typename axes_t::edges_container_type bin_edges(&mng_mr);

  axis_data.reserve(2u);
  axis_data.insert(axis_data.end(),
                   {dsized_index_range{0u, 4u}, dsized_index_range{5u, 2u}});
  bin_edges.reserve(7u);
  bin_edges.insert(bin_edges.end(), {1.f, 3.f, 9.f, 27.f, 81.f, -2.f, 4.f});

  axes_t axes(std::move(axis_data), std::move(bin_edges));

  // build host grid
  host_grid2_single_ci::bin_container_type bin_data(&mng_mr);
  bin_data.resize(4u * 2u, bin_t{});

  host_grid2_single_ci g2(std::move(bin_data), std::move(axes));

  const auto& axis_r = g2.template get_axis<axis::label::e_r>();
  const auto& axis_phi = g2.template get_axis<axis::label::e_phi>();

  // pre-check
  for (unsigned int i_r = 0u; i_r < axis_r.nbins(); i_r++) {
    for (unsigned int i_phi = 0u; i_phi < axis_phi.nbins(); i_phi++) {
      const auto& bin = g2.bin(i_r, i_phi);
      auto invalid_bin = bin_t{};

      test_content(bin.value(), invalid_bin.value());
    }
  }

  // run device side test, which populates the grid
  grid_replace_ci_test(get_data(g2), axis_r.nbins(), axis_phi.nbins());

  // post-check
  for (unsigned int i_r = 0u; i_r < axis_r.nbins(); i_r++) {
    for (unsigned int i_phi = 0u; i_phi < axis_phi.nbins(); i_phi++) {
      const dindex gbin = g2.serialize({i_r, i_phi});
      const auto& bin = g2.bin(gbin);

      const scalar gbin_f{static_cast<scalar>(gbin)};
      const point3 tp{axis_r.min() + gbin_f * axis_r.bin_width(i_r),
                      axis_phi.min() + gbin_f * axis_phi.bin_width(), 0.5f};

      test_content(bin.value(), tp);
    }
  }

  // Print the grid on device
  // print_grid<device_grid2_single_ci>(get_data(g2), axis_r.nbins(),
  // axis_phi.nbins());
}

TEST(grids_cuda, grid2_complete_populator) {
  // memory resource
  vecmem::cuda::managed_memory_resource mng_mr;

  // Build multi-axis
  using axes_t = host_grid2_array::axes_type;
  using bin_t = host_grid2_array::bin_type;

  typename axes_t::edge_offset_container_type axis_data(&mng_mr);
  typename axes_t::edges_container_type bin_edges(&mng_mr);

  axis_data.reserve(2u);
  axis_data.insert(axis_data.begin(),
                   {dsized_index_range{0u, 3u}, dsized_index_range{2u, 7u}});
  bin_edges.reserve(4);
  bin_edges.insert(bin_edges.begin(), {0.f, 3.f, -1.f, 6.f});

  axes_t axes(std::move(axis_data), std::move(bin_edges));

  // build host grid
  const point3 first_tp{3.f, 3.f, 3.f};

  host_grid2_array::bin_container_type bin_data(&mng_mr);
  bin_data.resize(3u * 7u, bin_t{}.init(first_tp));

  host_grid2_array g2(std::move(bin_data), std::move(axes));

  const auto& axis_r = g2.template get_axis<axis::label::e_r>();
  const auto& axis_phi = g2.template get_axis<axis::label::e_phi>();

  auto width_r = axis_r.bin_width();
  auto width_phi = axis_phi.bin_width();

  // pre-check
  for (unsigned int i_r = 0u; i_r < axis_r.nbins(); i_r++) {
    for (unsigned int i_phi = 0u; i_phi < axis_phi.nbins(); i_phi++) {
      const auto& bin = g2.bin(i_r, i_phi);
      auto invalid_bin = bin_t{}.init(first_tp);

      test_entry_collection(bin, invalid_bin);
    }
  }

  // run device side test, which populates the grid
  grid_complete_test(get_data(g2), axis_r.nbins(), axis_phi.nbins());

  // post-check
  for (unsigned int i_r = 0u; i_r < axis_r.nbins(); i_r++) {
    for (unsigned int i_phi = 0u; i_phi < axis_phi.nbins(); i_phi++) {
      const dindex gbin = g2.serialize({i_r, i_phi});
      const auto& bin = g2.bin(gbin);

      // Other point with which the bin has been completed
      const scalar gbin_f{static_cast<scalar>(gbin)};
      const point3 tp{axis_r.min() + gbin_f * width_r,
                      axis_phi.min() + gbin_f * width_phi,
                      static_cast<scalar>(0.5)};

      // Go through all points and compare
      int pt_idx{0};
      for (const auto& pt : bin) {
        if (pt_idx == 0) {
          EXPECT_EQ(pt, first_tp);
        } else {
          EXPECT_EQ(pt, tp);
        }
        pt_idx++;
      }
    }
  }

  // Print the grid on device
  // print_grid<device_grid2_array>(get_data(g2), axis_r.nbins(),
  // axis_phi.nbins());
}

// Test both the attach population and the const grid reading
TEST(grids_cuda, grid2_attach_populator) {
  // memory resource
  vecmem::cuda::managed_memory_resource mng_mr;

  // Build multi-axis
  using axes_t = host_grid2_array::axes_type;
  using bin_t = host_grid2_array::bin_type;

  typename axes_t::edge_offset_container_type axis_data(&mng_mr);
  typename axes_t::edges_container_type bin_edges(&mng_mr);

  axis_data.reserve(2);
  axis_data.insert(axis_data.begin(),
                   {dsized_index_range{0u, 2u}, dsized_index_range{2u, 65u}});
  bin_edges.reserve(4);
  bin_edges.insert(bin_edges.begin(),
                   {0.f, 6.f, -constant<scalar>::pi, constant<scalar>::pi});

  axes_t axes(std::move(axis_data), std::move(bin_edges));

  // build host grid
  const point3 first_tp{3.f, 3.f, 3.f};
  const point3 invalid_tp{0.f, 0.f, 0.f};

  host_grid2_array::bin_container_type bin_data(&mng_mr);
  bin_data.resize(2u * 65u, bin_t{}.init(first_tp));

  host_grid2_array g2(std::move(bin_data), std::move(axes));

  const auto& axis_r = g2.template get_axis<axis::label::e_r>();
  const auto& axis_phi = g2.template get_axis<axis::label::e_phi>();

  auto width_r = axis_r.bin_width();
  auto width_phi = axis_phi.bin_width();

  // pre-check
  for (unsigned int i_r = 0u; i_r < axis_r.nbins(); i_r++) {
    for (unsigned int i_phi = 0u; i_phi < axis_phi.nbins(); i_phi++) {
      auto bin = g2.bin(i_r, i_phi);
      auto invalid_bin = bin_t{}.init(first_tp);

      test_entry_collection(bin, invalid_bin);
    }
  }

  // run device side test, which populates the grid
  grid_attach_test(get_data(g2), axis_r.nbins(), axis_phi.nbins());

  // post-check
  for (unsigned int i_r = 0u; i_r < axis_r.nbins(); i_r++) {
    for (unsigned int i_phi = 0u; i_phi < axis_phi.nbins(); i_phi++) {
      const dindex gbin = g2.serialize({i_r, i_phi});
      const auto& bin = g2.bin(gbin);

      // Other point with which the bin has been completed
      const scalar gbin_f{static_cast<scalar>(gbin)};
      const point3 tp{axis_r.min() + gbin_f * width_r,
                      axis_phi.min() + gbin_f * width_phi,
                      static_cast<scalar>(0.5)};

      // Go through all points and compare
      int pt_idx{0};
      for (const auto& pt : bin) {
        if (pt_idx == 0) {
          EXPECT_POINT3_NEAR(pt, first_tp, 1e-6);
        } else if (pt_idx == 1) {
          EXPECT_POINT3_NEAR(pt, tp, 1e-6);
        } else {
          EXPECT_POINT3_NEAR(pt, invalid_tp, 1e-6);
        }
        pt_idx++;
      }
    }
  }

  // Print the grid on device
  // print_grid<device_grid2_array>(get_data(g2), axis_r.nbins(),
  // axis_phi.nbins());
}

// Test the attach population on a grid with dynamic bin capacities
TEST(grids_cuda, grid2_dynamic_attach_populator) {
  // memory resource
  vecmem::cuda::managed_memory_resource mng_mr;

  // Build multi-axis
  using axes_t = host_grid2_dynamic_array::axes_type;
  using bin_t = host_grid2_dynamic_array::bin_type;

  typename axes_t::edge_offset_container_type axis_data(&mng_mr);
  typename axes_t::edges_container_type bin_edges(&mng_mr);

  axis_data.reserve(2);
  axis_data.insert(axis_data.begin(),
                   {dsized_index_range{0u, 2u}, dsized_index_range{2u, 65u}});
  bin_edges.reserve(4);
  bin_edges.insert(bin_edges.begin(),
                   {0.f, 6.f, -constant<scalar>::pi, constant<scalar>::pi});

  axes_t axes(std::move(axis_data), std::move(bin_edges));

  // build host grid
  const point3 first_tp{3.f, 3.f, 3.f};
  const point3 invalid_tp{0.f, 0.f, 0.f};

  host_grid2_dynamic_array::bin_container_type bin_data{&mng_mr};
  vecmem::vector<bin_t::entry_type> entries{&mng_mr};
  bin_data.bins.resize(2 * 65);
  bin_data.entries.resize(4 * bin_data.bins.size());

  int i{0};
  dindex offset{0u};
  attach<> attacher{};

  // bin content with different bin capacities
  for (auto& data : bin_data.bins) {
    data.offset = offset;
    // Every second bin holds one element, otherwise three
    data.capacity = (i % 2) ? 1u : 3u;

    detray::bins::dynamic_array bin{bin_data.entries.data(), data};

    ASSERT_TRUE(bin.capacity() == (i % 2 ? 1u : 3u));
    ASSERT_TRUE(bin.size() == 0);

    offset += bin.capacity();

    // Populate the bin
    attacher(bin, first_tp);
    ++i;
  }

  host_grid2_dynamic_array g2(std::move(bin_data), std::move(axes));

  const auto& axis_r = g2.template get_axis<axis::label::e_r>();
  const auto& axis_phi = g2.template get_axis<axis::label::e_phi>();

  auto width_r = axis_r.bin_width();
  auto width_phi = axis_phi.bin_width();

  // pre-check
  for (unsigned int i_phi = 0u; i_phi < axis_phi.nbins(); i_phi++) {
    for (unsigned int i_r = 0u; i_r < axis_r.nbins(); i_r++) {
      int pt_idx{0};
      for (auto e : g2.bin(i_r, i_phi)) {
        if (pt_idx == 0) {
          EXPECT_EQ(e, first_tp);
        } else {
          EXPECT_EQ(e, invalid_tp);
        }
      }
    }
  }

  // run device side test, which populates the grid
  grid_dynamic_attach_test(get_data(g2), axis_r.nbins(), axis_phi.nbins());

  // post-check
  for (unsigned int i_r = 0u; i_r < axis_r.nbins(); i_r++) {
    for (unsigned int i_phi = 0u; i_phi < axis_phi.nbins(); i_phi++) {
      const dindex gbin = g2.serialize({i_r, i_phi});
      const auto& bin = g2.bin(gbin);

      // Other point with which the bin has been completed
      const scalar gbin_f{static_cast<scalar>(gbin)};
      const point3 tp{axis_r.min() + gbin_f * width_r,
                      axis_phi.min() + gbin_f * width_phi,
                      static_cast<scalar>(0.5)};

      // Go through all points and compare
      int pt_idx{0};
      for (const auto& e : bin) {
        if (pt_idx == 0) {
          EXPECT_POINT3_NEAR(e, first_tp, 1e-6);
        } else if (pt_idx == 1) {
          EXPECT_POINT3_NEAR(e, tp, 1e-6);
        } else {
          EXPECT_POINT3_NEAR(e, invalid_tp, 1e-6);
        }
        pt_idx++;
      }
    }
  }

  // Print the grid on device
  // print_grid<device_grid2_dynamic_array>(get_data(g2), axis_r.nbins(),
  // axis_phi.nbins());
}

TEST(grids_cuda, cylindrical3D_collection) {
  // Data-owning grid collection
  vecmem::cuda::managed_memory_resource mng_mr;

  using grid_collection_t = grid_collection<n_own_host_grid3_array>;
  using bin_t = grid_collection_t::value_type::bin_type;

  vecmem::vector<typename grid_collection_t::size_type> grid_offsets(&mng_mr);
  typename grid_collection_t::bin_container_type bin_data(&mng_mr);
  typename grid_collection_t::edge_offset_container_type edge_ranges(&mng_mr);
  typename grid_collection_t::edges_container_type bin_edges(&mng_mr);

  // Offsets for the grids into the bin storage
  grid_offsets.reserve(3u);
  grid_offsets.insert(grid_offsets.begin(), {0u, 48u, 72u});

  // Offsets into edges container and #bins for all axes
  edge_ranges.reserve(9u);
  edge_ranges.insert(edge_ranges.begin(),
                     {dsized_index_range{0u, 2u}, dsized_index_range{2u, 4u},
                      dsized_index_range{4u, 6u}, dsized_index_range{6u, 1u},
                      dsized_index_range{8u, 3u}, dsized_index_range{10u, 8u},
                      dsized_index_range{12u, 5u}, dsized_index_range{14u, 5u},
                      dsized_index_range{16u, 5u}});

  // Bin edges for all axes (two boundaries for regular binned axes)
  bin_edges.reserve(18u);
  bin_edges.insert(bin_edges.begin(),
                   {-10.f, 10.f, -20.f, 20.f, 0.f, 120.f, -5.f, 5.f, -15.f,
                    15.f, 0.f, 50.f, -15.f, 15.f, -35.f, 35.f, 0.f, 550.f});

  // Bin test entries
  bin_data.resize(197u, bin_t{});

  // Read the number of bins and bin entries on device
  vecmem::vector<unsigned int> n_bins(9u, &mng_mr);
  vecmem::vector<std::array<dindex, 3>> result_bins(bin_data.size(), &mng_mr);

  grid_collection_t grid_coll(std::move(grid_offsets), std::move(bin_data),
                              std::move(edge_ranges), std::move(bin_edges));

  // Call test function
  const auto& axis_r = grid_coll[2].template get_axis<axis::label::e_r>();
  const auto& axis_phi = grid_coll[2].template get_axis<axis::label::e_phi>();
  const auto& axis_z = grid_coll[2].template get_axis<axis::label::e_z>();

  grid_collection_test(get_data(grid_coll), vecmem::get_data(n_bins),
                       vecmem::get_data(result_bins), grid_coll.size(),
                       axis_r.nbins(), axis_phi.nbins(), axis_z.nbins());

  // Compare results
  EXPECT_EQ(4u, n_bins[0]);
  EXPECT_EQ(4u, n_bins[1]);
  EXPECT_EQ(8u, n_bins[2]);
  EXPECT_EQ(3u, n_bins[3]);
  EXPECT_EQ(3u, n_bins[4]);
  EXPECT_EQ(10u, n_bins[5]);
  EXPECT_EQ(7u, n_bins[6]);
  EXPECT_EQ(5u, n_bins[7]);
  EXPECT_EQ(7u, n_bins[8]);

  for (unsigned int i{0u}; i < bin_data.size(); ++i) {
    test_content(bin_data[i], result_bins[i]);
  }
}
