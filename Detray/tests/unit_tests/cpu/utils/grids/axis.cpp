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

// detray core
#include "detray/utils/grid/axis.hpp"

#include "detray/definitions/indexing.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/coordinates/coordinates.hpp"
#include "detray/utils/grid/concepts.hpp"
#include "detray/utils/grid/detail/axis_binning.hpp"
#include "detray/utils/grid/detail/axis_bounds.hpp"

// Detray test include(s)
#include "detray/test/framework/types.hpp"

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <limits>

using namespace detray;
using namespace detray::axis;

namespace {

using test_algebra = test::algebra;
using scalar = test::scalar;
using point3 = test::point3;

// Alias for testing
template <bool ownership, typename containers>
using cartesian_3D =
    multi_axis<ownership, cartesian3D<test_algebra>,
               single_axis<closed<label::e_x>, regular<scalar, containers>>,
               single_axis<closed<label::e_y>, regular<scalar, containers>>,
               single_axis<closed<label::e_z>, regular<scalar, containers>>>;

// Floating point comparison
constexpr scalar tol{1e-5f};

}  // anonymous namespace

GTEST_TEST(detray_grid, open_regular_axis) {
  // Lower bin edges: min and max bin edge for the regular axis
  vecmem::vector<scalar> bin_edges = {-10.f, -5.f, -3.f, 7.f, 7.f, 14.f, 20.f};
  // Regular axis: first entry is the offset in the bin edges vector (2), the
  // second entry is the number of bins (10): Lower and upper bin edges of
  // the max and min bin are -3 and 7 => stepsize is (7 - (-3)) / 10 = 1
  // ... -3, -2, -1,  0,  1,  2,  3,  4,  5,  6,  7 ...
  //   [0] [1] [2] [3] [4] [5] [6] [7] [8] [9] [10] [11]
  const dsized_index_range edge_range = {2u, 10u};

  // An open regular x-axis
  using x_axis_t = single_axis<open<label::e_x>, regular<scalar>>;
  x_axis_t or_axis{edge_range, &bin_edges};

  static_assert(concepts::axis<x_axis_t>);

  // Test axis bounds
  EXPECT_EQ(or_axis.label(), axis::label::e_x);
  EXPECT_EQ(or_axis.bounds(), axis::bounds::e_open);
  EXPECT_EQ(or_axis.binning(), axis::binning::e_regular);

  // N bins
  EXPECT_EQ(or_axis.nbins(), 10u + 2u);
  EXPECT_EQ(or_axis.m_binning.nbins(), 10u);
  EXPECT_NEAR(or_axis.span()[0], -3.f, tol);
  EXPECT_NEAR(or_axis.span()[1], 7.f, tol);
  // Axis bin access
  // Axis is open: Every value smaller than -3 is mapped into bin 0:
  // which is (-inf, -3]
  EXPECT_EQ(or_axis.bin(-4.f), 0u);
  EXPECT_EQ(or_axis.bin(2.5f), 6u);
  EXPECT_EQ(or_axis.bin(1.5f), 5u);
  EXPECT_EQ(or_axis.bin(7.5f), 11u);
  EXPECT_EQ(or_axis.bin(-3.5f), 0u);
  // Axis is open: Every value greater than 7 is mapped into bin 11:
  // which is [7, inf)
  EXPECT_EQ(or_axis.bin(8.f), 11u);

  // Axis range access - binned (symmetric & asymmetric)
  const darray<dindex, 2> nhood00i = {0u, 0u};
  const darray<dindex, 2> nhood01i = {0u, 1u};
  const darray<dindex, 2> nhood11i = {1u, 1u};
  const darray<dindex, 2> nhood44i = {4u, 4u};
  const darray<dindex, 2> nhood55i = {5u, 5u};

  axis::bin_range expected_range = {6u, 7u};
  EXPECT_EQ(or_axis.range(2.5f, nhood00i), expected_range);
  expected_range = {6u, 8u};
  EXPECT_EQ(or_axis.range(2.5f, nhood01i), expected_range);
  expected_range = {5u, 8u};
  EXPECT_EQ(or_axis.range(2.5f, nhood11i), expected_range);
  expected_range = {2u, 11u};
  EXPECT_EQ(or_axis.range(2.5f, nhood44i), expected_range);
  expected_range = {1u, 12u};
  EXPECT_EQ(or_axis.range(2.5f, nhood55i), expected_range);
  expected_range = {1u, 10u};
  EXPECT_EQ(or_axis.range(1.5f, nhood44i), expected_range);
  expected_range = {4u, 12u};
  EXPECT_EQ(or_axis.range(5.5f, nhood55i), expected_range);
  expected_range = {11u, 12u};
  EXPECT_EQ(or_axis.range(7.5f, nhood00i), expected_range);
  expected_range = {0u, 1u};
  EXPECT_EQ(or_axis.range(-3.5f, nhood00i), expected_range);

  // Axis range access - scalar (symmteric & asymmetric)
  const darray<scalar, 2> nhood00s = {0.f, 0.f};
  const darray<scalar, 2> nhood_tol = {0.01f, 0.01f};
  const darray<scalar, 2> nhood11s = {1.f, 1.f};
  const darray<scalar, 2> nhoodAlls = {10.f, 10.f};

  expected_range = {6u, 7u};
  EXPECT_EQ(or_axis.range(2.5f, nhood00s), expected_range);
  EXPECT_EQ(or_axis.range(2.5f, nhood_tol), expected_range);
  expected_range = {5u, 8u};
  EXPECT_EQ(or_axis.range(2.5f, nhood11s), expected_range);
  expected_range = {0u, 12u};
  EXPECT_EQ(or_axis.range(2.5f, nhoodAlls), expected_range);
}

GTEST_TEST(detray_grid, closed_regular_axis) {
  // Lower bin edges: min and max bin edge for the regular axis
  vecmem::vector<scalar> bin_edges = {-10.f, -3.f, -3.f, 7.f, 7.f, 14.f};
  // Regular axis: first entry is the offset in the bin edges vector (2), the
  // second entry is the number of bins (10): Lower and upper bin edges of
  // the max and min bin are -3 and 7 => stepsize is (7 - (-3)) / 10 = 1
  // -3, -2, -1,  0,  1,  2,  3,  4,  5,  6,  7
  //   [0] [1] [2] [3] [4] [5] [6] [7] [8] [9]
  const dsized_index_range edge_range = {2u, 10u};

  // A closed regular r-axis
  using r_axis_t = single_axis<closed<label::e_r>, regular<scalar>>;
  r_axis_t cr_axis{edge_range, &bin_edges};

  static_assert(concepts::axis<r_axis_t>);

  // Test axis bounds
  EXPECT_EQ(cr_axis.label(), axis::label::e_r);
  EXPECT_EQ(cr_axis.bounds(), axis::bounds::e_closed);
  EXPECT_EQ(cr_axis.binning(), axis::binning::e_regular);

  // N bins
  EXPECT_EQ(cr_axis.nbins(), 10u);
  EXPECT_NEAR(cr_axis.span()[0], -3.f, tol);
  EXPECT_NEAR(cr_axis.span()[1], 7.f, tol);
  // Axis bin access
  // Axis is closed: Every value smaller than -3 is mapped into bin 0:
  // which is [-3, -2)
  EXPECT_EQ(cr_axis.bin(-4.f), 0u);
  EXPECT_EQ(cr_axis.bin(2.5f), 5u);
  // Axis is closed: Every value greater than 7 is mapped into bin 9:
  // which is (6, 7]
  EXPECT_EQ(cr_axis.bin(8.f), 9u);

  // Axis range access - binned (symmetric & asymmetric)
  const darray<dindex, 2> nhood00i = {0u, 0u};
  const darray<dindex, 2> nhood01i = {0u, 1u};
  const darray<dindex, 2> nhood11i = {1u, 1u};
  const darray<dindex, 2> nhood44i = {4u, 4u};
  const darray<dindex, 2> nhood55i = {5u, 5u};

  axis::bin_range expected_range = {5u, 6u};
  EXPECT_EQ(cr_axis.range(2.5f, nhood00i), expected_range);
  expected_range = {5u, 7u};
  EXPECT_EQ(cr_axis.range(2.5f, nhood01i), expected_range);
  expected_range = {4u, 7u};
  EXPECT_EQ(cr_axis.range(2.5f, nhood11i), expected_range);
  expected_range = {1u, 10u};
  EXPECT_EQ(cr_axis.range(2.5f, nhood44i), expected_range);
  expected_range = {0u, 10u};
  EXPECT_EQ(cr_axis.range(2.5f, nhood55i), expected_range);
  expected_range = {0u, 9u};
  EXPECT_EQ(cr_axis.range(1.5f, nhood44i), expected_range);
  expected_range = {3u, 10u};
  EXPECT_EQ(cr_axis.range(5.5f, nhood55i), expected_range);
  expected_range = {9u, 10u};
  EXPECT_EQ(cr_axis.range(7.5f, nhood00i), expected_range);
  expected_range = {0u, 1u};
  EXPECT_EQ(cr_axis.range(-3.5f, nhood00i), expected_range);

  // Axis range access - scalar (symmteric & asymmetric)
  const darray<scalar, 2> nhood00s = {0.f, 0.f};
  const darray<scalar, 2> nhood_tol = {0.01f, 0.01f};
  const darray<scalar, 2> nhood11s = {1.f, 1.f};
  const darray<scalar, 2> nhoodAlls = {10.f, 10.f};

  expected_range = {5u, 6u};
  EXPECT_EQ(cr_axis.range(2.5f, nhood00s), expected_range);
  EXPECT_EQ(cr_axis.range(2.5f, nhood_tol), expected_range);
  expected_range = {4u, 7u};
  EXPECT_EQ(cr_axis.range(2.5f, nhood11s), expected_range);
  expected_range = {0u, 10u};
  EXPECT_EQ(cr_axis.range(2.5f, nhoodAlls), expected_range);
}

GTEST_TEST(detray_grid, circular_regular_axis) {
  // Let's say 36 modules, but with 4 directly at 0, pi/2, pi, -pi2
  const scalar half_module{constant<scalar>::pi / 72.f};
  const scalar phi_min{-constant<scalar>::pi + half_module};
  const scalar phi_max{constant<scalar>::pi - half_module};

  // Lower bin edges: min and max bin edge for the regular axis
  vecmem::vector<scalar> bin_edges = {-10.f, phi_min, phi_max, 56.f};
  // Regular axis: first entry is the offset in the bin edges vector (1), the
  // second entry is the number of bins (36)
  const dsized_index_range edge_range = {1u, 36u};

  // A closed regular x-axis
  using axis_t = single_axis<circular<>, regular<scalar>>;
  axis_t cr_axis(edge_range, &bin_edges);

  static_assert(concepts::axis<axis_t>);

  // Test axis bounds
  EXPECT_EQ(cr_axis.label(), axis::label::e_phi);
  EXPECT_EQ(cr_axis.bounds(), axis::bounds::e_circular);
  EXPECT_EQ(cr_axis.binning(), axis::binning::e_regular);

  // N bins
  EXPECT_EQ(cr_axis.nbins(), 36u);
  // Axis bin access
  // overflow
  EXPECT_EQ(cr_axis.bin(phi_max + tol), 0u);
  // underflow
  EXPECT_EQ(cr_axis.bin(phi_min - tol), 35u);
  // middle of the axis
  EXPECT_EQ(cr_axis.bin(0), 18u);

  // Bin wrapping test
  typename single_axis<circular<>, regular<scalar>>::bounds_type circ_bounds{};
  EXPECT_EQ(circ_bounds.wrap(4, 36u), 4u);
  EXPECT_EQ(circ_bounds.wrap(0, 36u), 0u);
  EXPECT_EQ(circ_bounds.wrap(-1, 36u), 35u);
  EXPECT_EQ(circ_bounds.wrap(36, 36u), 0u);
  EXPECT_EQ(circ_bounds.wrap(40, 36u), 4u);

  // Axis range access - binned (symmetric & asymmetric)
  const darray<dindex, 2> nhood00i = {0u, 0u};
  const darray<dindex, 2> nhood01i = {0u, 1u};
  const darray<dindex, 2> nhood11i = {1u, 1u};
  const darray<dindex, 2> nhood22i = {2u, 2u};

  axis::bin_range expected_range = {0u, 1u};
  EXPECT_EQ(circ_bounds.wrap(
                cr_axis.range(constant<scalar>::pi + tol, nhood00i), 36u),
            expected_range);
  expected_range = {0u, 2u};
  EXPECT_EQ(circ_bounds.wrap(
                cr_axis.range(constant<scalar>::pi + tol, nhood01i), 36u),
            expected_range);
  expected_range = {35u, 2u};
  EXPECT_EQ(circ_bounds.wrap(
                cr_axis.range(constant<scalar>::pi + tol, nhood11i), 36u),
            expected_range);
  expected_range = {34u, 3u};
  EXPECT_EQ(circ_bounds.wrap(
                cr_axis.range(constant<scalar>::pi + tol, nhood22i), 36u),
            expected_range);

  // Axis range access - scalar (symmetric & asymmteric)
  const darray<scalar, 2> nhood00s = {0.f, 0.f};
  const darray<scalar, 2> nhood_tol = {0.5f * tol, 0.5f * tol};
  const scalar bin_step{cr_axis.bin_width()};
  const darray<scalar, 2> nhood22s = {2.f * bin_step, 2.f * bin_step};

  expected_range = {0u, 1u};
  EXPECT_EQ(circ_bounds.wrap(
                cr_axis.range(constant<scalar>::pi + tol, nhood00s), 36u),
            expected_range);
  EXPECT_EQ(circ_bounds.wrap(
                cr_axis.range(constant<scalar>::pi + tol, nhood_tol), 36u),
            expected_range);
  expected_range = {34u, 3u};
  EXPECT_EQ(circ_bounds.wrap(
                cr_axis.range(constant<scalar>::pi + tol, nhood22s), 36u),
            expected_range);
}

GTEST_TEST(detray_grid, closed_irregular_axis) {
  // Lower bin edges: all lower bin edges for irregular binning, plus the
  // final upper bin edge
  vecmem::vector<scalar> bin_edges = {-100.f, -3.f, 1.f,  2.f, 4.f,
                                      8.f,    12.f, 15.f, 18.f};
  // Index offset and number of bins for the bin edges [-3, 15]
  dsized_index_range edge_range = {1u, 6u};

  // A closed irregular z-axis
  single_axis<closed<label::e_z>, irregular<scalar>> cir_axis(edge_range,
                                                              &bin_edges);

  // Test axis bounds
  EXPECT_EQ(cir_axis.label(), axis::label::e_z);
  EXPECT_EQ(cir_axis.bounds(), axis::bounds::e_closed);
  EXPECT_EQ(cir_axis.binning(), axis::binning::e_irregular);

  // Axis bin access
  //
  // N bins
  EXPECT_EQ(cir_axis.nbins(), 6u);
  // Bin tests
  EXPECT_EQ(cir_axis.bin(-2), 0u);
  EXPECT_EQ(cir_axis.bin(10), 4u);
  EXPECT_EQ(cir_axis.bin(5.8f), 3u);
  // Underflow test
  EXPECT_EQ(cir_axis.bin(-4), 0u);
  // Overflow test
  EXPECT_EQ(cir_axis.bin(17), 5u);

  // Axis range access - binned  (symmetric & asymmetric)
  const darray<dindex, 2> nhood00i = {0u, 0u};
  const darray<dindex, 2> nhood01i = {0u, 1u};
  const darray<dindex, 2> nhood11i = {1u, 1u};
  const darray<dindex, 2> nhood22i = {2u, 2u};

  axis::bin_range expected_range = {2u, 3u};
  EXPECT_EQ(cir_axis.range(3.f, nhood00i), expected_range);
  expected_range = {1u, 4u};
  EXPECT_EQ(cir_axis.range(3.f, nhood11i), expected_range);
  expected_range = {2u, 4u};
  EXPECT_EQ(cir_axis.range(3.f, nhood01i), expected_range);
  expected_range = {0u, 2u};
  EXPECT_EQ(cir_axis.range(0.f, nhood11i), expected_range);
  expected_range = {2u, 6u};
  EXPECT_EQ(cir_axis.range(10.f, nhood22i), expected_range);

  // Axis range access - scalar
  const darray<scalar, 2> nhood00s = {0.f, 0.f};
  const darray<scalar, 2> nhood10s = {1.5f, 0.2f};
  const darray<scalar, 2> nhood11s = {4.f, 5.5f};

  expected_range = {2u, 3u};
  EXPECT_EQ(cir_axis.range(3.f, nhood00s), expected_range);
  expected_range = {1u, 3u};
  EXPECT_EQ(cir_axis.range(3.f, nhood10s), expected_range);
  expected_range = {0u, 5u};
  EXPECT_EQ(cir_axis.range(3.f, nhood11s), expected_range);
}

template <typename axis_type>
void test_axis(const axis_type& /*unused*/) {
  // Lower bin edges: min and max bin edge for the regular axis
  vecmem::vector<scalar> bin_edges = {-10.f, -5.f, -3.f,  7.f,   7.f,  14.f,
                                      20.f,  50.f, 100.f, 120.f, 150.f};

  vecmem::vector<scalar> different_bin_edges = {
      -10.f, -5.f, -4.f, 7.f, 7.f, 14.f, 20.f, 50.f, 100.f, 120.f, 150.f};

  // Regular axis: first entry is the offset in the bin edges vector (2),
  // the second entry is the number of bins (10): Lower and upper bin
  // edges of the max and min bin are -3 and 7 => stepsize is (7 - (-3)) /
  // 10 = 1
  // ... -3, -2, -1,  0,  1,  2,  3,  4,  5,  6,  7 ...
  //   [0] [1] [2] [3] [4] [5] [6] [7] [8] [9] [10] [11]
  const dsized_index_range edge_range = {2u, 10u};

  const dsized_index_range different_edge_range = {2u, 8u};

  axis_type t_axis{edge_range, &bin_edges};

  axis_type another_t_axis{edge_range, &bin_edges};

  axis_type differnt_t_axis{edge_range, &different_bin_edges};

  axis_type differnt_edges_range_t_axis{different_edge_range, &bin_edges};

  // Test if the axes are equal
  EXPECT_EQ(t_axis, another_t_axis);

  // Confirm that the axes are not equal
  EXPECT_NE(t_axis, differnt_t_axis);

  // Confirm that the axes are not equal
  EXPECT_NE(t_axis, differnt_edges_range_t_axis);
};

GTEST_TEST(detray_grid, axis_comparison) {
  // Should be sufficient to test the comparison operators
  using x_axis_op_r_t = single_axis<open<label::e_x>, regular<scalar>>;
  using x_axis_cl_r_t = single_axis<closed<label::e_x>, regular<scalar>>;
  using x_axis_ci_r_t = single_axis<circular<label::e_x>, regular<scalar>>;
  using x_axis_op_ir_t = single_axis<open<label::e_x>, irregular<scalar>>;
  using x_axis_cl_ir_t = single_axis<closed<label::e_x>, irregular<scalar>>;
  using x_axis_ci_orr_t = single_axis<circular<label::e_x>, irregular<scalar>>;

  std::tuple<x_axis_op_r_t, x_axis_cl_r_t, x_axis_ci_r_t, x_axis_op_ir_t,
             x_axis_cl_ir_t, x_axis_ci_orr_t>
      axes = {};

  // Test them all
  std::apply([&](auto&... axis) { (test_axis(axis), ...); }, axes);
}

GTEST_TEST(detray_grid, multi_axis) {
  // readable axis ownership definition
  bool constexpr is_owning = true;
  bool constexpr is_not_owning = false;

  // Lower bin edges for all axes
  vecmem::vector<scalar> bin_edges = {-10.f, 10.f, -20.f, 20.f, 0.f, 100.f};
  // Offsets into edges container and #bins for all axes
  vecmem::vector<dsized_index_range> edge_ranges = {
      {0u, 20u}, {2u, 40u}, {4u, 50u}};

  // Non-owning multi axis test
  cartesian_3D<is_not_owning, host_container_types> axes(edge_ranges,
                                                         bin_edges);

  EXPECT_EQ(axes.dim, 3u);

  // Get single axis objects
  auto x_axis = axes.get_axis<label::e_x>();
  EXPECT_EQ(x_axis.nbins(), 20u);
  auto y_axis = axes.get_axis<label::e_y>();
  EXPECT_EQ(y_axis.nbins(), 40u);
  auto z_axis = axes.get_axis<label::e_z>();
  EXPECT_EQ(z_axis.nbins(), 50u);

  auto ax = axes.get_axis<decltype(y_axis)>();
  EXPECT_EQ(ax.nbins(), 40u);

  // Test bin search
  point3 p3{0.f, 0.f, 0.f};  // origin
  dmulti_index<dindex, 3> expected_bins{10u, 20u, 0u};
  EXPECT_EQ(axes.bins(p3), expected_bins);
  p3 = {-5.5f, -3.2f, 24.1f};
  expected_bins = {4u, 16u, 12u};
  EXPECT_EQ(axes.bins(p3), expected_bins);
  p3 = {-1.f, 21.f, 12.f};
  expected_bins = {9u, 39u, 6u};
  EXPECT_EQ(axes.bins(p3), expected_bins);

  // Test bin range search
  p3 = {1.f, 1.f, 10.f};
  // Axis range access - binned  (symmetric & asymmetric)
  const darray<dindex, 2> nhood00i = {0u, 0u};
  const darray<dindex, 2> nhood01i = {0u, 1u};
  const darray<dindex, 2> nhood22i = {2u, 2u};

  axis::multi_bin_range<3> expected_ranges{};
  expected_ranges[0] = {11u, 12u};
  expected_ranges[1] = {21u, 22u};
  expected_ranges[2] = {5u, 6u};
  EXPECT_EQ(axes.bin_ranges(p3, nhood00i), expected_ranges);
  expected_ranges[0] = {11u, 13u};
  expected_ranges[1] = {21u, 23u};
  expected_ranges[2] = {5u, 7u};
  EXPECT_EQ(axes.bin_ranges(p3, nhood01i), expected_ranges);
  expected_ranges[0] = {9u, 14u};
  expected_ranges[1] = {19u, 24u};
  expected_ranges[2] = {3u, 8u};
  EXPECT_EQ(axes.bin_ranges(p3, nhood22i), expected_ranges);

  // Axis range access - scalar
  const darray<scalar, 2> nhood00s = {0.f, 0.f};
  const darray<scalar, 2> nhood10s = {1.5f, 0.2f};
  const darray<scalar, 2> nhood11s = {4.f, 5.5f};

  expected_ranges[0] = {11u, 12u};
  expected_ranges[1] = {21u, 22u};
  expected_ranges[2] = {5u, 6u};
  EXPECT_EQ(axes.bin_ranges(p3, nhood00s), expected_ranges);
  expected_ranges[0] = {9u, 12u};
  expected_ranges[1] = {19u, 22u};
  expected_ranges[2] = {4u, 6u};
  EXPECT_EQ(axes.bin_ranges(p3, nhood10s), expected_ranges);
  expected_ranges[0] = {7u, 17u};
  expected_ranges[1] = {17u, 27u};
  expected_ranges[2] = {3u, 8u};
  EXPECT_EQ(axes.bin_ranges(p3, nhood11s), expected_ranges);

  // Owning multi axis test
  vecmem::vector<scalar> bin_edges_cp(bin_edges);
  // Offsets into edges container and #bins for all axes
  vecmem::vector<dsized_index_range> edge_ranges_cp(edge_ranges);

  cartesian_3D<is_owning, host_container_types> axes_own(
      std::move(edge_ranges_cp), std::move(bin_edges_cp));

  EXPECT_EQ(axes_own.bins(p3), axes.bins(p3));
  EXPECT_EQ(axes_own.bin_ranges(p3, nhood00i), axes.bin_ranges(p3, nhood00i));
  EXPECT_EQ(axes_own.bin_ranges(p3, nhood01i), axes.bin_ranges(p3, nhood01i));
  EXPECT_EQ(axes_own.bin_ranges(p3, nhood22i), axes.bin_ranges(p3, nhood22i));
  EXPECT_EQ(axes_own.bin_ranges(p3, nhood00s), axes.bin_ranges(p3, nhood00s));
  EXPECT_EQ(axes_own.bin_ranges(p3, nhood10s), axes.bin_ranges(p3, nhood10s));
  EXPECT_EQ(axes_own.bin_ranges(p3, nhood11s), axes.bin_ranges(p3, nhood11s));

  // Transfer to an owning multi-axis, which uses vecmem::device_vector
  cartesian_3D<is_owning, host_container_types>::view_type axes_view =
      get_data(axes_own);
  cartesian_3D<is_owning, device_container_types> axes_device(axes_view);

  // Get single axis objects
  auto x_axis_device = axes_device.get_axis<label::e_x>();
  EXPECT_EQ(x_axis_device.nbins(), 20u);
  auto y_axis_device = axes_device.get_axis<label::e_y>();
  EXPECT_EQ(y_axis_device.nbins(), 40u);
  auto z_axis_device = axes_device.get_axis<label::e_z>();
  EXPECT_EQ(z_axis_device.nbins(), 50u);
}

GTEST_TEST(detray_grid, multi_axis_comparison) {
  // Non-ownging test
  bool constexpr is_not_owning = false;

  // Lower bin edges for all axes
  vecmem::vector<scalar> bin_edges = {-10.f, 10.f, -20.f, 20.f, 0.f, 100.f};

  vecmem::vector<scalar> diff_bin_edges = {-10.f, 11.f, -20.f,
                                           20.f,  0.f,  100.f};

  // Offsets into edges container and #bins for all axes
  vecmem::vector<dsized_index_range> edge_ranges = {
      {0u, 20u}, {2u, 40u}, {4u, 50u}};

  // Offsets into edges container and #bins for all axes
  vecmem::vector<dsized_index_range> diff_edge_ranges = {
      {0u, 20u}, {3u, 40u}, {4u, 50u}};

  // Non-owning multi axis test : reference
  cartesian_3D<is_not_owning, host_container_types> ref_no_axes(edge_ranges,
                                                                bin_edges);

  // Non-owning multi axis test : test
  cartesian_3D<is_not_owning, host_container_types> test_no_axes(edge_ranges,
                                                                 bin_edges);

  // Owning multi axis test : reference
  EXPECT_EQ(ref_no_axes, test_no_axes);

  // Non-owning multi axis test : different bins
  cartesian_3D<is_not_owning, host_container_types> diff_edges_no_axes(
      edge_ranges, diff_bin_edges);
  // Confirm that the axes are not equal
  EXPECT_NE(ref_no_axes, diff_edges_no_axes);

  // Non-owning multi axis test : different bin ranges
  cartesian_3D<is_not_owning, host_container_types> diff_edge_ranges_no_axes(
      diff_edge_ranges, bin_edges);
  // Confirm that the axes are not equal
  EXPECT_NE(ref_no_axes, diff_edge_ranges_no_axes);

  // Owning test
  bool constexpr is_owning = true;

  vecmem::vector<scalar> bin_ref_o_edges = {-10.f, 10.f, -20.f,
                                            20.f,  0.f,  100.f};

  vecmem::vector<scalar> bin_same_o_edges = {-10.f, 10.f, -20.f,
                                             20.f,  0.f,  100.f};

  vecmem::vector<dsized_index_range> edge_ref_o_ranges = {
      {0u, 20u}, {2u, 40u}, {4u, 50u}};

  vecmem::vector<dsized_index_range> edge_diff_o_ranges = {
      {0u, 20u}, {2u, 39u}, {4u, 50u}};

  vecmem::vector<dsized_index_range> edge_same_o_ranges = {
      {0u, 20u}, {2u, 40u}, {4u, 50u}};

  vecmem::vector<scalar> bin_test_o_edges = {-10.f, 10.f, -20.f,
                                             20.f,  0.f,  100.f};

  vecmem::vector<scalar> same_test_o_edges = {-10.f, 10.f, -20.f,
                                              20.f,  0.f,  100.f};

  vecmem::vector<scalar> diff_test_o_edges = {-10.f, 11.f, -20.f,
                                              20.f,  0.f,  100.f};

  vecmem::vector<dsized_index_range> edge_test_o_ranges = {
      {0u, 20u}, {2u, 40u}, {4u, 50u}};

  // Non-owning multi axis test : reference
  cartesian_3D<is_owning, host_container_types> ref_o_axes(
      std::move(edge_ref_o_ranges), std::move(bin_ref_o_edges));

  // Non-owning multi axis test : test
  cartesian_3D<is_owning, host_container_types> test_o_axes(
      std::move(edge_test_o_ranges), std::move(bin_test_o_edges));

  // Owning multi axis test : reference
  EXPECT_EQ(ref_o_axes, test_o_axes);

  // Non-owning multi axis test : diff endes
  cartesian_3D<is_owning, host_container_types> diff_edges_o_axes(
      std::move(edge_same_o_ranges), std::move(diff_test_o_edges));

  // Confirm that the axes are not equal
  EXPECT_NE(ref_o_axes, diff_edges_o_axes);

  // Non-owning multi axis test : diff ranges
  cartesian_3D<is_owning, host_container_types> diff_edge_ranges_o_axes(
      std::move(edge_diff_o_ranges), std::move(bin_same_o_edges));

  // Confirm that the axes are not equal
  EXPECT_NE(ref_o_axes, diff_edge_ranges_o_axes);
}
