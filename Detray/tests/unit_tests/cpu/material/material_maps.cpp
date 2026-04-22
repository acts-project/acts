// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Detray include(s)
#include "detray/definitions/indexing.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes.hpp"
#include "detray/material/concepts.hpp"
#include "detray/material/material_map.hpp"

// Detray test include(s)
#include "detray/test/framework/types.hpp"

// GTest include(s)
#include <gtest/gtest.h>

using namespace detray;
using namespace detray::axis;

using test_algebra = test::algebra;
using scalar = test::scalar;

using material_t =
    typename material_grid_factory<test_algebra>::bin_type::entry_type;

namespace {

material_grid_factory<test_algebra> mat_map_factory{};

}  // anonymous namespace

/// Unittest: Test the construction of an annulus shaped material map
GTEST_TEST(detray_material, annulus_map) {
  constexpr scalar minR{7.2f * unit<scalar>::mm};
  constexpr scalar maxR{12.0f * unit<scalar>::mm};
  constexpr scalar minPhi{0.74195f};
  constexpr scalar maxPhi{1.33970f};

  mask<annulus2D, test_algebra> ann2{0u,     minR, maxR, minPhi,
                                     maxPhi, 0.f,  -2.f, 2.f};

  auto annulus_map = mat_map_factory.new_grid(ann2, {10u, 20u});

  static_assert(concepts::material_map<decltype(annulus_map)>);
  static_assert(!concepts::homogeneous_material<decltype(annulus_map)>);
  static_assert(concepts::surface_material<decltype(annulus_map)>);
  static_assert(!concepts::volume_material<decltype(annulus_map)>);

  EXPECT_EQ(annulus_map.dim, 2u);
  EXPECT_EQ(annulus_map.nbins(), 200u);

  auto r_axis = annulus_map.get_axis<label::e_r>();
  EXPECT_EQ(r_axis.nbins(), 10u);
  EXPECT_EQ(r_axis.min(), minR);
  EXPECT_EQ(r_axis.max(), maxR);

  auto phi_axis = annulus_map.get_axis<label::e_phi>();
  EXPECT_EQ(phi_axis.nbins(), 20u);
  EXPECT_EQ(phi_axis.min(), -minPhi);
  EXPECT_EQ(phi_axis.max(), maxPhi);

  // Add some material entries
  scalar thickness = 2.f * unit<scalar>::mm;
  for (dindex gbin = 0; gbin < annulus_map.nbins(); ++gbin) {
    annulus_map.template populate<replace<>>(
        gbin, material_t(silicon_tml<scalar>{}, thickness));
    thickness += 1.f * unit<scalar>::mm;
  }

  EXPECT_EQ(annulus_map.at(0, 0),
            material_t(silicon_tml<scalar>{}, 2.f * unit<scalar>::mm));
  EXPECT_EQ(annulus_map.at(1, 0),
            material_t(silicon_tml<scalar>{}, 3.f * unit<scalar>::mm));
  EXPECT_EQ(annulus_map.at(22, 0),
            material_t(silicon_tml<scalar>{}, 24.f * unit<scalar>::mm));
  EXPECT_EQ(annulus_map.at(199, 0),
            material_t(silicon_tml<scalar>{}, 201.f * unit<scalar>::mm));
  EXPECT_FALSE(annulus_map.at(199, 0) ==
               material_t(gold<scalar>{}, 201.f * unit<scalar>::mm));
}

/// Unittest: Test the construction of a cylindrical material map
GTEST_TEST(detray_material, cylinder_map) {
  constexpr scalar r{3.f * unit<scalar>::mm};
  constexpr scalar hz{4.f * unit<scalar>::mm};

  mask<cylinder2D, test_algebra> cyl{0u, r, -hz, hz};

  auto cylinder_map = mat_map_factory.new_grid(cyl, {10u, 20u});

  static_assert(concepts::material_map<decltype(cylinder_map)>);
  static_assert(!concepts::homogeneous_material<decltype(cylinder_map)>);
  static_assert(concepts::surface_material<decltype(cylinder_map)>);
  static_assert(!concepts::volume_material<decltype(cylinder_map)>);

  EXPECT_EQ(cylinder_map.dim, 2u);
  EXPECT_EQ(cylinder_map.nbins(), 200u);

  auto rphi_axis = cylinder_map.get_axis<label::e_rphi>();
  EXPECT_EQ(rphi_axis.nbins(), 10u);
  EXPECT_EQ(rphi_axis.min(), -constant<scalar>::pi);
  EXPECT_EQ(rphi_axis.max(), constant<scalar>::pi);

  auto z_axis = cylinder_map.get_axis<label::e_cyl_z>();
  EXPECT_EQ(z_axis.nbins(), 20u);
  EXPECT_EQ(z_axis.min(), -hz);
  EXPECT_EQ(z_axis.max(), hz);

  scalar thickness = 2.f * unit<scalar>::mm;
  for (dindex gbin = 0; gbin < cylinder_map.nbins(); ++gbin) {
    cylinder_map.template populate<replace<>>(
        gbin, material_t(vacuum<scalar>{}, thickness));
    thickness += 1.f * unit<scalar>::mm;
  }

  EXPECT_EQ(cylinder_map.at(0, 0),
            material_t(vacuum<scalar>{}, 2.f * unit<scalar>::mm));
  EXPECT_EQ(cylinder_map.at(1, 0),
            material_t(vacuum<scalar>{}, 3.f * unit<scalar>::mm));
  EXPECT_EQ(cylinder_map.at(22, 0),
            material_t(vacuum<scalar>{}, 24.f * unit<scalar>::mm));
  EXPECT_EQ(cylinder_map.at(199, 0),
            material_t(vacuum<scalar>{}, 201.f * unit<scalar>::mm));
  EXPECT_FALSE(cylinder_map.at(199, 0) ==
               material_t(gold<scalar>{}, 201.f * unit<scalar>::mm));
}

/// Unittest: Test the construction of a rectangular material map
GTEST_TEST(detray_material, rectangle_map) {
  constexpr scalar hx{1.f * unit<scalar>::mm};
  constexpr scalar hy{9.3f * unit<scalar>::mm};

  mask<rectangle2D, test_algebra> r2{0u, hx, hy};

  auto rectangle_map = mat_map_factory.new_grid(r2, {10u, 20u});

  static_assert(concepts::material_map<decltype(rectangle_map)>);
  static_assert(!concepts::homogeneous_material<decltype(rectangle_map)>);
  static_assert(concepts::surface_material<decltype(rectangle_map)>);
  static_assert(!concepts::volume_material<decltype(rectangle_map)>);

  EXPECT_EQ(rectangle_map.dim, 2u);
  EXPECT_EQ(rectangle_map.nbins(), 200u);

  auto x_axis = rectangle_map.get_axis<label::e_x>();
  EXPECT_EQ(x_axis.nbins(), 10u);
  EXPECT_EQ(x_axis.min(), -hx);
  EXPECT_EQ(x_axis.max(), hx);

  auto y_axis = rectangle_map.get_axis<label::e_y>();
  EXPECT_EQ(y_axis.nbins(), 20u);
  EXPECT_EQ(y_axis.min(), -hy);
  EXPECT_EQ(y_axis.max(), hy);

  scalar thickness = 2.f * unit<scalar>::mm;
  for (dindex gbin = 0; gbin < rectangle_map.nbins(); ++gbin) {
    rectangle_map.template populate<replace<>>(
        gbin, material_t(oxygen_gas<scalar>{}, thickness));
    thickness += 1.f * unit<scalar>::mm;
  }

  EXPECT_EQ(rectangle_map.at(0, 0),
            material_t(oxygen_gas<scalar>{}, 2.f * unit<scalar>::mm));
  EXPECT_EQ(rectangle_map.at(1, 0),
            material_t(oxygen_gas<scalar>{}, 3.f * unit<scalar>::mm));
  EXPECT_EQ(rectangle_map.at(22, 0),
            material_t(oxygen_gas<scalar>{}, 24.f * unit<scalar>::mm));
  EXPECT_EQ(rectangle_map.at(199, 0),
            material_t(oxygen_gas<scalar>{}, 201.f * unit<scalar>::mm));
  EXPECT_FALSE(rectangle_map.at(199, 0) ==
               material_t(gold<scalar>{}, 201.f * unit<scalar>::mm));
}

/// Unittest: Test the construction of a ring shaped material map
GTEST_TEST(detray_material, disc_map) {
  constexpr scalar inner_r{0.f * unit<scalar>::mm};
  constexpr scalar outer_r{3.5f * unit<scalar>::mm};

  mask<ring2D, test_algebra> r2{0u, inner_r, outer_r};

  auto disc_map = mat_map_factory.new_grid(r2, {10u, 20u});

  static_assert(concepts::material_map<decltype(disc_map)>);
  static_assert(!concepts::homogeneous_material<decltype(disc_map)>);
  static_assert(concepts::surface_material<decltype(disc_map)>);
  static_assert(!concepts::volume_material<decltype(disc_map)>);

  EXPECT_EQ(disc_map.dim, 2u);
  EXPECT_EQ(disc_map.nbins(), 200u);

  auto r_axis = disc_map.get_axis<label::e_r>();
  EXPECT_EQ(r_axis.nbins(), 10u);
  EXPECT_EQ(r_axis.min(), inner_r);
  EXPECT_EQ(r_axis.max(), outer_r);

  auto phi_axis = disc_map.get_axis<label::e_phi>();
  EXPECT_EQ(phi_axis.nbins(), 20u);
  EXPECT_EQ(phi_axis.min(), -constant<scalar>::pi);
  EXPECT_EQ(phi_axis.max(), constant<scalar>::pi);

  scalar thickness = 2.f * unit<scalar>::mm;
  for (dindex gbin = 0; gbin < disc_map.nbins(); ++gbin) {
    disc_map.template populate<replace<>>(
        gbin, material_t(aluminium<scalar>{}, thickness));
    thickness += 1.f * unit<scalar>::mm;
  }

  EXPECT_EQ(disc_map.at(0, 0),
            material_t(aluminium<scalar>{}, 2.f * unit<scalar>::mm));
  EXPECT_EQ(disc_map.at(1, 0),
            material_t(aluminium<scalar>{}, 3.f * unit<scalar>::mm));
  EXPECT_EQ(disc_map.at(22, 0),
            material_t(aluminium<scalar>{}, 24.f * unit<scalar>::mm));
  EXPECT_EQ(disc_map.at(199, 0),
            material_t(aluminium<scalar>{}, 201.f * unit<scalar>::mm));
  EXPECT_FALSE(disc_map.at(199, 0) ==
               material_t(gold<scalar>{}, 201.f * unit<scalar>::mm));
}

/// Unittest: Test the construction of a trapezoid material map
GTEST_TEST(detray_material, trapezoid_map) {
  constexpr scalar hx_miny{1.f * unit<scalar>::mm};
  constexpr scalar hx_maxy{3.f * unit<scalar>::mm};
  constexpr scalar hy{2.f * unit<scalar>::mm};
  constexpr scalar divisor{1.f / (2.f * hy)};

  mask<trapezoid2D, test_algebra> t2{0u, hx_miny, hx_maxy, hy, divisor};

  auto trapezoid_map = mat_map_factory.new_grid(t2, {10u, 20u});

  static_assert(concepts::material_map<decltype(trapezoid_map)>);
  static_assert(!concepts::homogeneous_material<decltype(trapezoid_map)>);
  static_assert(concepts::surface_material<decltype(trapezoid_map)>);
  static_assert(!concepts::volume_material<decltype(trapezoid_map)>);

  EXPECT_EQ(trapezoid_map.dim, 2u);
  EXPECT_EQ(trapezoid_map.nbins(), 200u);

  auto x_axis = trapezoid_map.get_axis<label::e_x>();
  EXPECT_EQ(x_axis.nbins(), 10u);
  EXPECT_EQ(x_axis.min(), -hx_maxy);
  EXPECT_EQ(x_axis.max(), hx_maxy);

  auto y_axis = trapezoid_map.get_axis<label::e_y>();
  EXPECT_EQ(y_axis.nbins(), 20u);
  EXPECT_EQ(y_axis.min(), -hy);
  EXPECT_EQ(y_axis.max(), hy);

  scalar thickness = 2.f * unit<scalar>::mm;
  for (dindex gbin = 0; gbin < trapezoid_map.nbins(); ++gbin) {
    trapezoid_map.template populate<replace<>>(
        gbin, material_t(gold<scalar>{}, thickness));
    thickness += 1.f * unit<scalar>::mm;
  }

  EXPECT_EQ(trapezoid_map.at(0, 0),
            material_t(gold<scalar>{}, 2.f * unit<scalar>::mm));
  EXPECT_EQ(trapezoid_map.at(1, 0),
            material_t(gold<scalar>{}, 3.f * unit<scalar>::mm));
  EXPECT_EQ(trapezoid_map.at(22, 0),
            material_t(gold<scalar>{}, 24.f * unit<scalar>::mm));
  EXPECT_EQ(trapezoid_map.at(199, 0),
            material_t(gold<scalar>{}, 201.f * unit<scalar>::mm));
  EXPECT_FALSE(trapezoid_map.at(199, 0) ==
               material_t(aluminium<scalar>{}, 201.f * unit<scalar>::mm));
}

/// Unittest: Test the material grid comparisons
GTEST_TEST(detray_material, material_grid_comparison) {
  /** Allows to create a regular grid to check the equality operator
   * grids can differ in:
   * - type (will never be compared)
   * - bins/axes
   * - entries (i.e. material data)
   */
  auto createGrid = [](const scalar hx, const scalar hy, unsigned int bx,
                       unsigned int by, bool distort_entries) {
    mask<rectangle2D, test_algebra> r2{0u, hx, hy};
    auto material_grid = mat_map_factory.new_grid(r2, {bx, by});

    // Fill the material grid with some data
    scalar thickness = 2.f * unit<scalar>::mm;
    for (dindex gbin = 0; gbin < material_grid.nbins(); ++gbin) {
      material_grid.template populate<replace<>>(
          gbin, material_t(oxygen_gas<scalar>{}, thickness));
      thickness += 1.f * unit<scalar>::mm;
      if (distort_entries) {
        thickness += 1.f * unit<scalar>::mm;
      }
    }
    // Return it
    return material_grid;
  };

  // Two equal grids
  auto grid_ref = createGrid(10.f, 20.f, 10u, 20u, false);
  auto grid_eq = createGrid(10.f, 20.f, 10u, 20u, false);
  EXPECT_EQ(grid_ref, grid_eq);

  // One grid with different size
  auto grid_neq_size = createGrid(11.f, 21.f, 10u, 20u, false);
  EXPECT_NE(grid_ref, grid_neq_size);

  // One grid with different binning
  auto grid_neq_bins = createGrid(10.f, 20.f, 11u, 21u, false);
  EXPECT_NE(grid_ref, grid_neq_bins);

  // One grid with different entries
  auto grid_neq_entries = createGrid(10.f, 20.f, 10u, 20u, true);
  EXPECT_NE(grid_ref, grid_neq_entries);
}
