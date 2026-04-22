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

// Project include(s)
#include "detray/definitions/units.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes/line.hpp"
#include "detray/geometry/surface_descriptor.hpp"
#include "detray/material/concepts.hpp"
#include "detray/material/material.hpp"
#include "detray/material/material_rod.hpp"
#include "detray/material/material_slab.hpp"
#include "detray/material/mixture.hpp"
#include "detray/material/predefined_materials.hpp"
#include "detray/navigation/intersection/ray_intersector.hpp"

// Detray test include(s)
#include "detray/test/framework/types.hpp"

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <limits>

using namespace detray;

using test_algebra = test::algebra;
using scalar = test::scalar;
using point2 = test::point2;
using point3 = test::point3;
using transform3 = test::transform3;
using vector3 = test::vector3;

namespace {

constexpr scalar tol{1e-7f};
constexpr auto max_val{std::numeric_limits<scalar>::max()};

}  // anonymous namespace

// This tests the material functionalities
GTEST_TEST(detray_material, materials) {
  // vacuum
  constexpr vacuum<scalar> vac;

  static_assert(concepts::material_params<decltype(vac)>);
  static_assert(concepts::homogeneous_material<decltype(vac)>);
  static_assert(!concepts::surface_material<decltype(vac)>);
  static_assert(concepts::volume_material<decltype(vac)>);
  static_assert(!concepts::material_map<decltype(vac)>);

  EXPECT_EQ(vac.X0(), max_val);
  EXPECT_EQ(vac.L0(), max_val);
  EXPECT_EQ(vac.Ar(), scalar{0});
  EXPECT_EQ(vac.Z(), scalar{0});
  EXPECT_EQ(vac.molar_density(), scalar{0});
  EXPECT_EQ(vac.molar_electron_density(), scalar{0});

  // beryllium
  constexpr beryllium_tml<scalar> beryll;

  static_assert(concepts::material_params<decltype(beryll)>);
  static_assert(concepts::homogeneous_material<decltype(beryll)>);
  static_assert(!concepts::surface_material<decltype(beryll)>);
  static_assert(concepts::volume_material<decltype(beryll)>);
  static_assert(!concepts::material_map<decltype(beryll)>);

  EXPECT_NEAR(beryll.X0(), static_cast<scalar>(352.8 * unit<double>::mm), 1e-4);
  EXPECT_NEAR(beryll.L0(), static_cast<scalar>(407.0 * unit<double>::mm), tol);
  EXPECT_NEAR(beryll.Ar(), 9.012f, tol);
  EXPECT_NEAR(beryll.Z(), 4.0f, tol);

  // @note molar density is obtained by the following equation:
  // molar_density = mass_density [GeV/c²] / molar_mass [GeV/c²] * mol / cm³
  EXPECT_NEAR(beryll.molar_density(),
              static_cast<scalar>(1.848 / beryllium_tml<double>().Ar() *
                                  unit<double>::mol / unit<double>::cm3),
              tol);

  // silicon
  constexpr silicon_tml<scalar> silicon;

  static_assert(concepts::material_params<decltype(silicon)>);
  static_assert(concepts::homogeneous_material<decltype(silicon)>);
  static_assert(!concepts::surface_material<decltype(silicon)>);
  static_assert(concepts::volume_material<decltype(silicon)>);
  static_assert(!concepts::material_map<decltype(silicon)>);

  EXPECT_NEAR(silicon.X0(), static_cast<scalar>(95.7 * unit<double>::mm), 1e-5);
  EXPECT_NEAR(silicon.L0(), static_cast<scalar>(465.2 * unit<double>::mm),
              1e-4);
  EXPECT_NEAR(silicon.Ar(), 28.03f, tol);
  EXPECT_NEAR(silicon.Z(), 14.f, tol);
  EXPECT_NEAR(silicon.molar_density(),
              static_cast<scalar>(2.32 / silicon_tml<double>().Ar() *
                                  unit<double>::mol / unit<double>::cm3),
              tol);
}

GTEST_TEST(detray_material, mixture) {
  // Check if material property doesn't change after mixing with other
  // material of 0 ratio
  constexpr oxygen_gas<scalar> pure_oxygen;
  constexpr mixture<scalar, oxygen_gas<scalar>,
                    aluminium<scalar, std::ratio<0, 1>>>
      oxygen_mix;

  static_assert(concepts::material_params<decltype(oxygen_mix)>);
  static_assert(concepts::homogeneous_material<decltype(oxygen_mix)>);
  static_assert(!concepts::surface_material<decltype(oxygen_mix)>);
  static_assert(concepts::volume_material<decltype(oxygen_mix)>);
  static_assert(!concepts::material_map<decltype(oxygen_mix)>);

  EXPECT_NEAR(pure_oxygen.X0(), oxygen_mix.X0(), 0.02f);
  EXPECT_NEAR(pure_oxygen.L0(), oxygen_mix.L0(), tol);
  EXPECT_NEAR(pure_oxygen.Ar(), oxygen_mix.Ar(), tol);
  EXPECT_NEAR(pure_oxygen.Z(), oxygen_mix.Z(), tol);
  EXPECT_NEAR(pure_oxygen.mass_density(), oxygen_mix.mass_density(), tol);
  EXPECT_NEAR(pure_oxygen.molar_density(), oxygen_mix.molar_density(), tol);

  // Air mixture check
  constexpr mixture<scalar, carbon_gas<scalar, std::ratio<0, 100>>,
                    nitrogen_gas<scalar, std::ratio<76, 100>>,
                    oxygen_gas<scalar, std::ratio<23, 100>>,
                    argon_gas<scalar, std::ratio<1, 100>>>
      air_mix;
  constexpr air<scalar> pure_air;

  static_assert(concepts::material_params<decltype(pure_air)>);
  static_assert(concepts::homogeneous_material<decltype(pure_air)>);
  static_assert(!concepts::surface_material<decltype(pure_air)>);
  static_assert(concepts::volume_material<decltype(pure_air)>);
  static_assert(!concepts::material_map<decltype(pure_air)>);

  EXPECT_TRUE(std::abs(air_mix.X0() - pure_air.X0()) / pure_air.X0() < 0.01f);
  EXPECT_TRUE(std::abs(air_mix.L0() - pure_air.L0()) / pure_air.L0() < 0.01f);
  EXPECT_TRUE(std::abs(air_mix.Ar() - pure_air.Ar()) / pure_air.Ar() < 0.01f);
  EXPECT_TRUE(std::abs(air_mix.Z() - pure_air.Z()) / pure_air.Z() < 0.01f);
  EXPECT_TRUE(std::abs(air_mix.mass_density() - pure_air.mass_density()) /
                  pure_air.mass_density() <
              0.01f);

  // Vector check
  material_slab<scalar> slab1(air_mix, 5.5f);
  material_slab<scalar> slab2(pure_air, 2.3f);
  material_slab<scalar> slab3(oxygen_gas<scalar>(), 2.f);

  std::vector<material_slab<scalar>> slab_vec;

  slab_vec.push_back(slab1);
  slab_vec.push_back(slab2);
  slab_vec.push_back(slab3);

  EXPECT_NEAR(slab_vec[0].thickness_in_X0(),
              slab1.thickness() / slab1.get_material().X0(), tol);
  EXPECT_NEAR(slab_vec[1].thickness_in_X0(),
              slab2.thickness() / slab2.get_material().X0(), tol);
  EXPECT_NEAR(slab_vec[2].thickness_in_X0(),
              slab3.thickness() / slab3.get_material().X0(), tol);
}

GTEST_TEST(detray_material, mixture2) {
  constexpr mixture<scalar, silicon_with_ded<scalar, std::ratio<1, 3>>,
                    cesium_iodide_with_ded<scalar, std::ratio<2, 3>>>
      mix;

  static_assert(concepts::material_params<decltype(mix)>);
  static_assert(concepts::homogeneous_material<decltype(mix)>);
  static_assert(!concepts::surface_material<decltype(mix)>);
  static_assert(concepts::volume_material<decltype(mix)>);
  static_assert(!concepts::material_map<decltype(mix)>);

  // For the moment, the mixture is not supposed to have density effect data
  // in any case
  EXPECT_EQ(mix.has_density_effect_data(), false);
  EXPECT_EQ(mix.Ar(), silicon<scalar>().Ar() * 1.f / 3.f +
                          cesium_iodide<scalar>().Ar() * 2.f / 3.f);
  EXPECT_EQ(mix.Z(), silicon<scalar>().Z() * 1.f / 3.f +
                         cesium_iodide<scalar>().Z() * 2.f / 3.f);
}

// This tests the material slab functionalities
GTEST_TEST(detray_material, material_slab) {
  constexpr material_slab<scalar> slab(oxygen_gas<scalar>(),
                                       2.f * unit<scalar>::mm);

  static_assert(concepts::material_slab<decltype(slab)>);
  static_assert(concepts::homogeneous_material<decltype(slab)>);
  static_assert(concepts::surface_material<decltype(slab)>);
  static_assert(!concepts::volume_material<decltype(slab)>);
  static_assert(!concepts::material_map<decltype(slab)>);

  const scalar cos_inc_ang{0.3f};

  EXPECT_NEAR(slab.path_segment(cos_inc_ang), 2.f * unit<scalar>::mm / 0.3f,
              tol);
  EXPECT_NEAR(slab.path_segment_in_X0(cos_inc_ang),
              slab.path_segment(cos_inc_ang) / slab.get_material().X0(), tol);
  EXPECT_NEAR(slab.path_segment_in_L0(cos_inc_ang),
              slab.path_segment(cos_inc_ang) / slab.get_material().L0(), tol);
}

// This tests the material rod functionalities
GTEST_TEST(detray_material, material_rod) {
  // Rod with 1 mm radius
  constexpr material_rod<scalar> rod(oxygen_gas<scalar>(),
                                     1.f * unit<scalar>::mm);

  static_assert(concepts::material_rod<decltype(rod)>);
  static_assert(concepts::homogeneous_material<decltype(rod)>);
  static_assert(concepts::surface_material<decltype(rod)>);
  static_assert(!concepts::volume_material<decltype(rod)>);
  static_assert(!concepts::material_map<decltype(rod)>);

  // tf3 with Identity rotation and no translation
  const vector3 x{1.f, 0.f, 0.f};
  const vector3 z{0.f, 0.f, 1.f};
  const vector3 t{0.f, 0.f, 0.f};
  const transform3 tf{t, vector::normalize(z), vector::normalize(x)};

  // Create a track
  const point3 pos{-1.f / 6.f, -10.f, 0.f};
  const vector3 dir{0.f, 1.f, 3.f};
  const free_track_parameters<test_algebra> trk(pos, 0.f, dir, -1.f);

  // Infinite wire with 1 mm radial cell size
  const mask<line_circular, test_algebra> ln{
      0u, 1.f * unit<scalar>::mm, std::numeric_limits<scalar>::infinity()};

  auto is = ray_intersector<line_circular, test_algebra, true>{}(
      detail::ray<test_algebra>(trk), surface_descriptor<>{}, ln, tf);

  const scalar cos_inc_ang{std::abs(vector::dot(
      line2D<test_algebra>::normal(tf, is.local()), vector::normalize(dir)))};
  const scalar approach{is.local()[0]};

  EXPECT_NEAR(rod.path_segment(cos_inc_ang, approach),
              2.f * std::sqrt(10.f - 10.f / 36.f), 1e-5f);
  EXPECT_NEAR(rod.path_segment_in_X0(cos_inc_ang, approach),
              rod.path_segment(cos_inc_ang, approach) / rod.get_material().X0(),
              tol);
  EXPECT_NEAR(rod.path_segment_in_L0(cos_inc_ang, approach),
              rod.path_segment(cos_inc_ang, approach) / rod.get_material().L0(),
              tol);
}
