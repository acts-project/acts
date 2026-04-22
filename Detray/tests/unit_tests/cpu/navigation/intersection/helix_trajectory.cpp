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
#include "detray/tracks/helix.hpp"
#include "detray/tracks/tracks.hpp"

// Detray test include(s)
#include "detray/test/framework/types.hpp"

// GTest include(s)
#include <gtest/gtest.h>

using namespace detray;

using scalar = test::scalar;

constexpr const scalar tol{1e-5f};

// This tests the base functionality of the Helix Gun
GTEST_TEST(detray_intersection, helix_trajectory) {
  using test_algebra = test::algebra;
  using vector3 = test::vector3;
  using point3 = test::point3;

  const point3 pos{0.f, 0.f, 0.f};
  const scalar time{0.f};
  const vector3 mom{static_cast<scalar>(1.f), static_cast<scalar>(0.f),
                    1.f * unit<scalar>::GeV};
  const scalar q{static_cast<scalar>(-1.) * unit<scalar>::e};

  // vertex
  free_track_parameters<test_algebra> vertex(pos, time, mom, q);

  // magnetic field
  const vector3 B{static_cast<scalar>(0.f), static_cast<scalar>(0.f),
                  1.f * unit<scalar>::T};

  const scalar p_mag{vector::norm(mom)};
  const scalar B_mag{vector::norm(B)};
  const scalar pz_along{vector::dot(mom, vector::normalize(B))};
  const scalar pt{std::sqrt(p_mag * p_mag - pz_along * pz_along)};

  // helix trajectory
  detail::helix helix_traj(vertex, B);
  EXPECT_NEAR(helix_traj.time(), 0.f, tol);
  EXPECT_NEAR(helix_traj.qop(), -constant<scalar>::inv_sqrt2, tol);

  // radius of helix
  scalar R{helix_traj.radius()};
  EXPECT_NEAR(R, pt / B_mag, tol);

  // (1 T of b_z, 1 GeV/c of p_T) ==> 3.33564095 m of radius
  EXPECT_NEAR(R, 3.33564095f * unit<scalar>::m, 10.f * tol);

  // Path length for one loop
  scalar S = 2.f * p_mag / B_mag * constant<scalar>::pi;

  // After half turn
  point3 half_loop_pos = helix_traj(S / 2.f);
  EXPECT_NEAR(half_loop_pos[0], 0.f, R * tol);
  EXPECT_NEAR(half_loop_pos[1], 2.f * R, R * tol);
  EXPECT_NEAR(half_loop_pos[2], pz_along / B_mag * constant<scalar>::pi,
              R * tol);

  point3 half_loop_dir = helix_traj.dir(S / 2.f);
  EXPECT_NEAR(half_loop_dir[0], -vertex.dir()[0], R * tol);
  EXPECT_NEAR(half_loop_dir[1], -vertex.dir()[1], R * tol);
  EXPECT_NEAR(half_loop_dir[2], vertex.dir()[2], R * tol);

  // After half turn in the opposite direction
  half_loop_pos = helix_traj(-S / 2.f);
  EXPECT_NEAR(half_loop_pos[0], 0.f, R * tol);
  EXPECT_NEAR(half_loop_pos[1], 2.f * R, R * tol);
  EXPECT_NEAR(half_loop_pos[2], -pz_along / B_mag * constant<scalar>::pi,
              R * tol);

  half_loop_dir = helix_traj.dir(-S / 2.f);
  EXPECT_NEAR(half_loop_dir[0], -vertex.dir()[0], R * tol);
  EXPECT_NEAR(half_loop_dir[1], -vertex.dir()[1], R * tol);
  EXPECT_NEAR(half_loop_dir[2], vertex.dir()[2], R * tol);

  // After one full turn
  point3 one_loop_pos = helix_traj(S);
  EXPECT_NEAR(one_loop_pos[0], 0.f, R * tol);
  EXPECT_NEAR(one_loop_pos[1], 0.f, R * tol);
  EXPECT_NEAR(one_loop_pos[2], 2.f * pz_along / B_mag * constant<scalar>::pi,
              R * tol);

  point3 one_loop_dir = helix_traj.dir(S);
  EXPECT_NEAR(one_loop_dir[0], vertex.dir()[0], R * tol);
  EXPECT_NEAR(one_loop_dir[1], vertex.dir()[1], R * tol);
  EXPECT_NEAR(one_loop_dir[2], vertex.dir()[2], R * tol);

  // After one full turn in the opposite direction
  one_loop_pos = helix_traj(-S);
  EXPECT_NEAR(one_loop_pos[0], 0.f, R * tol);
  EXPECT_NEAR(one_loop_pos[1], 0.f, R * tol);
  EXPECT_NEAR(one_loop_pos[2], -2.f * pz_along / B_mag * constant<scalar>::pi,
              R * tol);

  one_loop_dir = helix_traj.dir(-S);
  EXPECT_NEAR(one_loop_dir[0], vertex.dir()[0], R * tol);
  EXPECT_NEAR(one_loop_dir[1], vertex.dir()[1], R * tol);
  EXPECT_NEAR(one_loop_dir[2], vertex.dir()[2], R * tol);

  /*********************************
   * Same test with oppsite charge
   *********************************/

  free_track_parameters<test_algebra> vertex2(pos, time, mom, -q);

  // helix trajectory
  detail::helix helix_traj2(vertex2, B);

  EXPECT_NEAR(R, helix_traj2.radius(), tol);

  // After half turn
  half_loop_pos = helix_traj2(S / 2.f);
  EXPECT_NEAR(half_loop_pos[0], 0.f, R * tol);
  EXPECT_NEAR(half_loop_pos[1], -2.f * R, R * tol);
  EXPECT_NEAR(half_loop_pos[2], pz_along / B_mag * constant<scalar>::pi,
              R * tol);

  half_loop_dir = helix_traj2.dir(S / 2.f);
  EXPECT_NEAR(half_loop_dir[0], -vertex2.dir()[0], R * tol);
  EXPECT_NEAR(half_loop_dir[1], -vertex2.dir()[1], R * tol);
  EXPECT_NEAR(half_loop_dir[2], vertex2.dir()[2], R * tol);

  // After half turn in the opposite direction
  half_loop_pos = helix_traj2(-S / 2.f);
  EXPECT_NEAR(half_loop_pos[0], 0.f, R * tol);
  EXPECT_NEAR(half_loop_pos[1], -2.f * R, R * tol);
  EXPECT_NEAR(half_loop_pos[2], -pz_along / B_mag * constant<scalar>::pi,
              R * tol);

  half_loop_dir = helix_traj.dir(-S / 2.f);
  EXPECT_NEAR(half_loop_dir[0], -vertex.dir()[0], R * tol);
  EXPECT_NEAR(half_loop_dir[1], -vertex.dir()[1], R * tol);
  EXPECT_NEAR(half_loop_dir[2], vertex.dir()[2], R * tol);

  // After one full turn
  one_loop_pos = helix_traj2(S);
  EXPECT_NEAR(one_loop_pos[0], 0.f, R * tol);
  EXPECT_NEAR(one_loop_pos[1], 0.f, R * tol);
  EXPECT_NEAR(one_loop_pos[2], 2.f * pz_along / B_mag * constant<scalar>::pi,
              R * tol);

  one_loop_dir = helix_traj2.dir(S);
  EXPECT_NEAR(one_loop_dir[0], vertex2.dir()[0], R * tol);
  EXPECT_NEAR(one_loop_dir[1], vertex2.dir()[1], R * tol);
  EXPECT_NEAR(one_loop_dir[2], vertex2.dir()[2], R * tol);

  // After one full turn in the opposite direction
  one_loop_pos = helix_traj2(-S);
  EXPECT_NEAR(one_loop_pos[0], 0.f, R * tol);
  EXPECT_NEAR(one_loop_pos[1], 0.f, R * tol);
  EXPECT_NEAR(one_loop_pos[2], -2.f * pz_along / B_mag * constant<scalar>::pi,
              R * tol);

  one_loop_dir = helix_traj2.dir(-S);
  EXPECT_NEAR(one_loop_dir[0], vertex2.dir()[0], R * tol);
  EXPECT_NEAR(one_loop_dir[1], vertex2.dir()[1], R * tol);
  EXPECT_NEAR(one_loop_dir[2], vertex2.dir()[2], R * tol);
}

GTEST_TEST(detray_intersection, helix_trajectory_small_pT) {
  using test_algebra = test::algebra;
  using vector3 = test::vector3;
  using point3 = test::point3;

  const point3 pos{0.f, 0.f, 0.f};
  const scalar time{0.f};
  const vector3 mom{static_cast<scalar>(0.f), tol, 1.f * unit<scalar>::GeV};
  const scalar q{static_cast<scalar>(-1.) * unit<scalar>::e};

  // vertex
  free_track_parameters<test_algebra> vertex(pos, time, mom, q);

  // magnetic field
  const vector3 B{static_cast<scalar>(0.f), static_cast<scalar>(0.f),
                  1.f * unit<scalar>::T};

  // helix trajectory
  detail::helix helix_traj(vertex, B);
  EXPECT_NEAR(helix_traj.time(), 0.f, tol);
  EXPECT_NEAR(helix_traj.qop(), -1.f, tol);

  // After 10 mm
  const scalar path_length{10.f * unit<scalar>::mm};
  const point3 helix_pos = helix_traj(path_length);
  const point3 true_pos = pos + path_length * vector::normalize(mom);

  EXPECT_NEAR(true_pos[0], helix_pos[0], tol);
  EXPECT_NEAR(true_pos[1], helix_pos[1], tol);
  EXPECT_NEAR(true_pos[2], helix_pos[2], tol);
}

GTEST_TEST(detray_intersection, helix_direction_stability) {
  using test_algebra = test::algebra;
  using vector3 = test::vector3;
  using point3 = test::point3;

  // magnetic field
  const vector3 B{static_cast<scalar>(0.f), static_cast<scalar>(0.f),
                  1.f * unit<scalar>::T};

  const point3 pos{0.f, 0.f, 0.f};
  const scalar time{0.f};
  const vector3 mom{1.f * unit<scalar>::GeV, 1.f * unit<scalar>::GeV,
                    1.f * unit<scalar>::GeV};
  const scalar q{static_cast<scalar>(-1.) * unit<scalar>::e};

  // vertex
  free_track_parameters<test_algebra> vertex(pos, time, mom, q);

  // helix trajectory
  detail::helix hlx(vertex, B);

  for (int i = 0; i < 100; i++) {
    const auto d = hlx.dir(scalar(i) * 10.f);
    ASSERT_FLOAT_EQ(static_cast<float>(vector::theta(d)),
                    static_cast<float>(vector::theta(mom)));
  }
}
