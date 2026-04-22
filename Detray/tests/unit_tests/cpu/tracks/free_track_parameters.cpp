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
#include "detray/tracks/free_track_parameters.hpp"

// Detray test include(s)
#include "detray/test/framework/types.hpp"

// Google Test include(s)
#include <gtest/gtest.h>

using namespace detray;

using test_algebra = test::algebra;
using scalar = test::scalar;
using vector3 = test::vector3;
using point3 = test::point3;

using free_vector_t = typename free_track_parameters<test_algebra>::vector_type;

constexpr scalar tol{1e-5f};

GTEST_TEST(detray_tracks, free_track_parameters) {
  constexpr scalar charge = -1.f;

  point3 pos = {4.f, 10.f, 2.f};
  scalar time = 0.1f;
  vector3 mom = {10.f, 20.f, 30.f};
  vector3 norm_mom = vector::normalize(mom);

  auto free_vec = matrix::zero<free_vector_t>();
  getter::element(free_vec, e_free_pos0, 0u) = pos[0];
  getter::element(free_vec, e_free_pos1, 0u) = pos[1];
  getter::element(free_vec, e_free_pos2, 0u) = pos[2];
  getter::element(free_vec, e_free_time, 0u) = time;
  getter::element(free_vec, e_free_dir0, 0u) = norm_mom[0];
  getter::element(free_vec, e_free_dir1, 0u) = norm_mom[1];
  getter::element(free_vec, e_free_dir2, 0u) = norm_mom[2];
  getter::element(free_vec, e_free_qoverp, 0u) = charge / vector::norm(mom);

  // first constructor
  free_track_parameters<test_algebra> free_param1(free_vec);
  EXPECT_NEAR(free_param1.pos()[0], pos[0], tol);
  EXPECT_NEAR(free_param1.pos()[1], pos[1], tol);
  EXPECT_NEAR(free_param1.pos()[2], pos[2], tol);
  EXPECT_NEAR(free_param1.dir()[0], norm_mom[0], tol);
  EXPECT_NEAR(free_param1.dir()[1], norm_mom[1], tol);
  EXPECT_NEAR(free_param1.dir()[2], norm_mom[2], tol);
  EXPECT_NEAR(free_param1.time(), time, tol);
  EXPECT_NEAR(free_param1.qop(), charge / vector::norm(mom), tol);

  EXPECT_NEAR(vector::norm(free_param1.mom(charge)), vector::norm(mom), tol);
  EXPECT_NEAR(free_param1.pT(charge), vector::perp(mom), tol);
  EXPECT_NEAR(free_param1.qopT(), charge / free_param1.pT(charge), tol);
  EXPECT_NEAR(free_param1.pz(charge), mom[2], tol);
  EXPECT_NEAR(free_param1.qopz(), charge / free_param1.pz(charge), tol);
  EXPECT_NEAR(free_param1.mom(charge)[0],
              free_param1.p(charge) * free_param1.dir()[0], tol);
  EXPECT_NEAR(free_param1.mom(charge)[1],
              free_param1.p(charge) * free_param1.dir()[1], tol);
  EXPECT_NEAR(free_param1.mom(charge)[2],
              free_param1.p(charge) * free_param1.dir()[2], tol);

  // second constructor
  free_track_parameters<test_algebra> free_param2(pos, time, mom, charge);
  EXPECT_NEAR(free_param2.pos()[0], pos[0], tol);
  EXPECT_NEAR(free_param2.pos()[1], pos[1], tol);
  EXPECT_NEAR(free_param2.pos()[2], pos[2], tol);
  EXPECT_NEAR(free_param2.dir()[0], norm_mom[0], tol);
  EXPECT_NEAR(free_param2.dir()[1], norm_mom[1], tol);
  EXPECT_NEAR(free_param2.dir()[2], norm_mom[2], tol);
  EXPECT_NEAR(free_param2.time(), time, tol);
  EXPECT_NEAR(free_param2.qop(), charge / vector::norm(mom), tol);

  EXPECT_NEAR(vector::norm(free_param2.mom(charge)), vector::norm(mom), tol);
  EXPECT_NEAR(free_param2.pT(charge), vector::perp(mom), tol);
  EXPECT_NEAR(free_param2.mom(charge)[0],
              free_param2.p(charge) * free_param2.dir()[0], tol);
  EXPECT_NEAR(free_param2.mom(charge)[1],
              free_param2.p(charge) * free_param2.dir()[1], tol);
  EXPECT_NEAR(free_param2.mom(charge)[2],
              free_param2.p(charge) * free_param2.dir()[2], tol);

  EXPECT_TRUE(free_param2 == free_param1);

  // Test the setters and subscript operator
  pos = {1.f, 2.f, 3.f};
  time = 0.5f;
  const vector3 dir{vector::normalize(vector3{40.f, 50.f, 60.f})};
  const scalar qop{0.634f};

  free_param2.set_pos(pos);
  free_param2.set_dir(dir);
  free_param2.set_time(time);
  free_param2.set_qop(qop);

  EXPECT_NEAR(free_param2[e_free_pos0], pos[0], tol);
  EXPECT_NEAR(free_param2[e_free_pos1], pos[1], tol);
  EXPECT_NEAR(free_param2[e_free_pos2], pos[2], tol);
  EXPECT_NEAR(free_param2[e_free_time], time, tol);
  EXPECT_NEAR(free_param2[e_free_dir0], dir[0], tol);
  EXPECT_NEAR(free_param2[e_free_dir1], dir[1], tol);
  EXPECT_NEAR(free_param2[e_free_dir2], dir[2], tol);
  EXPECT_NEAR(free_param2[e_free_qoverp], qop, tol);
}
