// This file is part of the ACTS project.
//
// Copyright (C) 2017-2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE BoundaryCheck Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include <boost/test/data/test_case.hpp>
// leave blank line

#include <boost/test/output_test_stream.hpp>
// leave blank line

#include "ACTS/Surfaces/BoundaryCheck.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Units.hpp"

namespace Acts {
namespace Test {
  BOOST_AUTO_TEST_SUITE(Surfaces)
  // See: https://en.wikipedia.org/wiki/Bounding_volume
  //
  // Aligned box w/ simple check
  BOOST_AUTO_TEST_CASE(BoundaryCheckBoxSimple)
  {
    BoundaryCheck check(true);
    BOOST_TEST(check.isInside({0, 0}, -1, 1, -1, 1) == true);
    BOOST_TEST(check.isInside({2, 2}, -1, 1, -1, 1) == false);
    BOOST_TEST(check.isInside({0, 2}, -1, 1, -1, 1) == false);
    BOOST_TEST(check.isInside({2, 0}, -1, 1, -1, 1) == false);
  }
  // Aligned box w/ tolerance check along first axis
  BOOST_AUTO_TEST_CASE(BoundaryCheckBoxToleranceLoc0)
  {
    BoundaryCheck check(true, false, 1.5, 0.0);
    BOOST_TEST(check.isInside({0, 0}, -1, 1, -1, 1) == true);
    BOOST_TEST(check.isInside({2, 2}, -1, 1, -1, 1) == true);
    BOOST_TEST(check.isInside({4, 4}, -1, 1, -1, 1) == false);
    BOOST_TEST(check.isInside({0, 2}, -1, 1, -1, 1) == true);
    BOOST_TEST(check.isInside({2, 0}, -1, 1, -1, 1) == true);
  }
  // Aligned box w/ covariance check
  BOOST_AUTO_TEST_CASE(BoundaryCheckBoxCovariance)
  {
    ActsSymMatrixD<2> cov;
    cov << 1, 0.5, 0.5, 2;
    BoundaryCheck check(cov, 3.0);
    BOOST_TEST(check.isInside({0, 0}, -1, 1, -1, 1) == true);
    BOOST_TEST(check.isInside({2, 2}, -1, 1, -1, 1) == true);
    BOOST_TEST(check.isInside({4, 4}, -1, 1, -1, 1) == false);
    BOOST_TEST(check.isInside({0, 3}, -1, 1, -1, 1) == true);
    BOOST_TEST(check.isInside({3, 0}, -1, 1, -1, 1) == true);
  }
  // Triangle w/ simple check
  BOOST_AUTO_TEST_CASE(BoundaryCheckTriangleSimple)
  {
    Vector2D      vertices[] = {{-2, 0}, {2, 0}, {0, 2}};
    BoundaryCheck check(true);
    BOOST_TEST(check.isInside({0, 0}, vertices) == true);
    BOOST_TEST(check.isInside({0, 1}, vertices) == true);
    BOOST_TEST(check.isInside({2, 2}, vertices) == false);
    BOOST_TEST(check.isInside({0, -1}, vertices) == false);
  }
  // Triangle w/ covariance check
  BOOST_AUTO_TEST_CASE(BoundaryCheckTriangleCovariance)
  {
    Vector2D          vertices[] = {{-2, 0}, {2, 0}, {0, 2}};
    ActsSymMatrixD<2> cov;
    cov << 0.5, 0, 0, 0.5;
    BoundaryCheck check(cov, 4.1);
    BOOST_TEST(check.isInside({0, 0}, vertices) == true);
    BOOST_TEST(check.isInside({0, 1}, vertices) == true);
    BOOST_TEST(check.isInside({0, 2}, vertices) == true);
    BOOST_TEST(check.isInside({0, 3}, vertices) == true);
    BOOST_TEST(check.isInside({0, 4}, vertices) == true);
    BOOST_TEST(check.isInside({0, 5}, vertices) == false);
  }
  BOOST_AUTO_TEST_SUITE_END()
}  // end of namespace Test
}  // end of namespace Acts
