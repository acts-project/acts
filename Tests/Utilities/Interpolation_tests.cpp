// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE interpolation tests
#include <boost/test/included/unit_test.hpp>

#include "ACTS/Utilities/Interpolation.hpp"
#include "ACTS/Utilities/detail/interpolation_impl.hpp"

namespace Acts {

using namespace detail;

namespace Test {

  BOOST_AUTO_TEST_CASE(interpolation_1d)
  {
    typedef std::array<double, 1u> Point;
    typedef std::array<double, 2u> Values;

    Point  low  = {1.};
    Point  high = {2.};
    Values v    = {10., 20.};

    BOOST_TEST(interpolate(Point({0.5}), low, high, v) == 5.);
    BOOST_TEST(interpolate(Point({1.}), low, high, v) == 10.);
    BOOST_TEST(interpolate(Point({1.3}), low, high, v) == 13.);
    BOOST_TEST(interpolate(Point({1.5}), low, high, v) == 15.);
    BOOST_TEST(interpolate(Point({1.8}), low, high, v) == 18.);
    BOOST_TEST(interpolate(Point({2.}), low, high, v) == 20.);
    BOOST_TEST(interpolate(Point({2.3}), low, high, v) == 23.);
  }

  BOOST_AUTO_TEST_CASE(interpolation_2d)
  {
    typedef std::array<double, 2u> Point;
    typedef std::array<double, 4u> Values;

    Point  low  = {1., 1.};
    Point  high = {2., 3.};
    Values v    = {10., 20., 30., 40.};

    BOOST_TEST(interpolate(Point({1., 1.}), low, high, v) == 10.);
    BOOST_TEST(interpolate(Point({2., 1.}), low, high, v) == 20.);
    BOOST_TEST(interpolate(Point({1., 3.}), low, high, v) == 30.);
    BOOST_TEST(interpolate(Point({2., 3.}), low, high, v) == 40.);
    BOOST_TEST(interpolate(Point({1.3, 1.}), low, high, v) == 13.);
    BOOST_TEST(interpolate(Point({1.5, 1.}), low, high, v) == 15.);
    BOOST_TEST(interpolate(Point({1.8, 1.}), low, high, v) == 18.);
    BOOST_TEST(interpolate(Point({1.3, 3.}), low, high, v) == 33.);
    BOOST_TEST(interpolate(Point({1.5, 3.}), low, high, v) == 35.);
    BOOST_TEST(interpolate(Point({1.8, 3.}), low, high, v) == 38.);
    BOOST_TEST(interpolate(Point({1., 1.7}), low, high, v) == 17.);
    BOOST_TEST(interpolate(Point({1., 2.}), low, high, v) == 20.);
    BOOST_TEST(interpolate(Point({1., 2.5}), low, high, v) == 25.);
    BOOST_TEST(interpolate(Point({2., 1.7}), low, high, v) == 27.);
    BOOST_TEST(interpolate(Point({2., 2.}), low, high, v) == 30.);
    BOOST_TEST(interpolate(Point({2., 2.5}), low, high, v) == 35.);
    BOOST_TEST(interpolate(Point({1.5, 2.}), low, high, v) == 25.);
    BOOST_TEST(interpolate(Point({1.3, 1.7}), low, high, v) == 20.);
    BOOST_TEST(interpolate(Point({1.3, 2.5}), low, high, v) == 28.);
    BOOST_TEST(interpolate(Point({1.8, 1.7}), low, high, v) == 25.);
    BOOST_TEST(interpolate(Point({1.8, 2.5}), low, high, v) == 33.);
  }

  BOOST_AUTO_TEST_CASE(interpolation_3d)
  {
    typedef std::array<double, 3u> Point;
    typedef std::array<double, 8u> Values;

    Point  low  = {1., 1., 1.};
    Point  high = {2., 3., 4.};
    Values v    = {10., 20., 30., 40., 50., 60., 70., 80.};

    BOOST_TEST(interpolate(Point({1., 1., 1.}), low, high, v) == 10.);
    BOOST_TEST(interpolate(Point({2., 1., 1.}), low, high, v) == 20.);
    BOOST_TEST(interpolate(Point({1., 3., 1.}), low, high, v) == 30.);
    BOOST_TEST(interpolate(Point({2., 3., 1.}), low, high, v) == 40.);
    BOOST_TEST(interpolate(Point({1., 1., 4.}), low, high, v) == 50.);
    BOOST_TEST(interpolate(Point({2., 1., 4.}), low, high, v) == 60.);
    BOOST_TEST(interpolate(Point({1., 3., 4.}), low, high, v) == 70.);
    BOOST_TEST(interpolate(Point({2., 3., 4.}), low, high, v) == 80.);
    BOOST_TEST(interpolate(Point({1.5, 1., 1.}), low, high, v) == 15.);
    BOOST_TEST(interpolate(Point({1.5, 3., 1.}), low, high, v) == 35.);
    BOOST_TEST(interpolate(Point({1., 2., 1.}), low, high, v) == 20.);
    BOOST_TEST(interpolate(Point({2., 2., 1.}), low, high, v) == 30.);
    BOOST_TEST(interpolate(Point({1.5, 1., 4.}), low, high, v) == 55.);
    BOOST_TEST(interpolate(Point({1.5, 3., 4.}), low, high, v) == 75.);
    BOOST_TEST(interpolate(Point({1., 2., 4.}), low, high, v) == 60.);
    BOOST_TEST(interpolate(Point({2., 2., 4.}), low, high, v) == 70.);
    BOOST_TEST(interpolate(Point({1., 1., 2.5}), low, high, v) == 30.);
    BOOST_TEST(interpolate(Point({1., 3., 2.5}), low, high, v) == 50.);
    BOOST_TEST(interpolate(Point({2., 1., 2.5}), low, high, v) == 40.);
    BOOST_TEST(interpolate(Point({2., 3., 2.5}), low, high, v) == 60.);
    BOOST_TEST(interpolate(Point({1.5, 2., 2.5}), low, high, v) == 360. / 8);
    BOOST_TEST(interpolate(Point({1.3, 2.1, 1.6}), low, high, v) == 32.);
  }
}  // namespace Test

}  // namespace Acts
