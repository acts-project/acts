// This file is part of the Acts project.
//
// Copyright (C) 2017-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"

using namespace Acts;

BOOST_AUTO_TEST_SUITE(CommonDefinitions)

BOOST_AUTO_TEST_CASE(NavigationDirectionTests) {
  BOOST_CHECK(indexFromDirection(NavigationDirection::Backward) == 0u);
  BOOST_CHECK(indexFromDirection(NavigationDirection::Forward) == 1u);

  BOOST_CHECK(directionFromStepSize(-1.) == NavigationDirection::Backward);
  BOOST_CHECK(directionFromStepSize(1.) == NavigationDirection::Forward);

  BOOST_CHECK(invertDirection(NavigationDirection::Backward) ==
              NavigationDirection::Forward);
  BOOST_CHECK(invertDirection(NavigationDirection::Forward) ==
              NavigationDirection::Backward);

  NavigationDirection nFwd = NavigationDirection::Forward;
  NavigationDirection nBwd = NavigationDirection::Backward;

  BOOST_CHECK(2. * nFwd == 2.);
  BOOST_CHECK(7 * nFwd == 7);
  BOOST_CHECK(Vector3(1., 1., 1.) * nFwd == Vector3(1., 1., 1.));

  BOOST_CHECK(2. * nBwd == -2.);
  BOOST_CHECK(7 * nBwd == -7);
  BOOST_CHECK(Vector3(1., 1., 1.) * nBwd == Vector3(-1., -1., -1.));

  double a = 7.;
  a *= nFwd;
  BOOST_CHECK(a == 7.);
  a *= nBwd;
  BOOST_CHECK(a == -7.);

  float b = 8.;
  b *= nFwd;
  BOOST_CHECK(b == 8.);
  b *= nBwd;
  BOOST_CHECK(b == -8.);

  Vector3 c(9., 9., 9.);
  c *= nFwd;
  BOOST_CHECK(c == Vector3(9., 9., 9.));
  c *= nBwd;
  BOOST_CHECK(c == Vector3(-9., -9., -9.));
}

BOOST_AUTO_TEST_SUITE_END()
