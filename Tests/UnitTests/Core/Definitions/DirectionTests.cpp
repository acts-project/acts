// This file is part of the Acts project.
//
// Copyright (C) 2017-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Direction.hpp"

using namespace Acts;

BOOST_AUTO_TEST_SUITE(DefinitionsDirection)

BOOST_AUTO_TEST_CASE(DirectionTests) {
  constexpr Direction bwd = Direction::Backward;
  constexpr Direction fwd = Direction::Forward;

  BOOST_CHECK(bwd == Direction::Negative);
  BOOST_CHECK(fwd == Direction::Positive);

  BOOST_CHECK(Direction::fromScalar(-1.) == bwd);
  BOOST_CHECK(Direction::fromScalar(1.) == fwd);
  BOOST_CHECK(Direction::fromScalarZeroAsPositive(0) == fwd);

  BOOST_CHECK(Direction::fromIndex(0) == bwd);
  BOOST_CHECK(Direction::fromIndex(1) == fwd);

  BOOST_CHECK(bwd.index() == 0u);
  BOOST_CHECK(fwd.index() == 1u);

  BOOST_CHECK(bwd.sign() == -1);
  BOOST_CHECK(fwd.sign() == +1);

  BOOST_CHECK(bwd.invert() == fwd);
  BOOST_CHECK(fwd.invert() == bwd);

  BOOST_CHECK(bwd.toString() == "backward");
  BOOST_CHECK(fwd.toString() == "forward");

  BOOST_CHECK(2. * fwd == 2.);
  BOOST_CHECK(7 * fwd == 7);
  BOOST_CHECK(Vector3(1., 1., 1.) * fwd == Vector3(1., 1., 1.));

  BOOST_CHECK(2. * bwd == -2.);
  BOOST_CHECK(7 * bwd == -7);
  BOOST_CHECK(Vector3(1., 1., 1.) * bwd == Vector3(-1., -1., -1.));

  double a = 7.;
  a *= fwd;
  BOOST_CHECK(a == 7.);
  a *= bwd;
  BOOST_CHECK(a == -7.);

  float b = 8.;
  b *= fwd;
  BOOST_CHECK(b == 8.);
  b *= bwd;
  BOOST_CHECK(b == -8.);

  Vector3 c(9., 9., 9.);
  c *= fwd;
  BOOST_CHECK(c == Vector3(9., 9., 9.));
  c *= bwd;
  BOOST_CHECK(c == Vector3(-9., -9., -9.));
}

BOOST_AUTO_TEST_SUITE_END()
