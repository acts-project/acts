// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/TrackParameterHelpers.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

using namespace Acts;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(EventDataSuite)

BOOST_AUTO_TEST_CASE(isBoundVectorValiTest) {
  BOOST_CHECK(!isBoundVectorValid({1, 2, 3, 4, 5, 6}, true));
  BOOST_CHECK(isBoundVectorValid({1, 2, 1, 1, 5, 6}, true));
}

BOOST_AUTO_TEST_CASE(isFreeVectorValidTest) {
  BOOST_CHECK(!isFreeVectorValid({1, 2, 3, 4, 5, 6, 7, 8}));
  BOOST_CHECK(isFreeVectorValid({1, 2, 3, 4, 1, 0, 0, 8}));
}

BOOST_AUTO_TEST_CASE(normalizeBoundParametersTest) {
  CHECK_CLOSE_OR_SMALL(normalizeBoundParameters({1, 2, 3, 4, 5, 6}),
                       BoundVector(1, 2, -0.141593, 2.28319, 5, 6), 1e-3, 1e-3);
}

BOOST_AUTO_TEST_CASE(addBoundParametersTest) {
  CHECK_CLOSE_OR_SMALL(
      addBoundParameters({1, 2, 3, 4, 5, 6}, {0, 0, 0, 0, 0, 0}),
      normalizeBoundParameters({1, 2, 3, 4, 5, 6}), 1e-3, 1e-3);
  CHECK_CLOSE_OR_SMALL(
      addBoundParameters({1, 2, 3, 4, 5, 6}, {0, 0, 1, 1, 0, 0}),
      normalizeBoundParameters({1, 2, 4, 5, 5, 6}), 1e-3, 1e-3);
}

BOOST_AUTO_TEST_CASE(subtractBoundParametersTest) {
  CHECK_CLOSE_OR_SMALL(
      subtractBoundParameters({1, 2, 3, 4, 5, 6}, {1, 2, 3, 4, 5, 6}),
      BoundVector(0, 0, 0, 0, 0, 0), 1e-3, 1e-3);
  CHECK_CLOSE_OR_SMALL(
      addBoundParameters(
          subtractBoundParameters({1, 2, 3, 4, 5, 6}, {0, 0, 1, 1, 0, 0}),
          {0, 0, 1, 1, 0, 0}),
      normalizeBoundParameters({1, 2, 3, 4, 5, 6}), 1e-3, 1e-3);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
