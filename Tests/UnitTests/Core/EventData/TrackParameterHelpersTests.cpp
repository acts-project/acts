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

BOOST_AUTO_TEST_CASE(CanonicalDirectionMapsAndExpressions) {
  static_assert(eBoundDir0 == eBoundPhi);
  static_assert(eBoundDir1 == eBoundTheta);
  static_assert(eBoundSize == 6);
  static_assert(eFreeSize == 8);
  static_assert(sizeof(BoundVector) == 6 * sizeof(double));
  static_assert(sizeof(BoundMatrix) == 36 * sizeof(double));
  BoundMatrix covariance = BoundMatrix::Zero();
  Eigen::Map<BoundMatrix> covarianceMap(covariance.data());
  covarianceMap(eBoundDir0, eBoundDir1) = 0.25;
  BOOST_CHECK_EQUAL(covariance(eBoundPhi, eBoundTheta), 0.25);
  BoundVector parameters = BoundVector::Zero();
  Eigen::Map<BoundVector> map(parameters.data());
  map[eBoundDir0] = 0.4;
  map[eBoundDir1] = 1.2;
  const Eigen::Map<const BoundVector> constMap(parameters.data());
  const Vector3 expected = makeDirectionFromPhiTheta(0.4, 1.2);
  BOOST_CHECK_EQUAL(directionFromBoundParameters(map), expected);
  BOOST_CHECK_EQUAL(directionFromBoundParameters(constMap), expected);
  BOOST_CHECK_EQUAL(
      directionFromBoundParameters(parameters + BoundVector::Zero()), expected);
  BOOST_CHECK_EQUAL(parameters[eBoundPhi], 0.4);
  BOOST_CHECK_EQUAL(parameters[eBoundTheta], 1.2);
  map[eBoundDir1] = 0.;
  BOOST_CHECK_EQUAL(directionFromBoundParameters(map), Vector3::UnitZ());
}

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
