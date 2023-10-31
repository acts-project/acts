// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/AlgebraHelpers.hpp"

#include <Eigen/Dense>

namespace Acts {
namespace Test {

BOOST_AUTO_TEST_SUITE(AlgebraHelpers)

BOOST_AUTO_TEST_CASE(safeInverseSmallMatrix) {
  Eigen::Matrix<double, 2, 2> aDouble;
  aDouble << 1, 2, 3, 4;

  Eigen::Matrix<double, 2, 2> aDoubleInvRef;
  aDoubleInvRef << -2, 1, 1.5, -0.5;

  auto aDoubleInv = Acts::safeInverse(aDouble);

  BOOST_CHECK(aDoubleInv);
  BOOST_CHECK_EQUAL(*aDoubleInv, aDoubleInvRef);

  Eigen::Matrix<int, 2, 2> identityInt;
  identityInt << 1, 0, 0, 1;
  auto identityIntInv = Acts::safeInverse(identityInt);

  BOOST_CHECK(identityIntInv);
  BOOST_CHECK_EQUAL(*identityIntInv, identityInt);
}

BOOST_AUTO_TEST_CASE(safeInverseBadSmallMatrix) {
  Eigen::Matrix<double, 2, 2> aDouble;
  aDouble << 1, 1, 2, 2;

  auto aDoubleInv = Acts::safeInverse(aDouble);

  BOOST_CHECK(!aDoubleInv);
}

BOOST_AUTO_TEST_CASE(safeInverseLargeMatrix) {
  const Eigen::Matrix<double, 5, 5> identity{Eigen::MatrixXd::Identity(5, 5)};
  auto identityInv = Acts::safeInverse(identity);

  BOOST_CHECK(identityInv);
  BOOST_CHECK_EQUAL(*identityInv, identity);
}

BOOST_AUTO_TEST_CASE(safeInverseBadLargeMatrix) {
  const Eigen::Matrix<double, 5, 5> aDouble{Eigen::MatrixXd::Zero(5, 5)};
  auto aDoubleInv = Acts::safeInverse(aDouble);

  BOOST_CHECK(!aDoubleInv);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test
}  // namespace Acts
