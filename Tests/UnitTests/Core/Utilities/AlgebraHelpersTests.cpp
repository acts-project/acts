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
#include "Acts/Utilities/Logger.hpp"

#include <Eigen/Dense>

Acts::Logging::Level logLevel = Acts::Logging::VERBOSE;

namespace Acts::Test {

BOOST_AUTO_TEST_SUITE(AlgebraHelpers)

ACTS_LOCAL_LOGGER(Acts::getDefaultLogger("SafeInverse", logLevel))

BOOST_AUTO_TEST_CASE(SafeInverseSmallMatrix) {
  Eigen::Matrix<double, 2, 2> m;
  m << 1, 2, 3, 4;

  Eigen::Matrix<double, 2, 2> mInvRef;
  mInvRef << -2, 1, 1.5, -0.5;

  auto mInv = Acts::safeInverse(m);

  BOOST_CHECK(mInv);
  BOOST_CHECK_EQUAL(*mInv, mInvRef);

  Eigen::Matrix<int, 2, 2> identityInt;
  identityInt << 1, 0, 0, 1;
  auto identityIntInv = Acts::safeInverse(identityInt);

  BOOST_CHECK(identityIntInv);
  BOOST_CHECK_EQUAL(*identityIntInv, identityInt);
}

BOOST_AUTO_TEST_CASE(safeInverseLargeMatrix) {
  const Eigen::Matrix<double, 5, 5> identity{Eigen::MatrixXd::Identity(5, 5)};
  auto identityInv = Acts::safeInverse(identity);

  BOOST_CHECK(identityInv);
  BOOST_CHECK_EQUAL(*identityInv, identity);
}

BOOST_AUTO_TEST_CASE(SafeInverseBadSmallMatrix) {
  Eigen::Matrix<double, 2, 2> m;
  m << 1, 1, 2, 2;

  auto mInv = Acts::safeInverse(m);

  BOOST_CHECK(!mInv);
}

BOOST_AUTO_TEST_CASE(safeInverseBadLargeMatrix) {
  const Eigen::Matrix<double, 5, 5> m{Eigen::MatrixXd::Zero(5, 5)};
  auto mInv = Acts::safeInverse(m);

  BOOST_CHECK(!mInv);
}

BOOST_AUTO_TEST_CASE(SafeInverseFPESmallMatrix) {
  Eigen::Matrix<double, 4, 4> m = Eigen::MatrixXd::Identity(4, 4) * SIZE_MAX;
  m(1, 1) = 1;

  auto mInv = Acts::safeInverse(m);
  auto mInvInv = Acts::safeInverse(*mInv);

  BOOST_CHECK(mInv);
  BOOST_CHECK(!mInvInv);

  ACTS_VERBOSE("Test: SafeInverseFPESmallMatrix"
               << "\n"
               << "m:\n"
               << m << "\n"
               << "mInv:\n"
               << *mInv << "\n"
               << "mInvInv [garbage]:\n"
               << *mInvInv);
}

BOOST_AUTO_TEST_CASE(SafeInverseFPELargeMatrix) {
  Eigen::Matrix<double, 5, 5> m = Eigen::MatrixXd::Identity(5, 5) * SIZE_MAX;
  m(1, 1) = 1;

  auto mInv = Acts::safeInverse(m);

  BOOST_CHECK(!mInv);

  ACTS_VERBOSE("Test: SafeInverseFPELargeMatrix"
               << "\n"
               << "m:\n"
               << m << "\n"
               << "mInv [garbage]:\n"
               << *mInv);
}

/// This test should not compile
// BOOST_AUTO_TEST_CASE(SafeInverseNonsquareMatrix) {
//   Eigen::Matrix<double, 2, 3> m;
//   m << 1, 2, 3, 4, 5, 6;
//
//   auto mInv = Acts::safeInverse(m);
// }

/// This test should not compile
// BOOST_AUTO_TEST_CASE(SafeInverseDynamicMatrix) {
//   Eigen::MatrixXd m{Eigen::MatrixXd::Identity(2, 2)};
//
//   auto mInv = Acts::safeInverse(m);
// }

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Acts::Test
