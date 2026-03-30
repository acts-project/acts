// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/AlgebraHelpers.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <Eigen/Dense>

using namespace Acts;

Logging::Level logLevel = Logging::VERBOSE;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(UtilitiesSuite)

BOOST_AUTO_TEST_SUITE(SafeInverse)

ACTS_LOCAL_LOGGER(getDefaultLogger("SafeInverse", logLevel))

BOOST_AUTO_TEST_CASE(SafeInverseSmallMatrix) {
  Eigen::Matrix<double, 2, 2> m;
  m << 1, 2, 3, 4;

  Eigen::Matrix<double, 2, 2> mInvRef;
  mInvRef << -2, 1, 1.5, -0.5;

  auto mInv = safeInverse(m);

  BOOST_CHECK(mInv);
  BOOST_CHECK_EQUAL(*mInv, mInvRef);

  Eigen::Matrix<int, 2, 2> identityInt;
  identityInt << 1, 0, 0, 1;
  auto identityIntInv = safeInverse(identityInt);

  BOOST_CHECK(identityIntInv);
  BOOST_CHECK_EQUAL(*identityIntInv, identityInt);
}

BOOST_AUTO_TEST_CASE(safeInverseLargeMatrix) {
  const Eigen::Matrix<double, 5, 5> identity{Eigen::MatrixXd::Identity(5, 5)};
  auto identityInv = safeInverse(identity);

  BOOST_CHECK(identityInv);
  BOOST_CHECK_EQUAL(*identityInv, identity);
}

BOOST_AUTO_TEST_CASE(safeInverseDynamicMatrix) {
  Eigen::MatrixXd identity{Eigen::MatrixXd::Identity(2, 2)};

  auto identityInv = safeInverse(identity);

  BOOST_CHECK(identityInv);
  BOOST_CHECK_EQUAL(*identityInv, identity);
}

BOOST_AUTO_TEST_CASE(SafeInverseBadSmallMatrix) {
  Eigen::Matrix<double, 2, 2> m;
  m << 1, 1, 2, 2;

  auto mInv = safeInverse(m);

  BOOST_CHECK(!mInv);
}

BOOST_AUTO_TEST_CASE(safeInverseBadLargeMatrix) {
  const Eigen::Matrix<double, 5, 5> m{Eigen::MatrixXd::Zero(5, 5)};
  auto mInv = safeInverse(m);

  BOOST_CHECK(!mInv);
}

BOOST_AUTO_TEST_CASE(SafeInverseFPESmallMatrix) {
  Eigen::Matrix<double, 4, 4> m =
      Eigen::MatrixXd::Identity(4, 4) * std::numeric_limits<std::size_t>::max();
  m(1, 1) = 1;

  auto mInv = safeInverse(m);
  BOOST_REQUIRE(mInv.has_value());
  auto mInvInv = safeInverse(*mInv);
  BOOST_CHECK(!mInvInv);

  ACTS_VERBOSE("Test: SafeInverseFPESmallMatrix" << "\n"
                                                 << "m:\n"
                                                 << m << "\n"
                                                 << "mInv:\n"
                                                 << *mInv);
}

BOOST_AUTO_TEST_CASE(SafeInverseFPELargeMatrix) {
  Eigen::Matrix<double, 5, 5> m =
      Eigen::MatrixXd::Identity(5, 5) * std::numeric_limits<std::size_t>::max();
  m(1, 1) = 1;

  auto mInv = safeInverse(m);

  BOOST_CHECK(!mInv);

  ACTS_VERBOSE("Test: SafeInverseFPELargeMatrix" << "\n"
                                                 << "m:\n"
                                                 << m);
}

/// This test should not compile
// BOOST_AUTO_TEST_CASE(SafeInverseNonsquareMatrix) {
//   Eigen::Matrix<double, 2, 3> m;
//   m << 1, 2, 3, 4, 5, 6;
//
//   auto mInv = safeInverse(m);
// }

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SafeExp)

ACTS_LOCAL_LOGGER(getDefaultLogger("SafeExp", logLevel))

BOOST_AUTO_TEST_CASE(safeExpDouble) {
  using FloatType = double;

  // Values within the safe range
  BOOST_CHECK_CLOSE(safeExp<FloatType>(0.0), std::exp(0.0), 1e-8);
  BOOST_CHECK_CLOSE(safeExp<FloatType>(1.0), std::exp(1.0), 1e-8);
  BOOST_CHECK_CLOSE(safeExp<FloatType>(-1.0), std::exp(-1.0), 1e-8);

  // Values causing underflow
  BOOST_CHECK_EQUAL(safeExp<FloatType>(-600.0), 0.0);

  // Values causing overflow
  BOOST_CHECK_EQUAL(safeExp<FloatType>(600.0),
                    std::numeric_limits<FloatType>::infinity());
}

BOOST_AUTO_TEST_CASE(safeExpFloat) {
  using FloatType = float;

  // Values within the safe range
  BOOST_CHECK_CLOSE(safeExp<FloatType>(0.0f), std::exp(0.0f), 1e-8);
  BOOST_CHECK_CLOSE(safeExp<FloatType>(1.0f), std::exp(1.0f), 1e-8);
  BOOST_CHECK_CLOSE(safeExp<FloatType>(-1.0f), std::exp(-1.0f), 1e-8);

  // Values causing underflow
  BOOST_CHECK_EQUAL(safeExp<FloatType>(-60.0f), 0.0f);

  // Values causing overflow
  BOOST_CHECK_EQUAL(safeExp<FloatType>(60.0f),
                    std::numeric_limits<FloatType>::infinity());
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(MatrixIdxUnrolling)

template <std::size_t N>
bool testUnrolling()
  requires(N > 1)
{
  /// Fill a test matrix with unique values
  SquareMatrix<N> testMatrix{SquareMatrix<N>::Zero()};
  std::size_t counter{N};
  for (std::size_t i = 0; i < N; ++i) {
    for (std::size_t j = 0; j <= i; ++j) {
      testMatrix(j, i) = testMatrix(i, j) = counter;
      ++counter;
    }
  }
  /// Then unroll the matrix into a vector
  std::vector<double> vec{};
  for (std::size_t i = 0; i < N; ++i) {
    for (std::size_t j = 0; j <= i; ++j) {
      vec.push_back(testMatrix(i, j));
    }
  }

  if (vec.size() != sumUpToN(N)) {
    std::cout << "Compressed vector size mismatch: expected " << sumUpToN(N)
              << ", got " << vec.size() << std::endl;
    return false;
  }
  /// Next test whether the indices are correct
  for (std::size_t i = 0; i < N; ++i) {
    for (std::size_t j = 0; j <= i; ++j) {
      const auto idx = vecIdxFromSymMat<N>(i, j);
      if (idx >= vec.size()) {
        std::cout << "Index out of bounds: " << idx << " for size "
                  << vec.size() << std::endl;
        return false;
      }
      if (idx != vecIdxFromSymMat<N>(j, i)) {
        std::cout << "Index mismatch: expected " << vecIdxFromSymMat<N>(i, j)
                  << ", got " << idx << std::endl;
        return false;
      }
      if (vec[idx] != testMatrix(i, j)) {
        std::cout << "Value mismatch at index " << idx << ": expected "
                  << testMatrix(i, j) << ", got " << vec[idx] << std::endl;
        return false;
      }
      /// Finally, check whether the indices can be recovered
      const auto [iBack, jBack] = symMatIndices<N>(idx);
      if (iBack != i || jBack != j) {
        std::cout << "Index mismatch: expected (" << i << ", " << j
                  << "), got (" << iBack << ", " << jBack << ") for index "
                  << idx << std::endl;
        return false;
      }
    }
  }
  if constexpr (N > 2) {
    return testUnrolling<N - 1>();
  }
  return true;
}
BOOST_AUTO_TEST_CASE(matrixIdxUnrolling) {
  BOOST_CHECK_EQUAL(testUnrolling<20>(), true);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
