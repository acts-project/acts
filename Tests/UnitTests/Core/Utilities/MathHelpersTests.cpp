// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/ArrayHelpers.hpp"
#include "Acts/Utilities/MathHelpers.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <cmath>

namespace bdata = boost::unit_test::data;

const auto expDist = bdata::random(
    (bdata::engine = std::mt19937{}, bdata::seed = 0,
     bdata::distribution = std::uniform_real_distribution<double>(-4, 4)));

namespace ActsTetsts {

BOOST_AUTO_TEST_SUITE(UtilitiesSuite)

BOOST_DATA_TEST_CASE(fastHypot, expDist ^ expDist ^ bdata::xrange(100), xExp,
                     yExp, i) {
  static_cast<void>(i);

  const double x = std::pow(10, xExp);
  const double y = std::pow(10, yExp);

  const float fastFloat =
      Acts::fastHypot(static_cast<float>(x), static_cast<float>(y));
  const double fastDouble = Acts::fastHypot(x, y);

  const float stdFloat =
      std::hypot(static_cast<float>(x), static_cast<float>(y));
  const double stdDouble = std::hypot(x, y);

  CHECK_CLOSE_REL(stdFloat, fastFloat, 1e-6);
  CHECK_CLOSE_REL(stdDouble, fastDouble, 1e-6);
}

BOOST_AUTO_TEST_CASE(ProductTests) {
  static_assert(Acts::product(1u, 1u) == 1u);
  static_assert(Acts::product(2u, 2u) == 2u);
  static_assert(Acts::product(3u, 4u) == 12u);
  static_assert(Acts::product(1u, 5u) == 120u);
  static_assert(Acts::product(2u, 5u) == 120u);
  static_assert(Acts::product(3u, 5u) == 60u);

  static_assert(Acts::product(5u, 4u) == 1u);

  static_assert(Acts::product(0u, 0u) == 0u);
  static_assert(Acts::product(0u, 5u) == 0u);

  // These should fail at compile time
  // static_assert(Acts::product(std::uint8_t{1}, std::uint8_t{6}));
  // static_assert(Acts::product(std::uint16_t{1}, std::uint16_t{9}));
  // static_assert(Acts::product(std::uint32_t{1}, std::uint32_t{13}));
  // static_assert(Acts::product(std::uint64_t{1}, std::uint64_t{21}));

  // These expressions should fail at runtime
  // auto fail8  = Acts::product(std::uint8_t{1},  std::uint8_t{6});
  // auto fail16 = Acts::product(std::uint16_t{1}, std::uint16_t{9});
  // auto fail32 = Acts::product(std::uint32_t{1}, std::uint32_t{13});
  // auto fail64 = Acts::product(std::uint64_t{1}, std::uint64_t{21});

  const std::size_t zero{0};
  const std::size_t one{1};
  for (std::size_t k = 0; k <= 20; ++k) {
    BOOST_CHECK_EQUAL(Acts::product(zero, k), zero);
    BOOST_CHECK_EQUAL(Acts::product(one, k), Acts::factorial(k));
    BOOST_CHECK_EQUAL(Acts::product(k, k), k);

    for (std::size_t j = 1; j < k; ++j) {
      BOOST_CHECK_EQUAL(Acts::product(j, k), j * Acts::product(j + one, k));
      BOOST_CHECK_EQUAL(Acts::product(one, k),
                        Acts::product(one, j) * Acts::product(j + one, k));
    }
  }
}

BOOST_AUTO_TEST_CASE(FactorialTests) {
  static_assert(Acts::factorial(0u) == 1u);
  static_assert(Acts::factorial(1u) == 1u);
  static_assert(Acts::factorial(2u) == 2u);
  static_assert(Acts::factorial(5u) == 5u * Acts::factorial(4u));

  // These tests should fail at compile time
  // static_assert(Acts::factorial(std::uint8_t{6}));
  // static_assert(Acts::factorial(std::uint16_t{9}));
  // static_assert(Acts::factorial(std::uint32_t{13}));
  // static_assert(Acts::factorial(std::uint64_t{21}));

  // These expressions should fail at runtime
  // auto fail8 = Acts::factorial(std::uint8_t{6});
  // auto fail16 = Acts::factorial(std::uint16_t{9});
  // auto fail32 = Acts::factorial(std::uint32_t{13});
  // auto fail64 = Acts::factorial(std::uint64_t{21});

  const std::size_t one{1};
  for (std::size_t k = 1; k <= 20; ++k) {
    BOOST_CHECK_EQUAL(Acts::factorial(k), k * Acts::factorial(k - one));
    for (std::size_t j = 1; j <= k; ++j) {
      BOOST_CHECK_EQUAL(Acts::product(j, k),
                        Acts::factorial(k) / Acts::factorial(j - one));
    }
  }
}

BOOST_AUTO_TEST_CASE(PowerTests) {
  for (unsigned p = 0; p <= 15; ++p) {
    BOOST_CHECK_EQUAL(std::pow(2., p), Acts::pow(2., p));
    BOOST_CHECK_EQUAL(std::pow(0.5, p), Acts::pow(0.5, p));
    for (std::size_t k = 1; k <= 15; ++k) {
      BOOST_CHECK_EQUAL(std::pow(k, p), Acts::pow(k, p));
    }
  }
  for (int p = 0; p <= 15; ++p) {
    BOOST_CHECK_EQUAL(std::pow(2., p), Acts::pow(2., p));
    BOOST_CHECK_EQUAL(std::pow(2., -p), Acts::pow(2., -p));
    BOOST_CHECK_EQUAL(std::pow(0.5, p), Acts::pow(0.5, p));

    BOOST_CHECK_EQUAL(std::pow(0.5, -p), Acts::pow(0.5, -p));
  }
}

BOOST_AUTO_TEST_CASE(SumOfIntegers) {
  std::array<unsigned, 100> numberSeq{Acts::filledArray<unsigned, 100>(1)};
  std::iota(numberSeq.begin(), numberSeq.end(), 1);
  for (unsigned i = 1; i <= numberSeq.size(); ++i) {
    const unsigned sum =
        std::accumulate(numberSeq.begin(), numberSeq.begin() + i, 0);
    BOOST_CHECK_EQUAL(sum, Acts::sumUpToN(i));
  }
}

BOOST_AUTO_TEST_CASE(BinomialTests) {
  static_assert(Acts::binomial(0u, 0u) == 1u);

  static_assert(Acts::binomial(5u, 0u) == 1u);
  static_assert(Acts::binomial(5u, 1u) == 5u);
  static_assert(Acts::binomial(5u, 2u) == 10u);
  static_assert(Acts::binomial(5u, 3u) == 10u);
  static_assert(Acts::binomial(5u, 4u) == 5u);
  static_assert(Acts::binomial(5u, 5u) == 1u);

  static_assert(Acts::binomial(10u, 3u) == 120u);
  static_assert(Acts::binomial(10u, 7u) == 120u);
  static_assert(Acts::binomial(20ull, 10ull) == 184756ull);

  // This test should fail at compile time
  // static_assert(Acts::binomial(4u, 5u));

  BOOST_CHECK_EQUAL(Acts::binomial(1u, 1u), 1u);
  for (unsigned n = 1; n <= 10; ++n) {
    /// Check that the binomial of (n 0 is always 1)
    BOOST_CHECK_EQUAL(Acts::binomial(n, 0u), 1);
    BOOST_CHECK_EQUAL(Acts::binomial(n, n), 1);

    for (unsigned k = 1; k <= n - 1u; ++k) {
      /// Use recursive formula
      ///  n      n -1       n -1
      ///     =          +
      ///  k      k -1        k
      std::cout << "n: " << n << ", \tk: " << k
                << ", \tbinom(n, k): " << Acts::binomial(n, k)
                << ", \tbinom(n-1, k-1): " << Acts::binomial(n - 1, k - 1)
                << ", \tbinom(n-1, k): " << Acts::binomial(n - 1, k)
                << std::endl;
      BOOST_CHECK_EQUAL(Acts::binomial(n, k), Acts::binomial(n - 1, k - 1) +
                                                  Acts::binomial(n - 1, k));
      BOOST_CHECK_EQUAL(Acts::binomial(n, k), Acts::binomial(n, n - k));
    }
  }
}

BOOST_AUTO_TEST_CASE(CopySign) {
  BOOST_CHECK_EQUAL(Acts::copySign(5, -10), -5);
  BOOST_CHECK_EQUAL(Acts::copySign(5, 0), 5);
  BOOST_CHECK_EQUAL(Acts::copySign(5, 55), 5);

  static_assert(Acts::copySign(5, -10) == -5);
  static_assert(Acts::copySign(5, 0) == 5);
  static_assert(Acts::copySign(5, 55) == 5);

  BOOST_CHECK_EQUAL(Acts::copySign(5, -std::numeric_limits<double>::infinity()),
                    -5);
  BOOST_CHECK_EQUAL(Acts::copySign(5, std::numeric_limits<double>::infinity()),
                    5);

  static_assert(Acts::copySign(5., -std::numeric_limits<double>::infinity()) ==
                -5.);
  static_assert(Acts::copySign(5., 0) == 5);
  static_assert(Acts::copySign(5., std::numeric_limits<double>::infinity()) ==
                5.);

  BOOST_CHECK_EQUAL(Acts::copySign(5., -10.), -5.);
  BOOST_CHECK_EQUAL(Acts::copySign(5., 0.), 5.);
  BOOST_CHECK_EQUAL(Acts::copySign(5., 55.), 5.);

  enum class CopyEnum : int { b = 1, a = -1, c = 0 };
  BOOST_CHECK_EQUAL(Acts::copySign(5, CopyEnum::a), -5);
  BOOST_CHECK_EQUAL(Acts::copySign(5, CopyEnum::b), 5);
  BOOST_CHECK_EQUAL(Acts::copySign(5, CopyEnum::c), 5);

  const Acts::Vector3 v{Acts::Vector3::UnitZ()};
  CHECK_CLOSE_ABS(Acts::copySign(v, -1).dot(v), -1., 1.e-7);
  CHECK_CLOSE_ABS(Acts::copySign(v, 1).dot(v), 1., 1.e-7);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTetsts
