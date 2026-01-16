// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/MathHelpers.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <cmath>
#include <cstdint>

namespace bdata = boost::unit_test::data;

const auto expDist = bdata::random(
    (bdata::engine = std::mt19937{}, bdata::seed = 0,
     bdata::distribution = std::uniform_real_distribution<double>(-4, 4)));

namespace ActsTetsts {

BOOST_AUTO_TEST_SUITE(UtilitiesSuite)

BOOST_DATA_TEST_CASE(fastHypot, expDist ^ expDist ^ bdata::xrange(100), xExp,
                     yExp, i) {
  (void)i;

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

BOOST_AUTO_TEST_CASE(Factorial) {
  // Basic factorial tests
  BOOST_CHECK_EQUAL(Acts::factorial(0), 1);
  BOOST_CHECK_EQUAL(Acts::factorial(1), 1);
  BOOST_CHECK_EQUAL(Acts::factorial(5), 120);
  BOOST_CHECK_EQUAL(Acts::factorial(6), 720);

  // Partial factorial tests (n! / (lowerN-1)!)
  BOOST_CHECK_EQUAL(Acts::factorial(5, 3), 60);   // 5 * 4 * 3 = 60
  BOOST_CHECK_EQUAL(Acts::factorial(5, 5), 5);    // just 5
  BOOST_CHECK_EQUAL(Acts::factorial(5, 6), 1);    // lowerN > upperN returns 1

  // Compile-time overflow detection: the following would be a compile error:
  // static_assert(Acts::factorial(std::uint8_t{10}, std::uint8_t{1}) == 0);
  // because 10! overflows uint8_t and triggers the overflow assertion.

  // Verify maximum valid factorial for uint8_t (5! = 120 fits, 6! = 720 doesn't)
  static_assert(Acts::factorial(std::uint8_t{5}, std::uint8_t{1}) == 120);
  BOOST_CHECK_EQUAL(
      Acts::factorial(std::uint8_t{5}, std::uint8_t{1}), std::uint8_t{120});
}

BOOST_AUTO_TEST_CASE(Binomial) {
  // Basic binomial coefficient tests
  BOOST_CHECK_EQUAL(Acts::binomial(5, 0), 1);
  BOOST_CHECK_EQUAL(Acts::binomial(5, 1), 5);
  BOOST_CHECK_EQUAL(Acts::binomial(5, 2), 10);
  BOOST_CHECK_EQUAL(Acts::binomial(5, 5), 1);
  BOOST_CHECK_EQUAL(Acts::binomial(10, 3), 120);

  // Compile-time verification
  static_assert(Acts::binomial(5, 2) == 10);
  static_assert(Acts::binomial(10, 3) == 120);

  // Note: binomial with small integer types that would overflow in factorial
  // now triggers an assertion (compile error if constexpr, runtime assert otherwise)
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
