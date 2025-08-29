// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

#include <cmath>

namespace bdata = boost::unit_test::data;

namespace Acts::Test {

const auto expDist = bdata::random(
    (bdata::engine = std::mt19937{}, bdata::seed = 0,
     bdata::distribution = std::uniform_real_distribution<double>(-4, 4)));

BOOST_AUTO_TEST_SUITE(Utilities)

BOOST_DATA_TEST_CASE(TestFastHypot, expDist ^ expDist ^ bdata::xrange(100),
                     xExp, yExp, i) {
  (void)i;

  const double x = std::pow(10, xExp);
  const double y = std::pow(10, yExp);

  const float fastFloat =
      fastHypot(static_cast<float>(x), static_cast<float>(y));
  const double fastDouble = fastHypot(x, y);

  const float stdFloat =
      std::hypot(static_cast<float>(x), static_cast<float>(y));
  const double stdDouble = std::hypot(x, y);

  CHECK_CLOSE_REL(stdFloat, fastFloat, 1e-6);
  CHECK_CLOSE_REL(stdDouble, fastDouble, 1e-6);
}

BOOST_AUTO_TEST_CASE(TestAbs) {
  // Test signed integer types
  BOOST_CHECK_EQUAL(abs(-5), 5);
  BOOST_CHECK_EQUAL(abs(5), 5);
  BOOST_CHECK_EQUAL(abs(0), 0);
  BOOST_CHECK_EQUAL(abs(-1), 1);

  // Test unsigned integer types
  BOOST_CHECK_EQUAL(abs(5u), 5u);
  BOOST_CHECK_EQUAL(abs(0u), 0u);

  // Test floating point types
  BOOST_CHECK_EQUAL(abs(-3.14), 3.14);
  BOOST_CHECK_EQUAL(abs(3.14), 3.14);
  BOOST_CHECK_EQUAL(abs(0.0), 0.0);
  BOOST_CHECK_EQUAL(abs(-0.0), 0.0);

  // Edge cases: lowest signed wraps around to max int value
  BOOST_CHECK(abs(std::numeric_limits<int>::lowest()) ==
              std::numeric_limits<int>::max());

  BOOST_CHECK_EQUAL(
      abs(std::numeric_limits<int>::lowest() + 1),
      static_cast<unsigned int>(std::numeric_limits<int>::lowest()) - 1);

  // Test different types
  BOOST_CHECK_EQUAL(abs(static_cast<short>(-10)), static_cast<short>(10));
  BOOST_CHECK_EQUAL(abs(-10L), 10L);
  BOOST_CHECK_EQUAL(abs(-10.0f), 10.0f);
}

BOOST_AUTO_TEST_CASE(TestPow) {
  // Basic power calculations
  BOOST_CHECK_EQUAL(pow(2, 3u), 8);
  BOOST_CHECK_EQUAL(pow(5, 2u), 25);
  BOOST_CHECK_EQUAL(pow(3, 1u), 3);

  // Zero power should return 1
  BOOST_CHECK_EQUAL(pow(5, 0u), 1);
  BOOST_CHECK_EQUAL(pow(0, 0u), 1);
  BOOST_CHECK_EQUAL(pow(123, 0u), 1);

  // Zero base with positive power
  BOOST_CHECK_EQUAL(pow(0, 1u), 0);
  BOOST_CHECK_EQUAL(pow(0, 5u), 0);

  // Different type combinations
  BOOST_CHECK_EQUAL(pow(2.0, 3u), 8.0);
  BOOST_CHECK_EQUAL(pow(2, 3UL), 8);
  BOOST_CHECK_EQUAL(pow(2.0f, 2u), 4.0f);

  // Test with unsigned power type
  BOOST_CHECK_EQUAL(pow(3, 2u), 9);
  BOOST_CHECK_EQUAL(pow(2, 4u), 16);
}

BOOST_AUTO_TEST_CASE(TestSquare) {
  // Integer types
  BOOST_CHECK_EQUAL(square(5), 25);
  BOOST_CHECK_EQUAL(square(-3), 9);
  BOOST_CHECK_EQUAL(square(0), 0);

  // Floating point types
  BOOST_CHECK_EQUAL(square(3.0), 9.0);
  BOOST_CHECK_EQUAL(square(-2.5f), 6.25f);
  BOOST_CHECK_EQUAL(square(0.0), 0.0);

  // Different numeric types
  BOOST_CHECK_EQUAL(square(4L), 16L);
  BOOST_CHECK_EQUAL(square(static_cast<short>(3)), static_cast<short>(9));
  BOOST_CHECK_EQUAL(square(2u), 4u);
}

BOOST_AUTO_TEST_CASE(TestHypotSquare) {
  // Single argument
  BOOST_CHECK_EQUAL(hypotSquare(3), 9);
  BOOST_CHECK_EQUAL(hypotSquare(4.0), 16.0);

  // Two arguments (classic Pythagorean)
  BOOST_CHECK_EQUAL(hypotSquare(3, 4), 25);
  BOOST_CHECK_EQUAL(hypotSquare(-3, 4), 25);
  BOOST_CHECK_EQUAL(hypotSquare(0, 5), 25);

  // Three arguments
  BOOST_CHECK_EQUAL(hypotSquare(1, 2, 3), 14);
  BOOST_CHECK_EQUAL(hypotSquare(2, 3, 6), 49);

  // Mixed types
  BOOST_CHECK_EQUAL(hypotSquare(3, 4.0f), 25.0f);
  BOOST_CHECK_EQUAL(hypotSquare(2.0, 3), 13.0);

  // More arguments
  BOOST_CHECK_EQUAL(hypotSquare(1, 1, 1, 1), 4);
}

BOOST_AUTO_TEST_CASE(TestFastHypotExtended) {
  // Single argument
  CHECK_CLOSE_REL(fastHypot(3.0), 3.0, 1e-10);
  CHECK_CLOSE_REL(fastHypot(4.0f), 4.0f, 1e-6);

  // Two arguments (already tested in TestFastHypot, but add more cases)
  CHECK_CLOSE_REL(fastHypot(3, 4), 5.0, 1e-10);
  CHECK_CLOSE_REL(fastHypot(0.0, 5.0), 5.0, 1e-10);

  // Three arguments
  CHECK_CLOSE_REL(fastHypot(1.0, 2.0, 3.0), std::sqrt(14.0), 1e-10);
  CHECK_CLOSE_REL(fastHypot(2.0, 3.0, 6.0), 7.0, 1e-10);

  // Mixed types
  CHECK_CLOSE_REL(fastHypot(3.0f, 4.0), 5.0, 1e-6);

  // More arguments
  CHECK_CLOSE_REL(fastHypot(1.0, 1.0, 1.0, 1.0), 2.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(TestSumUpToN) {
  // Basic cases
  BOOST_CHECK_EQUAL(sumUpToN(0), 0);
  BOOST_CHECK_EQUAL(sumUpToN(1), 1);
  BOOST_CHECK_EQUAL(sumUpToN(2), 3);
  BOOST_CHECK_EQUAL(sumUpToN(3), 6);
  BOOST_CHECK_EQUAL(sumUpToN(5), 15);
  BOOST_CHECK_EQUAL(sumUpToN(10), 55);

  // Negative input (mathematical edge case)
  BOOST_CHECK_EQUAL(sumUpToN(-1), 0);  // N * (N + 1) / 2 = -1 * 0 / 2 = 0
  BOOST_CHECK_EQUAL(sumUpToN(-5), 10);
  BOOST_CHECK_EQUAL(sumUpToN(-10), 45);

  // Different integral types
  BOOST_CHECK_EQUAL(sumUpToN(5L), 15L);
  BOOST_CHECK_EQUAL(sumUpToN(5u), 15u);
  BOOST_CHECK_EQUAL(sumUpToN(static_cast<short>(4)), static_cast<short>(10));
}

BOOST_AUTO_TEST_CASE(TestFactorial) {
  // Basic factorial cases
  BOOST_CHECK_EQUAL(factorial(0), 1);
  BOOST_CHECK_EQUAL(factorial(1), 1);
  BOOST_CHECK_EQUAL(factorial(2), 2);
  BOOST_CHECK_EQUAL(factorial(3), 6);
  BOOST_CHECK_EQUAL(factorial(4), 24);
  BOOST_CHECK_EQUAL(factorial(5), 120);

  // Factorial with custom lower bound
  BOOST_CHECK_EQUAL(factorial(5, 3), 60);   // 5*4*3 = 60
  BOOST_CHECK_EQUAL(factorial(6, 4), 120);  // 6*5*4 = 120
  BOOST_CHECK_EQUAL(factorial(4, 4), 4);    // Just 4
  BOOST_CHECK_EQUAL(factorial(3, 5), 1);    // upperN < lowerN returns 1

  // Test with different integral types
  BOOST_CHECK_EQUAL(factorial<int>(4), 24);
  BOOST_CHECK_EQUAL(factorial<long>(4), 24L);
  BOOST_CHECK_EQUAL(factorial<unsigned int>(4u), 24u);

  // Test lowerN = 0 (should be treated as 1 due to std::max)
  BOOST_CHECK_EQUAL(factorial(5, 0), 120);  // Same as factorial(5, 1)
  BOOST_CHECK_EQUAL(factorial(3, 0), 6);    // Same as factorial(3, 1)
}

BOOST_AUTO_TEST_CASE(TestBinomial) {
  // Basic binomial coefficient cases
  BOOST_CHECK_EQUAL(binomial(5, 0), 1);   // C(5,0) = 1
  BOOST_CHECK_EQUAL(binomial(5, 1), 5);   // C(5,1) = 5
  BOOST_CHECK_EQUAL(binomial(5, 2), 10);  // C(5,2) = 10
  BOOST_CHECK_EQUAL(binomial(5, 3), 10);  // C(5,3) = 10
  BOOST_CHECK_EQUAL(binomial(5, 4), 5);   // C(5,4) = 5
  BOOST_CHECK_EQUAL(binomial(5, 5), 1);   // C(5,5) = 1

  // Larger binomial coefficients
  BOOST_CHECK_EQUAL(binomial(10, 3), 120);  // C(10,3) = 120
  BOOST_CHECK_EQUAL(binomial(8, 4), 70);    // C(8,4) = 70

  // Edge case: n = 0
  BOOST_CHECK_EQUAL(binomial(0, 0), 1);

  // Test with different integral types
  BOOST_CHECK_EQUAL(binomial<int>(6, 2), 15);
  BOOST_CHECK_EQUAL(binomial<long>(6L, 2L), 15L);
  BOOST_CHECK_EQUAL(binomial<unsigned int>(6u, 2u), 15u);

  BOOST_CHECK_EQUAL(binomial(3, 5), 0);
}

// Tests for undefined behavior and edge cases
BOOST_AUTO_TEST_CASE(TestUndefinedBehaviorCases) {
  // Test pow function division by zero cases
  // These should be handled gracefully or at least not crash

  // Let's test what happens with very small values instead of zero
  auto tiny = std::numeric_limits<double>::min();
  auto result_tiny = pow(tiny, 2u);
  BOOST_CHECK(std::isfinite(result_tiny));
  BOOST_CHECK_EQUAL(result_tiny,
                    0);  // This comes out to be zero due to multiplication

  // Test square overflow cases
  // For 32-bit int, max value is 2,147,483,647
  // sqrt(2,147,483,647) ≈ 46340.95
  int near_sqrt_max = 46340;
  BOOST_CHECK_GT(square(near_sqrt_max), 0);  // Should still be positive

  // Overflows for 32-bit int:
  int overflow_val = 46341;
  // Causes signed overflow => UB
  auto overflow_result = square(overflow_val);
  BOOST_CHECK_LT(overflow_result, 0);

  // Test sumUpToN overflow cases
  // For 32-bit int, N * (N + 1) overflows when N >= 46341
  // Test values just below the overflow threshold
  int safe_n = 46340;
  BOOST_CHECK_EQUAL(sumUpToN(safe_n), 1073720970);  // Should still be positive

  // This would overflow:
  int overflow_n = 46341;
  auto overflow_sum = sumUpToN(overflow_n);  // Would cause signed overflow UB
  BOOST_CHECK_LT(overflow_sum, 0);

  // Test factorial overflow cases
  // Factorial grows very fast:
  // 12! = 479,001,600 (fits in 32-bit int)
  // 13! = 6,227,020,800 (overflows 32-bit int)
  BOOST_CHECK_EQUAL(factorial(12), 479001600);

  // This overflows for 32-bit int => UB
  // auto overflow_factorial = factorial(13);
  // Same for 32-bit unsigned int
  // auto overflow_factorial = factorial(13u);

  // Test with 64-bit unsigned type to avoid signed overflow UB
  BOOST_CHECK_EQUAL(factorial(13ul), 6227020800);

  // Test hypotSquare overflow scenarios
  // Test values that individually are safe but together overflow
  int safe_individual = 30000;  // 30000^2 = 900,000,000 (fits in int)
  BOOST_CHECK_GT(square(safe_individual), 0);

  // But multiple large values together would overflow:
  auto overflow_hypot = hypotSquare(46000, 46000);
  BOOST_CHECK_LT(overflow_hypot, 0);  // Signed overflow => UB

  // Test with smaller values that are definitely safe
  BOOST_CHECK_EQUAL(hypotSquare(20000, 20000), 800000000);

  // Also works with 64-bit ints
  auto non_overflow_hypot = hypotSquare(46000ul, 46000ul);
  BOOST_CHECK_EQUAL(non_overflow_hypot, 4232000000);
}

// Test floating point edge cases
BOOST_AUTO_TEST_CASE(TestFloatingPointEdgeCases) {
  // Test with infinity
  double inf = std::numeric_limits<double>::infinity();
  BOOST_CHECK(std::isinf(square(inf)));
  BOOST_CHECK(std::isinf(fastHypot(inf, 1.0)));

  // Test with NaN
  double nan_val = std::numeric_limits<double>::quiet_NaN();
  BOOST_CHECK(std::isnan(square(nan_val)));
  BOOST_CHECK(std::isnan(fastHypot(nan_val, 1.0)));

  // Test very large floating point values that don't overflow
  double large = std::sqrt(std::numeric_limits<double>::max());
  BOOST_CHECK(std::isfinite(square(large)));  // Should not be infinite

  // Test very small values
  double tiny = std::numeric_limits<double>::min();
  // These are ZERO due to floating point arithmetic
  BOOST_CHECK_EQUAL(square(tiny), 0.0);
  BOOST_CHECK_EQUAL(fastHypot(tiny, tiny), 0.0);
}

// Test to demonstrate potential stack overflow issues
BOOST_AUTO_TEST_CASE(TestStackOverflowRisks) {
  // The factorial function uses recursion which can cause stack overflow
  // for large inputs. We test with progressively larger values.

  // These should be safe on most systems
  BOOST_CHECK_EQUAL(factorial(20ul), 2432902008176640000);
  // This overflows on 64-bits, but does not cause stack overflow
  factorial(50ul);

  // WARNING: These larger values could cause stack overflow
  // The exact limit depends on stack size, but typically:
  // - factorial(1000) might cause issues
  // - factorial(10000) would almost certainly cause stack overflow
  // DANGER: Do not uncomment these - they could crash the test
  // auto dangerous1 = factorial(1000);   // Likely stack overflow
  // auto dangerous2 = factorial(10000);  // Almost certain stack overflow

  // Test with a moderately large value that should still be safe
  // but demonstrates the recursive nature
  auto safe_large = factorial(100);
  static_cast<void>(safe_large);  // Avoid unused variable warning
  // We don't check the value since it will overflow, but we check it doesn't
  // crash
  BOOST_CHECK(true);  // If we get here, no stack overflow occurred
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Acts::Test
