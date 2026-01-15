// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cassert>
#include <cmath>
#include <type_traits>
namespace Acts {

/// @brief Returns the absolute of a number
///        (Can be removed for c++ 23)
/// @param n The number to take absolute value of
/// @return The absolute value of the input
template <typename T>
constexpr T abs(const T n) {
  if (std::is_constant_evaluated()) {
    if constexpr (std::is_signed_v<T>) {
      if (n < 0) {
        return -n;
      }
    }
    return n;
  } else {
    return std::abs(n);
  }
}
/// @brief Copies the sign of a signed variable onto the copyTo input object
///        Return type & magnitude remain unaffected by this method which allows
///        usage for Vectors & other types providing the - operator.
///        By convention, the zero is assigned to a positive sign.
/// @param copyTo: Variable to which the sign is copied to.
/// @param sign: Variable from which the sign is taken.
template <typename out_t, typename sign_t>
constexpr out_t copySign(const out_t& copyTo, const sign_t& sign) {
  if constexpr (std::is_enum_v<sign_t>) {
    return copySign(copyTo, static_cast<std::underlying_type_t<sign_t>>(sign));
  } else {
    constexpr sign_t zero = 0;
    return sign >= zero ? copyTo : -copyTo;
  }
}

/// @brief Calculates the ordinary power of the number x.
/// @param x: Number to take the power from
/// @param p: Power to take
/// @return x raised to the power p
template <typename T, std::integral P>
constexpr T pow(T x, P p) {
  if (std::is_constant_evaluated()) {
    constexpr T one = 1;
    if constexpr (std::is_signed_v<P>) {
      if (p < 0 && abs(x) > std::numeric_limits<T>::epsilon()) {
        x = one / x;
        p = -p;
      }
    }
    using unsigned_p = std::make_unsigned_t<P>;
    return p == 0 ? one : x * pow(x, static_cast<unsigned_p>(p) - 1);
  } else {
    return static_cast<T>(std::pow(x, static_cast<T>(p)));
  }
}

/// @brief Returns the square of the passed number
/// @param x The number to square
/// @return The square of the input
template <typename T>
constexpr auto square(T x) {
  return x * x;
}

/// @brief Calculates the sum of squares of arguments
/// @param args Variable number of arguments to square and sum
/// @return Sum of squares of all arguments
template <typename... T>
constexpr auto hypotSquare(T... args) {
  return (square(args) + ...);
}

/// @brief Fast hypotenuse calculation for multiple arguments
/// @param args Variable number of arguments
/// @return Square root of sum of squares of arguments
template <typename... T>
constexpr auto fastHypot(T... args) {
  return std::sqrt(hypotSquare(args...));
}

/// @brief Calculates the sum of 1 + 2 + 3+ ... + N using the
///        Gaussian sum formula
/// @param N: Number until which the sum runs
/// @return Sum of integers from 1 to N
template <std::integral T>
constexpr T sumUpToN(const T N) {
  return N * (N + 1) / 2;
}

/// @brief Calculates the factorial of a number
///        n!= n*(n-1)*...*3*2*1
/// @param n: Factor until which the factorial is calculated
/// @return Factorial result
template <std::unsigned_integral T>
constexpr T factorial(const T n) {
  constexpr unsigned bits = std::numeric_limits<T>::digits;

  // max n such that n! fits in T
  constexpr unsigned max_n = bits >= 64   ? 20
                             : bits >= 32 ? 12
                             : bits >= 16 ? 8
                             : bits >= 8  ? 5
                                          : 0;

  if (n > max_n) {
    throw std::overflow_error("factorial overflow");
  }

  T r = 1;
  for (T i = 2; i <= n; i++) {
    r *= i;
  }
  return r;
}

/// @brief Calculate the binomial coefficient
///              n        n!
///                 =  --------
///              k     k!(n-k)!
/// @param n Upper value in binomial coefficient
/// @param k Lower value in binomial coefficient
/// @return Binomial coefficient n choose k
template <std::unsigned_integral T>
constexpr T binomial(const T n, const T k) {
  if (k > n) {
    throw std::invalid_argument("k must be <= n");
  }

  return (factorial<T>(n) / factorial<T>(k)) / factorial<T>(n - k);
}

}  // namespace Acts
