// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>

namespace Acts {

/// @brief Returns the absolute of a number
///        (Can be removed for c++ 23)
/// @param n The number to take absolute value of
/// @return The absolute value of the input
template <typename T>
constexpr T abs(const T n) {
  if constexpr (std::is_signed_v<T>) {
    if (n < 0) {
      return -n;
    }
  }
  return n;
}

/// @brief Calculates the ordinary power of the number x.
/// @param x: Number to take the power from
/// @param p: Power to take
/// @return x raised to the power p
template <typename T, std::integral P>
constexpr T pow(T x, P p) {
  constexpr T one = 1;
  if constexpr (std::is_signed_v<P>) {
    if (p < 0 && abs(x) > std::numeric_limits<T>::epsilon()) {
      x = one / x;
      p = -p;
    }
  }
  using unsigned_p = std::make_unsigned_t<P>;
  return p == 0 ? one : x * pow(x, static_cast<unsigned_p>(p) - 1);
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
///        N!= N*(N-1)....*3*2*1
/// @param upperN: Upper factor until which the factorial is calculated
/// @param lowerN: Optional argument to remove the first factors from the calculation
/// @return Factorial result
template <std::integral T>
constexpr T factorial(const T upperN, const T lowerN = 1) {
  constexpr T one = 1;
  const T& limit = std::max(one, lowerN);
  return upperN >= limit ? upperN * factorial(upperN - 1, limit) : one;
}
/// @brief Calculate the binomial coefficient
///              n        n!
///                 =  --------
///              k     k!(n-k)!
/// @param n Upper value in binomial coefficient
/// @param k Lower value in binomial coefficient
/// @return Binomial coefficient n choose k
template <std::integral T>
constexpr T binomial(const T n, const T k) {
  return factorial<T>(n, n - k + 1) / factorial<T>(k);
}

}  // namespace Acts
