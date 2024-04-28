// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <bit>
#include <cstdint>
#include <limits>

namespace Acts {

/// @brief Fast inverse square root function
/// Taken from https://en.wikipedia.org/wiki/Fast_inverse_square_root
/// @param x the number to calculate the inverse square root of
constexpr float fastInverseSqrt(float x) noexcept {
  // enable only on IEEE 754
  static_assert(std::numeric_limits<double>::is_iec559);

  union {
    float f;
    std::uint32_t i;
  } u = {x};

  u.i = 0x5f3759df - (u.i >> 1);
  u.f *= 1.5f - (0.5f * x * u.f * u.f);
  // 2nd iteration, this can be removed
  // u.f *= 1.5f - (0.5f * x * u.f * u.f);
  return u.f;
}

/// @brief Fast inverse square root function
/// Taken from https://en.wikipedia.org/wiki/Fast_inverse_square_root
/// @param x the number to calculate the inverse square root of
constexpr double fastInverseSqrt(double x) noexcept {
  // enable only on IEEE 754
  static_assert(std::numeric_limits<double>::is_iec559);

  union {
    double f;
    std::uint64_t i;
  } u = {x};

  u.i = 0x5fe6eb50c7b537a9 - (u.i >> 1);
  u.f *= 1.5 - (0.5 * x * u.f * u.f);
  // 2nd iteration, this can be removed
  // u.f *= 1.5 - (0.5 * x * u.f * u.f);
  return u.f;
}

/// @brief Fast power function @see `std::pow`
/// Taken from
/// https://martin.ankerl.com/2012/01/25/optimized-approximative-pow-in-c-and-cpp
/// @param a the base
/// @param b the exponent
constexpr double fastPow(double a, double b) {
  // enable only on IEEE 754
  static_assert(std::numeric_limits<double>::is_iec559);

  union {
    double d;
    std::int32_t x[2];
  } u = {a};
  u.x[1] = static_cast<std::int32_t>(b * (u.x[1] - 1072632447) + 1072632447);
  u.x[0] = 0;
  return u.d;
}

}  // namespace Acts
