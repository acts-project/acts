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

template <typename T>
constexpr auto square(T x) {
  return x * x;
}

template <typename... T>
constexpr auto hypotSquare(T... args) {
  return (square(args) + ...);
}

template <typename... T>
constexpr auto fastHypot(T... args) {
  return std::sqrt(hypotSquare(args...));
}

/// @brief Helper struct to calculate sin/cos in one go
template <typename T>
struct sincos {
  /// @brief Standard constructor taking an angle
  /// @param alpha: Angle to consider
  sincos(const T alpha) : cs{std::cos(alpha)}, sn{std::sin(alpha)} {}
  T cs{0.};
  T sn{0.};
};


/// @brief Calculates the sum of 1 + 2 + 3+ ... + N using the
///        Gaussian sum formula
/// @param N: Upper limit of the sum
template<typename T,
        typename std::enable_if_t<std::is_integral_v<T>, bool> = true>
constexpr T sumUp (const T N) {
    return N *(N +1) / 2;
} 



}  // namespace Acts
