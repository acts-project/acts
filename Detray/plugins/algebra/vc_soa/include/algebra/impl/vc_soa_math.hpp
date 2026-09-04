// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "algebra/impl/vc_soa_concepts.hpp"

// Vc include(s).
#ifdef _MSC_VER
#pragma warning(push, 0)
#endif  // MSVC
#include <Vc/Vc>
#ifdef _MSC_VER
#pragma warning(pop)
#endif  // MSVC

// System include(s)
#include <algorithm>
#include <cmath>

namespace detray::algebra::math {

/// Math functions on single values
/// @{
using std::abs;
using std::acos;
using std::asin;
using std::atan;
using std::atan2;
using std::atanh;
using std::ceil;
using std::copysign;
using std::cos;
using std::cosh;
using std::exp;
using std::fabs;
using std::floor;
using std::fma;
using std::hypot;
using std::log;
using std::log10;
using std::max;
using std::min;
using std::pow;
using std::signbit;
using std::sin;
using std::sinh;
using std::sqrt;
using std::tan;
using std::tanh;
/// @}

/// Vc overloads of common math functions
/// @{
template <detray::concepts::vc_simd_vector T>
constexpr decltype(auto) abs(T &&vec) {
  return Vc::abs(std::forward<T>(vec));
}

template <detray::concepts::vc_simd_vector T>
constexpr decltype(auto) fabs(T &&vec) {
  return Vc::abs(std::forward<T>(vec));
}

template <detray::concepts::vc_simd_vector T>
constexpr decltype(auto) sqrt(T &&vec) {
  return Vc::sqrt(std::forward<T>(vec));
}

template <detray::concepts::vc_simd_vector T>
constexpr decltype(auto) exp(T &&vec) {
  return Vc::exp(std::forward<T>(vec));
}

template <detray::concepts::vc_simd_vector T>
constexpr decltype(auto) log(T &&vec) {
  return Vc::log(std::forward<T>(vec));
}

template <detray::concepts::vc_simd_vector T>
constexpr decltype(auto) sin(T &&vec) {
  return Vc::sin(std::forward<T>(vec));
}

template <detray::concepts::vc_simd_vector T>
constexpr decltype(auto) asin(T &&vec) {
  return Vc::asin(std::forward<T>(vec));
}

template <detray::concepts::vc_simd_vector T>
constexpr decltype(auto) cos(T &&vec) {
  return Vc::cos(std::forward<T>(vec));
}

template <detray::concepts::vc_simd_vector T>
constexpr decltype(auto) tan(T &&vec) {
  // It seems there is no dedicated @c Vc::tan function ?
  return Vc::sin(std::forward<T>(vec)) / Vc::cos(std::forward<T>(vec));
}

template <detray::concepts::vc_simd_vector T>
constexpr decltype(auto) atan(T &&vec) {
  return Vc::atan(std::forward<T>(vec));
}

template <detray::concepts::vc_simd_vector T,
          detray::concepts::vc_simd_vector S>
constexpr decltype(auto) copysign(T &&mag, S &&sgn) {
  return Vc::copysign(std::forward<T>(mag), std::forward<S>(sgn));
}

template <detray::concepts::vc_simd_vector T>
constexpr decltype(auto) min(T &&vec) {
  return Vc::min(std::forward<T>(vec));
}

template <detray::concepts::vc_simd_vector T>
constexpr decltype(auto) max(T &&vec) {
  return Vc::max(std::forward<T>(vec));
}

template <detray::concepts::vc_simd_vector T>
constexpr decltype(auto) signbit(T &&vec) {
  return Vc::isnegative(std::forward<T>(vec));
}

template <detray::concepts::vc_simd_vector T>
constexpr decltype(auto) fma(T &&x, T &&y, T &&z) {
  return Vc::fma(std::forward<T>(x), std::forward<T>(y), std::forward<T>(z));
}
/// @}

}  // namespace detray::algebra::math
