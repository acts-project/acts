// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "algebra/impl/array_soa_simd.hpp"
#include "detray/algebra/concepts.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

// System include(s)
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <type_traits>

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

/// Lane-wise overloads of common math functions
/// @{

/// Unary function applied to every lane
#define DETRAY_ARRAY_SOA_UNARY_MATH(FN)                  \
  template <concepts::value T, std::size_t W>            \
  DETRAY_HOST_DEVICE constexpr array_soa::simd<T, W> FN( \
      const array_soa::simd<T, W> &v) {                  \
    array_soa::simd<T, W> ret;                           \
    DETRAY_UNROLL_N(W)                                   \
    for (std::size_t i = 0u; i < W; ++i) {               \
      ret[i] = static_cast<T>(FN(v[i]));                 \
    }                                                    \
    return ret;                                          \
  }

/// Binary function applied to every lane, with a lane bundle or a single
/// value as second argument
#define DETRAY_ARRAY_SOA_BINARY_MATH(FN)                                \
  template <concepts::value T, std::size_t W>                           \
  DETRAY_HOST_DEVICE constexpr array_soa::simd<T, W> FN(                \
      const array_soa::simd<T, W> &a, const array_soa::simd<T, W> &b) { \
    array_soa::simd<T, W> ret;                                          \
    DETRAY_UNROLL_N(W)                                                  \
    for (std::size_t i = 0u; i < W; ++i) {                              \
      ret[i] = static_cast<T>(FN(a[i], b[i]));                          \
    }                                                                   \
    return ret;                                                         \
  }                                                                     \
  template <concepts::value T, std::size_t W>                           \
  DETRAY_HOST_DEVICE constexpr array_soa::simd<T, W> FN(                \
      const array_soa::simd<T, W> &a, std::type_identity_t<T> b) {      \
    array_soa::simd<T, W> ret;                                          \
    DETRAY_UNROLL_N(W)                                                  \
    for (std::size_t i = 0u; i < W; ++i) {                              \
      ret[i] = static_cast<T>(FN(a[i], b));                             \
    }                                                                   \
    return ret;                                                         \
  }                                                                     \
  template <concepts::value T, std::size_t W>                           \
  DETRAY_HOST_DEVICE constexpr array_soa::simd<T, W> FN(                \
      std::type_identity_t<T> a, const array_soa::simd<T, W> &b) {      \
    array_soa::simd<T, W> ret;                                          \
    DETRAY_UNROLL_N(W)                                                  \
    for (std::size_t i = 0u; i < W; ++i) {                              \
      ret[i] = static_cast<T>(FN(a, b[i]));                             \
    }                                                                   \
    return ret;                                                         \
  }

// clang-format off
DETRAY_ARRAY_SOA_UNARY_MATH(abs)
DETRAY_ARRAY_SOA_UNARY_MATH(acos)
DETRAY_ARRAY_SOA_UNARY_MATH(asin)
DETRAY_ARRAY_SOA_UNARY_MATH(atan)
DETRAY_ARRAY_SOA_UNARY_MATH(atanh)
DETRAY_ARRAY_SOA_UNARY_MATH(ceil)
DETRAY_ARRAY_SOA_UNARY_MATH(cos)
DETRAY_ARRAY_SOA_UNARY_MATH(cosh)
DETRAY_ARRAY_SOA_UNARY_MATH(exp)
DETRAY_ARRAY_SOA_UNARY_MATH(fabs)
DETRAY_ARRAY_SOA_UNARY_MATH(floor)
DETRAY_ARRAY_SOA_UNARY_MATH(log)
DETRAY_ARRAY_SOA_UNARY_MATH(log10)
DETRAY_ARRAY_SOA_UNARY_MATH(sin)
DETRAY_ARRAY_SOA_UNARY_MATH(sinh)
DETRAY_ARRAY_SOA_UNARY_MATH(sqrt)
DETRAY_ARRAY_SOA_UNARY_MATH(tan)
DETRAY_ARRAY_SOA_UNARY_MATH(tanh)

DETRAY_ARRAY_SOA_BINARY_MATH(atan2)
DETRAY_ARRAY_SOA_BINARY_MATH(copysign)
DETRAY_ARRAY_SOA_BINARY_MATH(hypot)
DETRAY_ARRAY_SOA_BINARY_MATH(max)
DETRAY_ARRAY_SOA_BINARY_MATH(min)
DETRAY_ARRAY_SOA_BINARY_MATH(pow)
// clang-format on

#undef DETRAY_ARRAY_SOA_UNARY_MATH
#undef DETRAY_ARRAY_SOA_BINARY_MATH

template <concepts::value T, std::size_t W>
DETRAY_HOST_DEVICE constexpr array_soa::mask<W> signbit(
    const array_soa::simd<T, W> &v) {
  array_soa::mask<W> ret;
  DETRAY_UNROLL_N(W)
  for (std::size_t i = 0u; i < W; ++i) {
    ret[i] = signbit(v[i]);
  }
  return ret;
}

template <concepts::value T, std::size_t W>
DETRAY_HOST_DEVICE constexpr array_soa::simd<T, W> fma(
    const array_soa::simd<T, W> &x, const array_soa::simd<T, W> &y,
    const array_soa::simd<T, W> &z) {
  array_soa::simd<T, W> ret;
  DETRAY_UNROLL_N(W)
  for (std::size_t i = 0u; i < W; ++i) {
    ret[i] = static_cast<T>(fma(x[i], y[i], z[i]));
  }
  return ret;
}
/// @}

}  // namespace detray::algebra::math
