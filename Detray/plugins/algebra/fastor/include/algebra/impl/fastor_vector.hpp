// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/algebra/concepts.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

// Fastor include(s).
#ifdef _MSC_VER
#pragma warning(disable : 4244 4701 4702)
#endif  // MSVC
#include <Fastor/Fastor.h>
#ifdef _MSC_VER
#pragma warning(default : 4244 4701 4702)
#endif  // MSVC

namespace detray::algebra::fastor::math {

// Note that for all `Fastor::AbstractTensors`, the `N` template parameters
// refers to the number of dimensions, not the number of elements.

/// This method retrieves phi from a vector @param v
template <typename Derived, auto N>
DETRAY_HOST_DEVICE constexpr auto phi(
    const Fastor::AbstractTensor<Derived, N> &a) {
  // we first force evaluation of whatever was passed in.
  auto v = Fastor::evaluate(a);
  // using `auto` relieves us from having to extract the exact dimension of
  // the vector from somewhere. For all intents and purposes, we can consider
  // the type to be `Fastor::Tensor<scalar_t, SIZE>`.

  // we use the cmath version of `atan2` because Fastor's `atan2` works on
  // `Fastor::Tensor`s element-wise, which we don't want.
  return algebra::math::atan2(v[1], v[0]);
}

/// This method retrieves theta from a vector, vector base with rows >= 3
///
/// @param v the input vector
template <typename Derived, auto N>
DETRAY_HOST constexpr auto theta(
    const Fastor::AbstractTensor<Derived, N> &a) noexcept {
  // we first force evaluation of whatever was passed in.
  auto v = Fastor::evaluate(a);
  // using `auto` relieves us from having to extract the exact dimension of
  // the vector from somewhere. For all intents and purposes, we can consider
  // the type to be `Fastor::Tensor<scalar_t, SIZE>`.

  // we use the cmath version of `atan2` because Fastor's `atan2` works on
  // `Fastor::Tensor`s element-wise, which we don't want.
  return algebra::math::atan2(Fastor::norm(v(Fastor::fseq<0, 2>())), v[2]);
}

/// This method retrieves the perpendicular magnitude of a vector with rows >= 2
///
/// @param v the input vector
template <typename Derived, auto N>
DETRAY_HOST constexpr auto perp(
    const Fastor::AbstractTensor<Derived, N> &a) noexcept {
  // we first force evaluation of whatever was passed in.
  // using `auto` relieves us from having to extract the exact dimension of
  // the vector from somewhere. For all intents and purposes, we can consider
  // the type to be `Fastor::Tensor<scalar_t, SIZE>`.
  auto v = Fastor::evaluate(a);

  // we use the cmath version of `sqrt` because Fastor's `sqrt` works on
  // `Fastor::Tensor`s element-wise, which we don't want.
  return algebra::math::sqrt(
      Fastor::inner(v(Fastor::fseq<0, 2>()), v(Fastor::fseq<0, 2>())));
}

/// This method retrieves the norm of a vector, no dimension restriction
///
/// @param v the input vector
template <typename Derived, auto N>
DETRAY_HOST constexpr auto norm(const Fastor::AbstractTensor<Derived, N> &v) {
  return Fastor::norm(v);
}

/// This method retrieves the pseudo-rapidity from a vector or vector base with
/// rows >= 3
///
/// @param v the input vector
template <typename Derived, auto N>
DETRAY_HOST constexpr auto eta(
    const Fastor::AbstractTensor<Derived, N> &a) noexcept {
  auto v = Fastor::evaluate(a);

  // we use the cmath version of `atanh` because Fastor's `atanh` works on
  // `Fastor::Tensor`s element-wise, which we don't want.
  return algebra::math::atanh(v[2] / Fastor::norm(v));
}

/// Get a normalized version of the input vector
///
/// @param v the input vector
template <typename Derived, auto N>
DETRAY_HOST constexpr auto normalize(
    const Fastor::AbstractTensor<Derived, N> &v) {
  return v / Fastor::norm(v);
}

/// Dot product between two input vectors
///
/// @param a the first input vector
/// @param b the second input vector
///
/// @return the scalar dot product value
template <typename Derived0, auto N0, typename Derived1, auto N1>
DETRAY_HOST constexpr auto dot(const Fastor::AbstractTensor<Derived0, N0> &a,
                               const Fastor::AbstractTensor<Derived1, N1> &b) {
  return Fastor::inner(a, b);
}

/// Dot product between two pure vectors (of type `Tensor<scalar_t, N>`)
///
/// @param a the first input vector
/// @param b the second input vector
///
/// @return the scalar dot product value
template <concepts::scalar scalar_t, auto N>
DETRAY_HOST_DEVICE constexpr scalar_t dot(
    const Fastor::Tensor<scalar_t, N> &a,
    const Fastor::Tensor<scalar_t, N> &b) {
  return Fastor::inner(a, b);
}

/// Dot product between a vector and a matrix slice
///
/// @param a the first input: a vector (`Tensor<scalar_t, N>`)
/// @param b the second input: a matrix (`Tensor<scalar_t, N, 1>`)
///
/// @return the scalar dot product value
template <concepts::scalar scalar_t, auto N>
DETRAY_HOST constexpr scalar_t dot(const Fastor::Tensor<scalar_t, N> &a,
                                   const Fastor::Tensor<scalar_t, N, 1> &b) {
  // We need to specify the type of the Tensor slice because Fastor by default
  // is lazy, so it returns an intermediate type which does not play well with
  // the Fastor::inner function.
  return Fastor::inner(a,
                       Fastor::Tensor<scalar_t, N>(b(Fastor::fseq<0, N>(), 0)));
}

/// Dot product between a matrix slice and a vector
///
/// @param a the first input: a matrix (`Tensor<scalar_t, N, 1>`)
/// @param b the second input: a vector (`Tensor<scalar_t, N>`)
///
/// @return the scalar dot product value
template <concepts::scalar scalar_t, auto N>
DETRAY_HOST constexpr scalar_t dot(const Fastor::Tensor<scalar_t, N, 1> &a,
                                   const Fastor::Tensor<scalar_t, N> &b) {
  return Fastor::inner(Fastor::Tensor<scalar_t, N>(a(Fastor::fseq<0, N>(), 0)),
                       b);
}

/// Dot product between two matrix slices
///
/// @param a the first input: a matrix (`Tensor<scalar_t, N, 1>`)
/// @param b the second input: a matrix (`Tensor<scalar_t, N, 1>`)
///
/// @return the scalar dot product value
template <concepts::scalar scalar_t, auto N>
DETRAY_HOST constexpr scalar_t dot(const Fastor::Tensor<scalar_t, N, 1> &a,
                                   const Fastor::Tensor<scalar_t, N, 1> &b) {
  return Fastor::inner(Fastor::Tensor<scalar_t, N>(a(Fastor::fseq<0, 3>(), 0)),
                       Fastor::Tensor<scalar_t, N>(b(Fastor::fseq<0, 3>(), 0)));
}

/// Cross product between two input vectors
///
/// @param a the first input vector
/// @param b the second input vector
///
/// @return a vector (expression) representing the cross product
template <typename Derived0, auto N0, typename Derived1, auto N1>
DETRAY_HOST constexpr auto cross(
    const Fastor::AbstractTensor<Derived0, N0> &a,
    const Fastor::AbstractTensor<Derived1, N1> &b) {
  return Fastor::cross(a, b);
}

/// Cross product between two pure vectors (of type `Tensor<scalar_t, 3>`)
///
/// @param a the first input vector
/// @param b the second input vector
///
/// @return a vector (expression) representing the cross product
template <concepts::scalar scalar_t>
DETRAY_HOST constexpr auto cross(const Fastor::Tensor<scalar_t, 3> &a,
                                 const Fastor::Tensor<scalar_t, 3> &b) {
  return Fastor::cross(a, b);
}

/// Cross product between a vector and a matrix slice
///
/// @param a the first input: a vector (`Tensor<scalar_t, 3>`)
/// @param b the second input: a matrix (`Tensor<scalar_t, 3, 1>`)
///
/// @return a vector (expression) representing the cross product
template <concepts::scalar scalar_t>
DETRAY_HOST constexpr auto cross(const Fastor::Tensor<scalar_t, 3> &a,
                                 const Fastor::Tensor<scalar_t, 3, 1> &b) {
  // We need to specify the type of the Tensor slice because Fastor by default
  // is lazy, so it returns an intermediate type which does not play well with
  // the Fastor::cross function.
  return Fastor::cross(a,
                       Fastor::Tensor<scalar_t, 3>(b(Fastor::fseq<0, 3>(), 0)));
}

/// Cross product between a matrix slice and a vector
///
/// @param a the first input: a matrix (`Tensor<scalar_t, 3, 1>`)
/// @param b the second input: a vector (`Tensor<scalar_t, 3>`)
///
/// @return a vector (expression) representing the cross product
template <concepts::scalar scalar_t>
DETRAY_HOST constexpr auto cross(const Fastor::Tensor<scalar_t, 3, 1> &a,
                                 const Fastor::Tensor<scalar_t, 3> &b) {
  return Fastor::cross(Fastor::Tensor<scalar_t, 3>(a(Fastor::fseq<0, 3>(), 0)),
                       b);
}

/// Cross product between two matrix slices
///
/// @param a the second input matrix (`Tensor<scalar_t, 3, 1>`)
/// @param b the first input matrix (`Tensor<scalar_t, 3, 1>`)
///
/// @return a vector (expression) representing the cross product
template <concepts::scalar scalar_t>
DETRAY_HOST constexpr auto cross(const Fastor::Tensor<scalar_t, 3, 1> &a,
                                 const Fastor::Tensor<scalar_t, 3, 1> &b) {
  return Fastor::cross(Fastor::Tensor<scalar_t, 3>(a(Fastor::fseq<0, 3>(), 0)),
                       Fastor::Tensor<scalar_t, 3>(b(Fastor::fseq<0, 3>(), 0)));
}

}  // namespace detray::algebra::fastor::math
