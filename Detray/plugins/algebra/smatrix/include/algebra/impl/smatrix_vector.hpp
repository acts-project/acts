// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/definitions/detail/qualifiers.hpp"

// ROOT/Smatrix include(s).
#include <Math/Expression.h>
#include <Math/Functions.h>
#include <Math/SVector.h>
#include <TMath.h>

namespace detray::algebra::smatrix::math {

/// This method retrieves phi from a vector, vector base with rows >= 2
///
/// @param v the input vector
template <concepts::scalar scalar_t, auto N>
  requires(N >= 2)
DETRAY_HOST constexpr scalar_t phi(
    const ROOT::Math::SVector<scalar_t, N> &v) noexcept {
  return static_cast<scalar_t>(TMath::ATan2(v[1], v[0]));
}

template <concepts::scalar scalar_t, class A, auto N>
  requires(N >= 2)
DETRAY_HOST constexpr scalar_t phi(
    const ROOT::Math::VecExpr<A, scalar_t, N> &v) noexcept {
  return static_cast<scalar_t>(TMath::ATan2(v.apply(1), v.apply(0)));
}

/// This method retrieves theta from a vector, vector base with rows >= 3
///
/// @param v the input vector
template <concepts::scalar scalar_t, auto N>
  requires(N >= 3)
DETRAY_HOST constexpr scalar_t theta(
    const ROOT::Math::SVector<scalar_t, N> &v) noexcept {
  return static_cast<scalar_t>(
      TMath::ATan2(TMath::Sqrt(v[0] * v[0] + v[1] * v[1]), v[2]));
}

template <concepts::scalar scalar_t, class A, auto N>
  requires(N >= 3)
DETRAY_HOST constexpr scalar_t theta(
    const ROOT::Math::VecExpr<A, scalar_t, N> &v) noexcept {
  return static_cast<scalar_t>(TMath::ATan2(
      TMath::Sqrt(v.apply(0) * v.apply(0) + v.apply(1) * v.apply(1)),
      v.apply(2)));
}

/// This method retrieves the norm of a vector, no dimension restriction
///
/// @param v the input vector
template <concepts::scalar scalar_t, auto N>
DETRAY_HOST constexpr scalar_t norm(const ROOT::Math::SVector<scalar_t, N> &v) {
  return static_cast<scalar_t>(TMath::Sqrt(ROOT::Math::Dot(v, v)));
}

template <concepts::scalar scalar_t, class A, auto N>
DETRAY_HOST constexpr scalar_t norm(
    const ROOT::Math::VecExpr<A, scalar_t, N> &v) {
  return static_cast<scalar_t>(TMath::Sqrt(ROOT::Math::Dot(v, v)));
}

/// This method retrieves the pseudo-rapidity from a vector or vector base with
/// rows >= 3
///
/// @param v the input vector
template <concepts::scalar scalar_t, auto N>
  requires(N >= 3)
DETRAY_HOST constexpr scalar_t eta(
    const ROOT::Math::SVector<scalar_t, N> &v) noexcept {
  return static_cast<scalar_t>(TMath::ATanH(v[2] / norm(v)));
}

template <concepts::scalar scalar_t, class A, auto N>
  requires(N >= 3)
DETRAY_HOST constexpr scalar_t eta(
    const ROOT::Math::VecExpr<A, scalar_t, N> &v) noexcept {
  return static_cast<scalar_t>(TMath::ATanH(v.apply(2) / norm(v)));
}

/// This method retrieves the perpendicular magnitude of a vector with rows >= 2
///
/// @param v the input vector
template <concepts::scalar scalar_t, auto N>
  requires(N >= 2)
DETRAY_HOST constexpr scalar_t perp(
    const ROOT::Math::SVector<scalar_t, N> &v) noexcept {
  return static_cast<scalar_t>(TMath::Sqrt(v[0] * v[0] + v[1] * v[1]));
}

template <concepts::scalar scalar_t, class A, auto N>
  requires(N >= 2)
DETRAY_HOST constexpr scalar_t perp(
    const ROOT::Math::VecExpr<A, scalar_t, N> &v) noexcept {
  return static_cast<scalar_t>(
      TMath::Sqrt(v.apply(0) * v.apply(0) + v.apply(1) * v.apply(1)));
}

/// Get a normalized version of the input vector
///
/// @param v the input vector
template <concepts::scalar scalar_t, auto N>
DETRAY_HOST constexpr ROOT::Math::SVector<scalar_t, N> normalize(
    const ROOT::Math::SVector<scalar_t, N> &v) {
  return ROOT::Math::Unit(v);
}

/// Get a normalized version of the input vector
///
/// @param v the input vector
template <concepts::scalar scalar_t, class A, auto N>
DETRAY_HOST constexpr ROOT::Math::SVector<scalar_t, N> normalize(
    const ROOT::Math::VecExpr<A, scalar_t, N> &v) {
  return ROOT::Math::Unit(v);
}

/// Dot product between two input vectors
///
/// @param a the first input vector
/// @param b the second input vector
///
/// @return the scalar dot product value
template <concepts::scalar scalar_t, auto N>
DETRAY_HOST constexpr scalar_t dot(const ROOT::Math::SVector<scalar_t, N> &a,
                                   const ROOT::Math::SVector<scalar_t, N> &b) {
  return ROOT::Math::Dot(a, b);
}

/// Dot product between two input vectors
///
/// @param a the first input vector
/// @param b the second input vector
///
/// @return the scalar dot product value
template <concepts::scalar scalar_t, class A, auto N>
DETRAY_HOST constexpr scalar_t dot(
    const ROOT::Math::SVector<scalar_t, N> &a,
    const ROOT::Math::VecExpr<A, scalar_t, N> &b) {
  return ROOT::Math::Dot(a, b);
}

/// Dot product between two input vectors
///
/// @param a the first input vector
/// @param b the second input vector
///
/// @return the scalar dot product value
template <concepts::scalar scalar_t, class A, auto N>
DETRAY_HOST constexpr scalar_t dot(const ROOT::Math::VecExpr<A, scalar_t, N> &a,
                                   const ROOT::Math::SVector<scalar_t, N> &b) {
  return ROOT::Math::Dot(a, b);
}

/// Dot product between two input vectors
///
/// @param a the first input vector
/// @param b the second input vector
///
/// @return the scalar dot product value
template <concepts::scalar scalar_t, class A, class B, auto N>
DETRAY_HOST constexpr scalar_t dot(
    const ROOT::Math::VecExpr<A, scalar_t, N> &a,
    const ROOT::Math::VecExpr<B, scalar_t, N> &b) {
  return ROOT::Math::Dot(a, b);
}

/// Dot product between two input vectors
///
/// @param a the first input vector
/// @param b the second input vector
///
/// @return the scalar dot product value
template <concepts::scalar scalar_t, class A, auto N>
DETRAY_HOST constexpr scalar_t dot(
    const ROOT::Math::SMatrix<scalar_t, N, 1> &a,
    const ROOT::Math::VecExpr<A, scalar_t, N> &b) {
  return ROOT::Math::Dot(a.Col(0), b);
}

/// Dot product between two input vectors
///
/// @param a the first input vector
/// @param b the second input vector
///
/// @return the scalar dot product value
template <concepts::scalar scalar_t, class A, auto N>
DETRAY_HOST constexpr scalar_t dot(
    const ROOT::Math::VecExpr<A, scalar_t, N> &a,
    const ROOT::Math::SMatrix<scalar_t, N, 1> &b) {
  return dot(b, a);
}

/// Dot product between two input vectors
///
/// @param a the first input vector
/// @param b the second input vector
///
/// @return the scalar dot product value
template <concepts::scalar scalar_t, auto N>
DETRAY_HOST constexpr scalar_t dot(const ROOT::Math::SMatrix<scalar_t, N, 1> &a,
                                   const ROOT::Math::SVector<scalar_t, N> &b) {
  return ROOT::Math::Dot(a.Col(0), b);
}

/// Dot product between two input vectors
///
/// @param a the first input vector
/// @param b the second input vector
///
/// @return the scalar dot product value
template <concepts::scalar scalar_t, auto N>
DETRAY_HOST constexpr scalar_t dot(
    const ROOT::Math::SVector<scalar_t, N> &a,
    const ROOT::Math::SMatrix<scalar_t, N, 1> &b) {
  return dot(b, a);
}

/// Cross product between two input vectors
///
/// @param a the first input vector
/// @param b the second input vector
///
/// @return a vector (expression) representing the cross product
template <concepts::scalar scalar_t>
DETRAY_HOST constexpr ROOT::Math::SVector<scalar_t, 3> cross(
    const ROOT::Math::SVector<scalar_t, 3> &a,
    const ROOT::Math::SVector<scalar_t, 3> &b) {
  return ROOT::Math::Cross(a, b);
}

/// Cross product between two input vectors
///
/// @param a the first input vector
/// @param b the second input vector
///
/// @return a vector (expression) representing the cross product
template <concepts::scalar scalar_t, class A>
DETRAY_HOST constexpr ROOT::Math::SVector<scalar_t, 3> cross(
    const ROOT::Math::SVector<scalar_t, 3> &a,
    const ROOT::Math::VecExpr<A, scalar_t, 3> &b) {
  return ROOT::Math::Cross(a, b);
}

/// Cross product between two input vectors
///
/// @param a the first input vector
/// @param b the second input vector
///
/// @return a vector (expression) representing the cross product
template <concepts::scalar scalar_t, class A>
DETRAY_HOST constexpr ROOT::Math::SVector<scalar_t, 3> cross(
    const ROOT::Math::VecExpr<A, scalar_t, 3> &a,
    const ROOT::Math::SVector<scalar_t, 3> &b) {
  return ROOT::Math::Cross(a, b);
}

/// Cross product between two input vectors
///
/// @param a the first input vector
/// @param b the second input vector
///
/// @return a vector (expression) representing the cross product
template <concepts::scalar scalar_t, class A, class B>
DETRAY_HOST constexpr ROOT::Math::SVector<scalar_t, 3> cross(
    const ROOT::Math::VecExpr<A, scalar_t, 3> &a,
    const ROOT::Math::VecExpr<B, scalar_t, 3> &b) {
  return ROOT::Math::Cross(a, b);
}

/// Cross product between vector3 and matrix<3,1>
///
/// @param a the first input vector
/// @param b the second input matrix<3,1>
///
/// @return a vector (expression) representing the cross product
template <concepts::scalar scalar_t>
DETRAY_HOST constexpr ROOT::Math::SVector<scalar_t, 3> cross(
    const ROOT::Math::SVector<scalar_t, 3> &a,
    const ROOT::Math::SMatrix<scalar_t, 3, 1> &b) {
  return ROOT::Math::Cross(a, b.Col(0));
}

/// Cross product between matrix<3,1> and vector3
///
/// @param a the second input matrix<3,1>
/// @param b the first input vector
///
/// @return a vector (expression) representing the cross product
template <concepts::scalar scalar_t>
DETRAY_HOST constexpr ROOT::Math::SVector<scalar_t, 3> cross(
    const ROOT::Math::SMatrix<scalar_t, 3, 1> &a,
    const ROOT::Math::SVector<scalar_t, 3> &b) {
  return ROOT::Math::Cross(a.Col(0), b);
}

/// Cross product between two matrix<3,1>
///
/// @param a the second input matrix<3,1>
/// @param b the first input matrix<3,1>
///
/// @return a vector (expression) representing the cross product
template <concepts::scalar scalar_t>
DETRAY_HOST constexpr ROOT::Math::SVector<scalar_t, 3> cross(
    const ROOT::Math::SMatrix<scalar_t, 3, 1> &a,
    const ROOT::Math::SMatrix<scalar_t, 3, 1> &b) {
  return ROOT::Math::Cross(a.Col(0), b.Col(0));
}

}  // namespace detray::algebra::smatrix::math
