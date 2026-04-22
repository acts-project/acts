// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/impl/detail/smatrix_errorcheck.hpp"
#include "detray/algebra/utils/approximately_equal.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

// ROOT/Smatrix include(s).
#include "Math/SMatrix.h"
#include "Math/SVector.h"

// System include(s)
#include <cassert>
#include <concepts>
#include <limits>

namespace detray::algebra::smatrix::math {

/// Transform wrapper class to ensure standard API towards the ROOT::SMatrix lib
template <concepts::scalar scalar_t>
struct transform3 {
  /// @name Type definitions for the struct
  /// @{

  /// Scalar type used by the transform
  using scalar_type = scalar_t;

  /// Array type used by the transform
  template <unsigned int N>
  using array_type = ROOT::Math::SVector<scalar_t, N>;

  /// 3-element "vector" type
  using vector3 = array_type<3>;
  /// Point in 3D space
  using point3 = vector3;
  /// Point in 2D space
  using point2 = array_type<2>;

  /// 4x4 matrix type
  using matrix44 = ROOT::Math::SMatrix<scalar_type, 4, 4>;

  /// Helper type to cast this to another floating point precision
  template <concepts::scalar o_scalar_t>
  using other_type = transform3<o_scalar_t>;

  /// @}

  /// @name Data objects
  /// @{

  matrix44 _data = ROOT::Math::SMatrixIdentity();
  matrix44 _data_inv = ROOT::Math::SMatrixIdentity();

  /// @}

  /// Constructor with arguments: t, x, y, z
  ///
  ///  @param t the translation (or origin of the new frame)
  ///  @param x the x axis of the new frame
  ///  @param y the y axis of the new frame
  ///  @param z the z axis of the new frame, normal vector for planes
  DETRAY_HOST_DEVICE
  transform3(const vector3 &t, const vector3 &x, const vector3 &y,
             const vector3 &z, bool get_inverse = true) {
    _data(0, 0) = x[0];
    _data(1, 0) = x[1];
    _data(2, 0) = x[2];
    _data(0, 1) = y[0];
    _data(1, 1) = y[1];
    _data(2, 1) = y[2];
    _data(0, 2) = z[0];
    _data(1, 2) = z[1];
    _data(2, 2) = z[2];
    _data(0, 3) = t[0];
    _data(1, 3) = t[1];
    _data(2, 3) = t[2];

    if (get_inverse) {
      int ifail = 0;
      _data_inv = _data.Inverse(ifail);
      SMATRIX_CHECK(ifail);
    }
  }

  /// Constructor with arguments: t, z, x
  ///
  /// @param t the translation (or origin of the new frame)
  /// @param z the z axis of the new frame, normal vector for planes
  /// @param x the x axis of the new frame
  DETRAY_HOST
  transform3(const vector3 &t, const vector3 &z, const vector3 &x,
             bool get_inverse = true)
      : transform3(t, x, ROOT::Math::Cross(z, x), z, get_inverse) {}

  /// Constructor with arguments: translation
  ///
  /// @param t is the translation
  DETRAY_HOST
  explicit transform3(const vector3 &t) {
    _data(0, 3) = t[0];
    _data(1, 3) = t[1];
    _data(2, 3) = t[2];

    int ifail = 0;
    _data_inv = _data.Inverse(ifail);
    SMATRIX_CHECK(ifail);
  }

  /// Constructor with arguments: matrix
  ///
  /// @param m is the full 4x4 matrix
  DETRAY_HOST
  explicit transform3(const matrix44 &m) {
    _data = m;

    int ifail = 0;
    _data_inv = _data.Inverse(ifail);
    SMATRIX_CHECK(ifail);
  }

  /// Constructor with arguments: matrix and its inverse
  ///
  /// @param m is the full 4x4 matrix
  /// @param m_inv is the inverse to m
  DETRAY_HOST
  transform3(const matrix44 &m, const matrix44 &m_inv)
      : _data{m}, _data_inv{m_inv} {
    // The assertion will not hold for (casts to) int
    if constexpr (std::floating_point<scalar_type>) {
      [[maybe_unused]] constexpr auto epsilon{
          std::numeric_limits<scalar_type>::epsilon()};
      assert(algebra::approx_equal(matrix44(m * m_inv),
                                   matrix44(ROOT::Math::SMatrixIdentity()),
                                   16.f * epsilon, 1e-6f));
    }
  }

  /// Constructor with arguments: matrix as ROOT::Math::SVector<scalar_t, 16>
  /// of scalars
  ///
  /// @param ma is the full 4x4 matrix as a 16-element array
  DETRAY_HOST
  explicit transform3(const array_type<16> &ma) {
    _data(0, 0) = ma[0];
    _data(1, 0) = ma[4];
    _data(2, 0) = ma[8];
    _data(3, 0) = ma[12];
    _data(0, 1) = ma[1];
    _data(1, 1) = ma[5];
    _data(2, 1) = ma[9];
    _data(3, 1) = ma[13];
    _data(0, 2) = ma[2];
    _data(1, 2) = ma[6];
    _data(2, 2) = ma[10];
    _data(3, 2) = ma[14];
    _data(0, 3) = ma[3];
    _data(1, 3) = ma[7];
    _data(2, 3) = ma[11];
    _data(3, 3) = ma[15];

    int ifail = 0;
    _data_inv = _data.Inverse(ifail);
    // Ignore failures here, since the unit test does manage to trigger an
    // error from ROOT in this place...
  }

  /// Default constructors
  transform3() = default;
  transform3(const transform3 &rhs) = default;
  ~transform3() = default;

  /// Equality operator
  DETRAY_HOST
  constexpr bool operator==(const transform3 &rhs) const {
    return _data == rhs._data;
  }

  /// This method retrieves the rotation of a transform
  DETRAY_HOST
  constexpr auto rotation() const {
    return (_data.template Sub<ROOT::Math::SMatrix<scalar_type, 3, 3> >(0, 0));
  }

  /// This method retrieves x axis
  DETRAY_HOST_DEVICE
  constexpr point3 x() const { return (_data.template SubCol<vector3>(0, 0)); }

  /// This method retrieves y axis
  DETRAY_HOST_DEVICE
  constexpr point3 y() const { return (_data.template SubCol<vector3>(1, 0)); }

  /// This method retrieves z axis
  DETRAY_HOST_DEVICE
  constexpr point3 z() const { return (_data.template SubCol<vector3>(2, 0)); }

  /// This method retrieves the translation of a transform
  DETRAY_HOST
  constexpr vector3 translation() const {
    return (_data.template SubCol<vector3>(3, 0));
  }

  /// This method retrieves the 4x4 matrix of a transform
  DETRAY_HOST
  constexpr matrix44 matrix() const { return _data; }

  /// This method retrieves the 4x4 matrix of an inverse transform
  DETRAY_HOST
  constexpr matrix44 matrix_inverse() const { return _data_inv; }

  /// This method transform from a point from the local 3D cartesian frame to
  /// the global 3D cartesian frame
  DETRAY_HOST
  constexpr point3 point_to_global(const point3 &v) const {
    ROOT::Math::SVector<scalar_type, 4> vector_4;
    vector_4.Place_at(v, 0);
    vector_4[3] = static_cast<scalar_type>(1);
    return ROOT::Math::SVector<scalar_type, 4>(_data * vector_4)
        .template Sub<point3>(0);
  }

  /// This method transform from a vector from the global 3D cartesian frame
  /// into the local 3D cartesian frame
  DETRAY_HOST
  constexpr point3 point_to_local(const point3 &v) const {
    ROOT::Math::SVector<scalar_type, 4> vector_4;
    vector_4.Place_at(v, 0);
    vector_4[3] = static_cast<scalar_type>(1);
    return ROOT::Math::SVector<scalar_type, 4>(_data_inv * vector_4)
        .template Sub<point3>(0);
  }

  /// This method transform from a vector from the local 3D cartesian frame to
  /// the global 3D cartesian frame
  DETRAY_HOST
  constexpr point3 vector_to_global(const vector3 &v) const {
    ROOT::Math::SVector<scalar_type, 4> vector_4;
    vector_4.Place_at(v, 0);
    return ROOT::Math::SVector<scalar_type, 4>(_data * vector_4)
        .template Sub<point3>(0);
  }

  /// This method transform from a vector from the global 3D cartesian frame
  /// into the local 3D cartesian frame
  DETRAY_HOST
  constexpr point3 vector_to_local(const vector3 &v) const {
    ROOT::Math::SVector<scalar_type, 4> vector_4;
    vector_4.Place_at(v, 0);
    return ROOT::Math::SVector<scalar_type, 4>(_data_inv * vector_4)
        .template Sub<point3>(0);
  }
};  // struct transform3

}  // namespace detray::algebra::smatrix::math
