// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

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

  /// Vector and point types
  using vector2 = array_type<2>;
  using vector3 = array_type<3>;
  using point2 = vector2;
  using point3 = vector3;

  /// 4x4 matrix type
  using matrix44 = ROOT::Math::SMatrix<scalar_type, 4, 4>;

  /// Helper type to cast this to another floating point precision
  template <concepts::scalar o_scalar_t>
  using other_type = transform3<o_scalar_t>;

  /// @}

  /// @name Data objects
  /// @{

  matrix44 _data = ROOT::Math::SMatrixIdentity();

  /// @}

  /// Constructor with arguments: t, x, y, z
  ///
  ///  @param t the translation (or origin of the new frame)
  ///  @param x the x axis of the new frame
  ///  @param y the y axis of the new frame
  ///  @param z the z axis of the new frame, normal vector for planes
  DETRAY_HOST_DEVICE
  transform3(const vector3 &t, const vector3 &x, const vector3 &y,
             const vector3 &z) {
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
  }

  /// Constructor with arguments: t, z, x
  ///
  /// @param t the translation (or origin of the new frame)
  /// @param z the z axis of the new frame, normal vector for planes
  /// @param x the x axis of the new frame
  DETRAY_HOST
  transform3(const vector3 &t, const vector3 &z, const vector3 &x)
      : transform3(t, x, ROOT::Math::Cross(z, x), z) {}

  /// Constructor with arguments: translation
  ///
  /// @param t is the translation
  DETRAY_HOST
  explicit transform3(const vector3 &t) {
    _data(0, 3) = t[0];
    _data(1, 3) = t[1];
    _data(2, 3) = t[2];
  }

  /// Constructor with arguments: matrix
  ///
  /// @param m is the full 4x4 matrix
  DETRAY_HOST
  explicit transform3(const matrix44 &m) { _data = m; }

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
  ///
  /// @note the rotation is assumed to be orthonormal, so that the inverse is
  /// the transpose of the rotation and @c -R^T*t as the translation.
  DETRAY_HOST
  constexpr matrix44 matrix_inverse() const {
    matrix44 ret = ROOT::Math::SMatrixIdentity();

    const vector3 t{translation()};

    for (unsigned int i = 0u; i < 3u; ++i) {
      for (unsigned int j = 0u; j < 3u; ++j) {
        ret(i, j) = _data(j, i);
      }
      ret(i, 3) = -ROOT::Math::Dot(_data.template SubCol<vector3>(i, 0), t);
    }

    return ret;
  }

  /// This method transform from a point from the local 3D cartesian frame to
  /// the global 3D cartesian frame
  DETRAY_HOST
  constexpr point3 point_to_global(const point2 &v) const {
    ROOT::Math::SVector<scalar_type, 4> vector_4;
    vector_4.Place_at(v, 0);
    vector_4[3] = static_cast<scalar_type>(1);
    return ROOT::Math::SVector<scalar_type, 4>(_data * vector_4)
        .template Sub<point3>(0);
  }

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
    const vector3 d{v - translation()};

    return point3{ROOT::Math::Dot(x(), d), ROOT::Math::Dot(y(), d),
                  ROOT::Math::Dot(z(), d)};
  }

  /// This method transform from a vector from the local 2D cartesian frame to
  /// the global 3D cartesian frame
  DETRAY_HOST
  constexpr point3 vector_to_global(const vector2 &v) const {
    ROOT::Math::SVector<scalar_type, 4> vector_4;
    vector_4.Place_at(v, 0);
    return ROOT::Math::SVector<scalar_type, 4>(_data * vector_4)
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
    return point3{ROOT::Math::Dot(x(), v), ROOT::Math::Dot(y(), v),
                  ROOT::Math::Dot(z(), v)};
  }
};  // struct transform3

}  // namespace detray::algebra::smatrix::math
