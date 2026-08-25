// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "algebra/impl/detail/eigen_array.hpp"
#include "algebra/impl/detail/eigen_includes.hpp"
#include "detray/algebra/utils/approximately_equal.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

// System include(s).
#include <cassert>
#include <concepts>
#include <limits>
#include <type_traits>

namespace detray::algebra::eigen::math {

/// Transform wrapper class to ensure standard API within different plugins
template <concepts::scalar scalar_t>
struct transform3 {
  /// @name Type definitions for the struct
  /// @{

  /// Scalar type used by the transform
  using scalar_type = scalar_t;

  /// Array type used by the transform
  template <int N>
  using array_type = eigen::array<scalar_type, N>;

  /// Vector and point types
  using vector2 = array_type<2>;
  /// Three-component vector type used for spatial quantities.
  using vector3 = array_type<3>;
  /// Two-dimensional point representation.
  using point2 = vector2;
  /// Three-dimensional point representation.
  using point3 = vector3;

  /// 4x4 matrix type
  using matrix44 =
      typename Eigen::Transform<scalar_type, 3, Eigen::Affine>::MatrixType;

  /// Helper type to cast this to another floating point precision
  template <concepts::scalar o_scalar_t>
  using other_type = transform3<o_scalar_t>;

  /// @}

  /// @name Data objects
  /// @{

  /// Local-to-global affine transformation.
  Eigen::Transform<scalar_type, 3, Eigen::Affine> _data;

  /// @}

  /// Constructor with arguments: t, x, y, z
  ///
  /// @param t the translation (or origin of the new frame)
  /// @param x the x axis of the new frame
  /// @param y the y axis of the new frame
  /// @param z the z axis of the new frame, normal vector for planes
  DETRAY_HOST_DEVICE
  transform3(const vector3 &t, const vector3 &x, const vector3 &y,
             const vector3 &z) {
    _data.setIdentity();

    auto &matrix = _data.matrix();
    matrix.template block<3, 1>(0, 0) = x;
    matrix.template block<3, 1>(0, 1) = y;
    matrix.template block<3, 1>(0, 2) = z;
    matrix.template block<3, 1>(0, 3) = t;
  }

  /// Constructor with arguments: t, z, x
  ///
  /// @param t the translation (or origin of the new frame)
  /// @param z the z axis of the new frame, normal vector for planes
  /// @param x the x axis of the new frame
  DETRAY_HOST_DEVICE
  transform3(const vector3 &t, const vector3 &z, const vector3 &x)
      : transform3(t, x, z.cross(x), z) {}

  /// Constructor with arguments: translation
  ///
  /// @param t is the transform
  DETRAY_HOST_DEVICE
  explicit transform3(const vector3 &t) {
    _data.setIdentity();

    auto &matrix = _data.matrix();
    matrix.template block<3, 1>(0, 3) = t;
  }

  /// Constructor with arguments: matrix
  ///
  /// @param m is the full 4x4 matrix
  DETRAY_HOST_DEVICE
  explicit transform3(const matrix44 &m) { _data.matrix() = m; }

  /// Constructor with arguments: matrix as std::array of scalar
  ///
  /// @param ma is the full 4x4 matrix as a 16 array
  DETRAY_HOST_DEVICE
  explicit transform3(const array_type<16> &ma) {
    _data.matrix() << ma[0], ma[1], ma[2], ma[3], ma[4], ma[5], ma[6], ma[7],
        ma[8], ma[9], ma[10], ma[11], ma[12], ma[13], ma[14], ma[15];
  }

  /// Default constructor: set contents to identity matrices
  DETRAY_HOST_DEVICE
  transform3() { _data.setIdentity(); }

  /// Copy an existing transform.
  /// @param rhs Transform whose matrices are copied.
  transform3(const transform3 &rhs) = default;
  ~transform3() = default;

  /// Equality operator
  /// @param rhs transform to compare with
  /// @return whether the transforms have equal spatial mappings
  DETRAY_HOST_DEVICE
  constexpr bool operator==(const transform3 &rhs) const {
    return (_data.isApprox(rhs._data));
  }

  /// Rotate a vector into / from a frame
  ///
  /// @param m is the rotation matrix
  /// @param v is the vector to be rotated
  /// @return the rotated vector
  template <typename derived_type>
    requires(Eigen::MatrixBase<derived_type>::RowsAtCompileTime == 2 &&
             Eigen::MatrixBase<derived_type>::ColsAtCompileTime == 1)
  DETRAY_HOST_DEVICE static constexpr auto rotate(
      const Eigen::Transform<scalar_type, 3, Eigen::Affine> &m,
      const Eigen::MatrixBase<derived_type> &v) {
    return m.matrix().template block<3, 2>(0, 0) * v;
  }

  /// Rotate a vector into / from a frame
  ///
  /// @param m is the rotation matrix
  /// @param v is the vector to be rotated
  /// @return the rotated vector
  template <typename derived_type>
    requires(Eigen::MatrixBase<derived_type>::RowsAtCompileTime == 3 &&
             Eigen::MatrixBase<derived_type>::ColsAtCompileTime == 1)
  DETRAY_HOST_DEVICE static constexpr auto rotate(
      const Eigen::Transform<scalar_type, 3, Eigen::Affine> &m,
      const Eigen::MatrixBase<derived_type> &v) {
    return m.matrix().template block<3, 3>(0, 0) * v;
  }

  /// This method retrieves the rotation of a transform
  /// @return rotation matrix
  DETRAY_HOST_DEVICE
  constexpr auto rotation() const {
    return _data.matrix().template block<3, 3>(0, 0);
  }

  /// This method retrieves x axis
  /// @return x-axis vector
  DETRAY_HOST_DEVICE
  constexpr point3 x() const {
    return _data.matrix().template block<3, 1>(0, 0);
  }

  /// This method retrieves y axis
  /// @return y-axis vector
  DETRAY_HOST_DEVICE
  constexpr point3 y() const {
    return _data.matrix().template block<3, 1>(0, 1);
  }

  /// This method retrieves z axis
  /// @return z-axis vector
  DETRAY_HOST_DEVICE
  constexpr point3 z() const {
    return _data.matrix().template block<3, 1>(0, 2);
  }

  /// This method retrieves the translation of a transform
  /// @return translation vector
  DETRAY_HOST_DEVICE
  constexpr point3 translation() const {
    return _data.matrix().template block<3, 1>(0, 3);
  }

  /// This method retrieves the 4x4 matrix of a transform
  /// @return forward transformation matrix
  DETRAY_HOST_DEVICE
  constexpr const matrix44 &matrix() const { return _data.matrix(); }

  /// This method retrieves the 4x4 matrix of an inverse transform
  ///
  /// @note the rotation is assumed to be orthonormal, so that the inverse is
  /// the transpose of the rotation and @c -R^T*t as the translation.
  /// @return inverse transformation matrix
  DETRAY_HOST_DEVICE
  constexpr matrix44 matrix_inverse() const {
    matrix44 ret{matrix44::Identity()};

    ret.template block<3, 3>(0, 0) = _data.linear().transpose();
    ret.template block<3, 1>(0, 3) =
        -(_data.linear().transpose() * _data.translation());

    return ret;
  }

  /// This method transform from a point from the local 2D cartesian frame to
  /// the global 3D cartesian frame
  /// @param v point in local plane coordinates
  /// @return point in global coordinates
  template <typename derived_type>
    requires(Eigen::MatrixBase<derived_type>::RowsAtCompileTime == 2 &&
             Eigen::MatrixBase<derived_type>::ColsAtCompileTime == 1)
  DETRAY_HOST_DEVICE constexpr auto point_to_global(
      const Eigen::MatrixBase<derived_type> &v) const {
    return (_data.matrix().template block<3, 2>(0, 0) * v +
            _data.matrix().template block<3, 1>(0, 3));
  }

  /// This method transform from a point from the local 3D cartesian frame to
  /// the global 3D cartesian frame
  /// @param v point in local coordinates
  /// @return point in global coordinates
  template <typename derived_type>
    requires(Eigen::MatrixBase<derived_type>::RowsAtCompileTime == 3 &&
             Eigen::MatrixBase<derived_type>::ColsAtCompileTime == 1)
  DETRAY_HOST_DEVICE constexpr auto point_to_global(
      const Eigen::MatrixBase<derived_type> &v) const {
    return (_data * v);
  }

  /// This method transform from a vector from the global 3D cartesian frame
  /// into the local 3D cartesian frame
  /// @param v vector in global coordinates
  /// @return vector in local coordinates
  template <typename derived_type>
    requires(Eigen::MatrixBase<derived_type>::RowsAtCompileTime == 3 &&
             Eigen::MatrixBase<derived_type>::ColsAtCompileTime == 1)
  DETRAY_HOST_DEVICE constexpr point3 point_to_local(
      const Eigen::MatrixBase<derived_type> &v) const {
    return _data.linear().transpose() * (v - _data.translation());
  }

  /// This method transform from a vector from the local 3D cartesian frame to
  /// the global 3D cartesian frame
  /// @param v vector in local coordinates
  /// @return vector in global coordinates
  template <typename derived_type>
    requires(Eigen::MatrixBase<derived_type>::RowsAtCompileTime == 2 &&
             Eigen::MatrixBase<derived_type>::ColsAtCompileTime == 1)
  DETRAY_HOST_DEVICE constexpr auto vector_to_global(
      const Eigen::MatrixBase<derived_type> &v) const {
    return (_data.linear().matrix().template block<3, 2>(0, 0) * v);
  }

  /// This method transform from a vector from the local 3D cartesian frame to
  /// the global 3D cartesian frame
  /// @param v vector in local coordinates
  /// @return vector in global coordinates
  template <typename derived_type>
    requires(Eigen::MatrixBase<derived_type>::RowsAtCompileTime == 3 &&
             Eigen::MatrixBase<derived_type>::ColsAtCompileTime == 1)
  DETRAY_HOST_DEVICE constexpr auto vector_to_global(
      const Eigen::MatrixBase<derived_type> &v) const {
    return (_data.linear() * v);
  }

  /// This method transform from a vector from the global 3D cartesian frame
  /// into the local 3D cartesian frame
  /// @param v vector in global coordinates
  /// @return vector in local coordinates
  template <typename derived_type>
    requires(Eigen::MatrixBase<derived_type>::RowsAtCompileTime == 3 &&
             Eigen::MatrixBase<derived_type>::ColsAtCompileTime == 1)
  DETRAY_HOST_DEVICE constexpr vector3 vector_to_local(
      const Eigen::MatrixBase<derived_type> &v) const {
    return _data.linear().transpose() * v;
  }
};  // struct transform3

}  // namespace detray::algebra::eigen::math
