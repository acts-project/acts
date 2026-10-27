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
#include "detray/algebra/utils/approximately_equal.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

// Eigen include(s).
#ifdef _MSC_VER
#pragma warning(push, 0)
#endif  // MSVC
#ifdef __NVCC_DIAG_PRAGMA_SUPPORT__
#pragma nv_diagnostic push
#pragma nv_diag_suppress 20011
#pragma nv_diag_suppress 20012
#pragma nv_diag_suppress 20014
#endif  // __NVCC_DIAG_PRAGMA_SUPPORT__
#include <Eigen/Core>
#include <Eigen/Geometry>
#ifdef _MSC_VER
#pragma warning(pop)
#endif  // MSVC
#ifdef __NVCC_DIAG_PRAGMA_SUPPORT__
#pragma nv_diagnostic pop
#endif  // __NVCC_DIAG_PRAGMA_SUPPORT__

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
  /// Global-to-local affine transformation.
  Eigen::Transform<scalar_type, 3, Eigen::Affine> _data_inv;

  /// @}

  /// Constructor with arguments: t, x, y, z
  ///
  /// @param t the translation (or origin of the new frame)
  /// @param x the x axis of the new frame
  /// @param y the y axis of the new frame
  /// @param z the z axis of the new frame, normal vector for planes
  /// @param get_inverse whether to compute and store the inverse
  DETRAY_HOST_DEVICE
  transform3(const vector3 &t, const vector3 &x, const vector3 &y,
             const vector3 &z, bool get_inverse = true) {
    _data.setIdentity();

    auto &matrix = _data.matrix();
    matrix.template block<3, 1>(0, 0) = x;
    matrix.template block<3, 1>(0, 1) = y;
    matrix.template block<3, 1>(0, 2) = z;
    matrix.template block<3, 1>(0, 3) = t;

    if (get_inverse) {
      _data_inv = _data.inverse();
    } else {
      _data_inv.setIdentity();
    }
  }

  /// Constructor with arguments: t, z, x
  ///
  /// @param t the translation (or origin of the new frame)
  /// @param z the z axis of the new frame, normal vector for planes
  /// @param x the x axis of the new frame
  /// @param get_inverse whether to compute and store the inverse
  DETRAY_HOST_DEVICE
  transform3(const vector3 &t, const vector3 &z, const vector3 &x,
             bool get_inverse = true)
      : transform3(t, x, z.cross(x), z, get_inverse) {}

  /// Constructor with arguments: translation
  ///
  /// @param t is the transform
  DETRAY_HOST_DEVICE
  explicit transform3(const vector3 &t) {
    _data.setIdentity();

    auto &matrix = _data.matrix();
    matrix.template block<3, 1>(0, 3) = t;

    _data_inv = _data.inverse();
  }

  /// Constructor with arguments: matrix
  ///
  /// @param m is the full 4x4 matrix
  DETRAY_HOST_DEVICE
  explicit transform3(const matrix44 &m) {
    _data.matrix() = m;

    _data_inv = _data.inverse();
  }

  /// Constructor with arguments: matrix and its inverse
  ///
  /// @param m is the full 4x4 matrix
  /// @param m_inv is the inverse to m
  DETRAY_HOST_DEVICE
  transform3(const matrix44 &m, const matrix44 &m_inv)
      : _data{m}, _data_inv{m_inv} {
    // The assertion will not hold for (casts to) int
    if constexpr (std::floating_point<scalar_type>) {
      [[maybe_unused]] constexpr auto epsilon{
          std::numeric_limits<scalar_type>::epsilon()};
      assert(algebra::approx_equal(matrix44(m * m_inv),
                                   matrix44(matrix44::Identity()),
                                   16.f * epsilon, 1e-6f));
    }
  }

  /// Constructor with arguments: matrix as std::array of scalar
  ///
  /// @param ma is the full 4x4 matrix as a 16 array
  DETRAY_HOST_DEVICE
  explicit transform3(const array_type<16> &ma) {
    _data.matrix() << ma[0], ma[1], ma[2], ma[3], ma[4], ma[5], ma[6], ma[7],
        ma[8], ma[9], ma[10], ma[11], ma[12], ma[13], ma[14], ma[15];

    _data_inv = _data.inverse();
  }

  /// Default constructor: set contents to identity matrices
  DETRAY_HOST_DEVICE
  transform3() {
    _data.setIdentity();
    _data_inv.setIdentity();
  }

  /// Copy an existing transform.
  /// @param rhs Transform whose matrices are copied.
  transform3(const transform3 &rhs) = default;
  ~transform3() = default;

  /// Compare the forward transforms.
  /// @param rhs transform to compare with
  /// @return whether the transforms have equal spatial mappings
  DETRAY_HOST_DEVICE
  constexpr bool operator==(const transform3 &rhs) const {
    return (_data.isApprox(rhs._data));
  }

  /// Rotate a vector into / from a frame
  ///
  /// @param m is the rotation matrix
  /// @param v vector to rotate
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
  /// @param v vector to rotate
  /// @return the rotated vector
  template <typename derived_type>
    requires(Eigen::MatrixBase<derived_type>::RowsAtCompileTime == 3 &&
             Eigen::MatrixBase<derived_type>::ColsAtCompileTime == 1)
  DETRAY_HOST_DEVICE static constexpr auto rotate(
      const Eigen::Transform<scalar_type, 3, Eigen::Affine> &m,
      const Eigen::MatrixBase<derived_type> &v) {
    return m.matrix().template block<3, 3>(0, 0) * v;
  }

  /// Return the upper-left 3x3 rotation block.
  /// @return rotation matrix
  DETRAY_HOST_DEVICE
  constexpr auto rotation() const {
    return _data.matrix().template block<3, 3>(0, 0);
  }

  /// Return the local x axis in global coordinates.
  /// @return x-axis vector
  DETRAY_HOST_DEVICE
  constexpr point3 x() const {
    return _data.matrix().template block<3, 1>(0, 0);
  }

  /// Return the local y axis in global coordinates.
  /// @return y-axis vector
  DETRAY_HOST_DEVICE
  constexpr point3 y() const {
    return _data.matrix().template block<3, 1>(0, 1);
  }

  /// Return the local z axis in global coordinates.
  /// @return z-axis vector
  DETRAY_HOST_DEVICE
  constexpr point3 z() const {
    return _data.matrix().template block<3, 1>(0, 2);
  }

  /// Return the frame origin in global coordinates.
  /// @return translation vector
  DETRAY_HOST_DEVICE
  constexpr point3 translation() const {
    return _data.matrix().template block<3, 1>(0, 3);
  }

  /// Return the local-to-global transformation matrix.
  /// @return forward transformation matrix
  DETRAY_HOST_DEVICE
  constexpr const matrix44 &matrix() const { return _data.matrix(); }

  /// Return the global-to-local transformation matrix.
  /// @return inverse transformation matrix
  DETRAY_HOST_DEVICE
  constexpr const matrix44 &matrix_inverse() const {
    return _data_inv.matrix();
  }

  /// Transform a point from local plane coordinates to global coordinates.
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

  /// Transform a point from local coordinates to global coordinates.
  /// @param v point in local coordinates
  /// @return point in global coordinates
  template <typename derived_type>
    requires(Eigen::MatrixBase<derived_type>::RowsAtCompileTime == 3 &&
             Eigen::MatrixBase<derived_type>::ColsAtCompileTime == 1)
  DETRAY_HOST_DEVICE constexpr auto point_to_global(
      const Eigen::MatrixBase<derived_type> &v) const {
    return (_data * v);
  }

  /// Transform a global vector to local coordinates.
  /// @param v vector in global coordinates
  /// @return vector in local coordinates
  template <typename derived_type>
    requires(Eigen::MatrixBase<derived_type>::RowsAtCompileTime == 3 &&
             Eigen::MatrixBase<derived_type>::ColsAtCompileTime == 1)
  DETRAY_HOST_DEVICE constexpr auto point_to_local(
      const Eigen::MatrixBase<derived_type> &v) const {
    return (_data_inv * v);
  }

  /// Transform a local vector to global coordinates.
  /// @param v vector in local coordinates
  /// @return vector in global coordinates
  template <typename derived_type>
    requires(Eigen::MatrixBase<derived_type>::RowsAtCompileTime == 2 &&
             Eigen::MatrixBase<derived_type>::ColsAtCompileTime == 1)
  DETRAY_HOST_DEVICE constexpr auto vector_to_global(
      const Eigen::MatrixBase<derived_type> &v) const {
    return (_data.linear().matrix().template block<3, 2>(0, 0) * v);
  }

  /// Transform a local vector to global coordinates.
  /// @param v vector in local coordinates
  /// @return vector in global coordinates
  template <typename derived_type>
    requires(Eigen::MatrixBase<derived_type>::RowsAtCompileTime == 3 &&
             Eigen::MatrixBase<derived_type>::ColsAtCompileTime == 1)
  DETRAY_HOST_DEVICE constexpr auto vector_to_global(
      const Eigen::MatrixBase<derived_type> &v) const {
    return (_data.linear() * v);
  }

  /// Transform a global vector to local coordinates.
  /// @param v vector in global coordinates
  /// @return vector in local coordinates
  template <typename derived_type>
    requires(Eigen::MatrixBase<derived_type>::RowsAtCompileTime == 3 &&
             Eigen::MatrixBase<derived_type>::ColsAtCompileTime == 1)
  DETRAY_HOST_DEVICE constexpr auto vector_to_local(
      const Eigen::MatrixBase<derived_type> &v) const {
    return (_data_inv.linear() * v);
  }
};  // struct transform3

}  // namespace detray::algebra::eigen::math
