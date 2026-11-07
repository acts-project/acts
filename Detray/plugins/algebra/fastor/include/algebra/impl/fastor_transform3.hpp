// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "algebra/impl/detail/fastor_matrix_wrapper.hpp"
#include "algebra/impl/detail/fastor_vector_wrapper.hpp"
#include "algebra/impl/fastor_matrix.hpp"
#include "detray/algebra/utils/approximately_equal.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

// Fastor include(s).
#ifdef _MSC_VER
#pragma warning(disable : 4244 4701 4702)
#endif  // MSVC
#include <Fastor/Fastor.h>
#ifdef _MSC_VER
#pragma warning(default : 4244 4701 4702)
#endif  // MSVC

// System include(s).
#include <cassert>
#include <concepts>
#include <cstddef>
#include <limits>

namespace detray::algebra::fastor::math {

/// Transform wrapper class to ensure standard API towards the Fastor library
template <concepts::scalar scalar_t>
struct transform3 {
  /// @name Type definitions for the struct
  /// @{

  /// Scalar type used by the transform
  using scalar_type = scalar_t;

  /// Array type used by the transform
  template <std::size_t N>
  using array_type = fastor::Vector<scalar_t, N>;

  /// Vector and point types
  using vector2 = array_type<2>;
  /// Three-component vector type used for spatial quantities.
  using vector3 = array_type<3>;
  /// Two-dimensional point representation.
  using point2 = vector2;
  /// Three-dimensional point representation.
  using point3 = vector3;

  /// 4x4 matrix type
  using matrix44 = fastor::Matrix<scalar_type, 4UL, 4UL>;

  /// Helper type to cast this to another floating point precision
  template <concepts::scalar o_scalar_t>
  using other_type = transform3<o_scalar_t>;

  /// @}

  /// @name Data objects
  /// @{

  matrix44 _data;      ///< Local-to-global transformation matrix.
  matrix44 _data_inv;  ///< Global-to-local transformation matrix.

  /// @}

  /// Default constructor: identity
  DETRAY_HOST_DEVICE
  constexpr transform3() {
    _data.eye();
    _data_inv.eye();
  }

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
    // The matrix needs to be initialized to the identity matrix first. We
    // only modify the top 4x3 portion of the matrix, so it doesn't matter
    // what values it initially had. However, the bottom row is required to
    // have the values of [0, 0, 0, 1], so that's why we set `_data` to the
    // identity matrix first.
    _data.eye2();

    _data(Fastor::fseq<0, 3>(), 0) = x;
    _data(Fastor::fseq<0, 3>(), 1) = y;
    _data(Fastor::fseq<0, 3>(), 2) = z;
    _data(Fastor::fseq<0, 3>(), 3) = t;

    if (get_inverse) {
      _data_inv = Fastor::inverse(_data);
    }
  }

  /// Constructor with arguments: t, z, x
  ///
  /// @param t the translation (or origin of the new frame)
  /// @param z the z axis of the new frame, normal vector for planes
  /// @param x the x axis of the new frame
  /// @param get_inverse whether to compute and store the inverse
  DETRAY_HOST
  transform3(const vector3 &t, const vector3 &z, const vector3 &x,
             bool get_inverse = true)
      : transform3(t, x, Fastor::cross(z, x), z, get_inverse) {}

  /// Constructor with arguments: translation
  ///
  /// @param t is the translation
  DETRAY_HOST
  explicit transform3(const vector3 &t) {
    // The matrix needs to be initialized to the identity matrix first. In
    // this case, the `transform3` requires `_data` to look just like an
    // identity matrix except for the third column, which is the one we are
    // modifying here.
    _data.eye2();

    _data(Fastor::fseq<0, 3>(), 3) = t;

    _data_inv = Fastor::inverse(_data);
  }

  /// Constructor with arguments: matrix
  ///
  /// @param m is the full 4x4 matrix
  DETRAY_HOST
  explicit transform3(const matrix44 &m) : _data{m} {
    _data_inv = Fastor::inverse(_data);
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

      [[maybe_unused]] matrix44 identity_matrix;
      identity_matrix.eye2();
      assert(algebra::approx_equal(m * m_inv, identity_matrix, 16.f * epsilon,
                                   1e-6f));
    }
  }

  /// Constructor with arguments: matrix as Fastor::Tensor<scalar_t, 16> of
  /// scalars
  ///
  /// @param ma is the full 4x4 matrix as a 16-element array
  DETRAY_HOST
  explicit transform3(const array_type<16> &ma) : _data{ma} {
    _data_inv = Fastor::inverse(_data);
  }

  /// Copy an existing transform.
  /// @param rhs Transform whose matrices are copied.
  transform3(const transform3 &rhs) = default;
  ~transform3() = default;

  /// Equality operator
  /// @param rhs transform to compare with
  /// @return whether the transforms have equal spatial mappings
  DETRAY_HOST
  constexpr bool operator==(const transform3 &rhs) const {
    return Fastor::isequal(_data, rhs._data);
  }

  /// This method retrieves the rotation of a transform
  /// @return rotation matrix
  DETRAY_HOST
  constexpr auto rotation() const {
    return Fastor::Tensor<scalar_t, 3, 3>(
        _data(Fastor::fseq<0, 3>(), Fastor::fseq<0, 3>()));
  }

  /// This method retrieves x axis
  /// @return x-axis vector
  DETRAY_HOST_DEVICE
  constexpr point3 x() const { return _data(Fastor::fseq<0, 3>(), 0); }

  /// This method retrieves y axis
  /// @return y-axis vector
  DETRAY_HOST_DEVICE
  constexpr point3 y() const { return _data(Fastor::fseq<0, 3>(), 1); }

  /// This method retrieves z axis
  /// @return z-axis vector
  DETRAY_HOST_DEVICE
  constexpr point3 z() const { return _data(Fastor::fseq<0, 3>(), 2); }

  /// This method retrieves the translation of a transform
  /// @return translation vector
  DETRAY_HOST
  constexpr vector3 translation() const {
    return _data(Fastor::fseq<0, 3>(), 3);
  }

  /// This method retrieves the 4x4 matrix of a transform
  /// @return forward transformation matrix
  DETRAY_HOST
  constexpr matrix44 matrix() const { return _data; }

  /// This method retrieves the 4x4 matrix of an inverse transform
  /// @return inverse transformation matrix
  DETRAY_HOST
  constexpr matrix44 matrix_inverse() const { return _data_inv; }

  /// This method transform from a point from the local 3D cartesian frame to
  /// the global 3D cartesian frame
  /// @param v point in local coordinates
  /// @return point in global coordinates
  DETRAY_HOST
  constexpr point3 point_to_global(const point2 &v) const {
    Fastor::Tensor<scalar_type, 4> vector_4;
    vector_4(Fastor::fseq<0, 2>()) = v;
    vector_4[2] = static_cast<scalar_type>(0);
    vector_4[3] = static_cast<scalar_type>(1);
    return Fastor::Tensor<scalar_type, 3>(
        Fastor::matmul(_data, vector_4)(Fastor::fseq<0, 3>()));
  }

  /// This method transform from a point from the local 3D cartesian frame to
  /// the global 3D cartesian frame
  /// @param v point in local coordinates
  /// @return point in global coordinates
  DETRAY_HOST
  constexpr point3 point_to_global(const point3 &v) const {
    Fastor::Tensor<scalar_type, 4> vector_4;
    vector_4(Fastor::fseq<0, 3>()) = v;
    vector_4[3] = static_cast<scalar_type>(1);
    return Fastor::Tensor<scalar_type, 3>(
        Fastor::matmul(_data, vector_4)(Fastor::fseq<0, 3>()));
  }

  /// Transform a global point to local coordinates.
  /// @param v point in global coordinates
  /// @return point in local coordinates
  DETRAY_HOST
  constexpr point3 point_to_local(const point3 &v) const {
    Fastor::Tensor<scalar_type, 4> vector_4;
    vector_4(Fastor::fseq<0, 3>()) = v;
    vector_4[3] = static_cast<scalar_type>(1);
    return Fastor::Tensor<scalar_type, 3>(
        Fastor::matmul(_data_inv, vector_4)(Fastor::fseq<0, 3>()));
  }

  /// This method transform from a vector from the local 2D cartesian frame to
  /// the global 3D cartesian frame
  /// @param v vector in local plane coordinates
  /// @return vector in global coordinates
  DETRAY_HOST
  constexpr point3 vector_to_global(const vector2 &v) const {
    Fastor::Tensor<scalar_type, 4> vector_4;
    vector_4(Fastor::fseq<0, 2>()) = v;
    vector_4[2] = static_cast<scalar_type>(0);
    vector_4[3] = static_cast<scalar_type>(0);
    return Fastor::Tensor<scalar_type, 3>(
        Fastor::matmul(_data, vector_4)(Fastor::fseq<0, 3>()));
  }

  /// This method transform from a vector from the local 3D cartesian frame to
  /// the global 3D cartesian frame
  /// @param v vector in local coordinates
  /// @return vector in global coordinates
  DETRAY_HOST
  constexpr point3 vector_to_global(const vector3 &v) const {
    Fastor::Tensor<scalar_type, 4> vector_4;
    vector_4(Fastor::fseq<0, 3>()) = v;
    vector_4[3] = static_cast<scalar_type>(0);
    return Fastor::Tensor<scalar_type, 3>(
        Fastor::matmul(_data, vector_4)(Fastor::fseq<0, 3>()));
  }

  /// This method transform from a vector from the global 3D cartesian frame
  /// into the local 3D cartesian frame
  /// @param v vector in global coordinates
  /// @return vector in local coordinates
  DETRAY_HOST
  constexpr point3 vector_to_local(const vector3 &v) const {
    Fastor::Tensor<scalar_type, 4> vector_4;
    vector_4(Fastor::fseq<0, 3>()) = v;
    vector_4[3] = static_cast<scalar_type>(0);
    return Fastor::Tensor<scalar_type, 3>(
        Fastor::matmul(_data_inv, vector_4)(Fastor::fseq<0, 3>()));
  }
};  // struct transform3

}  // namespace detray::algebra::fastor::math
