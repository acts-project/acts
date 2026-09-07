// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// for GNU: ignore this specific warning, otherwise just include Eigen/Dense
#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#if __GNUC__ >= 12
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
#include <Eigen/Core>
#include <Eigen/Geometry>
#pragma GCC diagnostic pop
#else
#include <Eigen/Core>
#include <Eigen/Geometry>
#endif

#include <cassert>

namespace Acts {

/// @defgroup algebra_types Algebra types
///
/// These are the default vector/matrix types that should be used throughout the
/// codebase. They all use the common ACTS scalar type but support variable size
/// either at compile- or runtime.
///
/// Eigen does not have a distinct type for symmetric matrices. A typedef for
/// fixed-size matrices is still defined to simplify definition (one template
/// size vs two template size for generic matrices) and to clarify semantic
/// meaning in interfaces. It also ensures that the matrix is square. However,
/// the user is responsible for ensuring that the values are symmetric.
///
/// Without a distinct type for symmetric matrices, there is no way to provide
/// any conditions e.g. square size, for the dynamic-sized case. Consequently,
/// no dynamic-sized symmetric matrix type is defined. Use the
/// @ref Acts::DynamicMatrix instead.
///
/// @{

/// @brief Fixed-size vector type for N-dimensional vectors
/// @tparam kSize The dimension of the vector
template <unsigned int kSize>
using Vector = Eigen::Matrix<double, kSize, 1>;

/// @brief Fixed-size matrix type for NxM matrices
/// @tparam kRows Number of rows
/// @tparam kCols Number of columns
template <unsigned int kRows, unsigned int kCols>
using Matrix = Eigen::Matrix<double, kRows, kCols>;

/// @brief Fixed-size square matrix type for NxN matrices
/// @tparam kSize The dimension of the square matrix
template <unsigned int kSize>
using SquareMatrix = Eigen::Matrix<double, kSize, kSize>;

/// @brief Dynamic-sized vector type
using DynamicVector = Eigen::Matrix<double, Eigen::Dynamic, 1>;

/// @brief Dynamic-sized matrix type
using DynamicMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

/// @brief 2-dimensional vector type for 2D coordinates
using Vector2 = Vector<2>;
/// @brief 3-dimensional vector type for e.g. spatial coordinates and momenta
using Vector3 = Vector<3>;
/// @brief 4-dimensional vector type for space-time coordinates
using Vector4 = Vector<4>;

/// @brief 2x2 square matrix type, typically used for 2D coordinate covariance
using SquareMatrix2 = SquareMatrix<2>;
/// @brief 3x3 square matrix type, typically used for 3D coordinate covariance
using SquareMatrix3 = SquareMatrix<3>;
/// @brief 4x4 square matrix type, typically used for 4D coordinate covariance
using SquareMatrix4 = SquareMatrix<4>;

/// @brief 2D translation transformation
using Translation2 = Eigen::Translation<double, 2>;
/// @brief 3D translation transformation
using Translation3 = Eigen::Translation<double, 3>;

/// @brief 2D rotation matrix
using RotationMatrix2 = SquareMatrix2;
/// @brief 3D rotation matrix
using RotationMatrix3 = SquareMatrix3;

/// @brief Rotation defined by an angle around a rotation axis in 3D
using AngleAxis3 = Eigen::AngleAxis<double>;

/// @brief 2D rigid transformation, see @ref Acts::Transform3
using Transform2 = Eigen::Transform<double, 2, Eigen::Isometry>;
/// @brief 3D rigid transformation (rotation/reflection plus translation)
///
/// The linear part is required to be orthogonal: geometry, navigation and
/// propagation invert it by transposition and read the local frame axes off it
/// directly. Scaling or shearing is therefore a compile error - use
/// @ref Acts::AffineTransform3 for a general transformation.
using Transform3 = Eigen::Transform<double, 3, Eigen::Isometry>;

/// @brief 3D general affine transformation, allowing scaling and shearing
///
/// Names the transformations arriving from external geometry sources. They
/// have to be converted to a @ref Acts::Transform3 explicitly before ACTS
/// geometry can be built from them.
using AffineTransform3 = Eigen::Transform<double, 3, Eigen::Affine>;

/// Tolerance for transform equivalence checks
constexpr double s_transformEquivalentTolerance = 1e-9;

/// @brief Check whether a matrix is orthogonal, i.e. a rotation or reflection
///
/// This is the invariant @ref Acts::Transform3 relies on. Internal call sites
/// use @ref Acts::makeTransform3, which asserts it; call this directly to
/// reject a matrix coming from outside ACTS.
///
/// @param rotation The matrix to check
/// @return Whether the matrix is orthogonal within
///         @ref Acts::s_transformEquivalentTolerance
inline bool isOrthogonal(const RotationMatrix3& rotation) {
  return (rotation * rotation.transpose())
      .isApprox(RotationMatrix3::Identity(), s_transformEquivalentTolerance);
}

/// @brief Build a @ref Acts::Transform3 from a rotation and a translation
///
/// Maps `p -> rotation * p + translation`, i.e. the columns of @p rotation are
/// the local frame axes and @p translation its origin, both in the target
/// frame. Replaces the Eigen product `Translation3(translation) * rotation`,
/// which is affine and cannot be assigned to a @ref Acts::Transform3.
///
/// @param rotation The orthogonal linear part, i.e. the local frame axes
/// @param translation The local frame origin, in the target frame
/// @return The combined rigid transformation
inline Transform3 makeTransform3(const RotationMatrix3& rotation,
                                 const Vector3& translation) {
  assert(isOrthogonal(rotation) &&
         "Transform3 requires an orthogonal rotation part");
  Transform3 transform = Transform3::Identity();
  transform.linear() = rotation;
  transform.translation() = translation;
  return transform;
}

/// @}

}  // namespace Acts
