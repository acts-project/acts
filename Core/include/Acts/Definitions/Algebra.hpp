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
/// The linear part is required to be orthogonal. This is what the geometry,
/// navigation and propagation code assumes throughout: inverses are taken as
/// transposes and local frame axes are read directly off the linear part.
/// Scaling or shearing a @c Transform3 is therefore a compile error - use
/// @ref Acts::AffineTransform3 where a general transformation is really meant.
using Transform3 = Eigen::Transform<double, 3, Eigen::Isometry>;

/// @brief 2D general affine transformation, see @ref Acts::AffineTransform3
using AffineTransform2 = Eigen::Transform<double, 2, Eigen::Affine>;
/// @brief 3D general affine transformation, allowing scaling and shearing
///
/// ACTS geometry cannot be built from these - they exist to name the
/// transformations that arrive from external geometry sources, which have to
/// be converted to a @ref Acts::Transform3 explicitly.
using AffineTransform3 = Eigen::Transform<double, 3, Eigen::Affine>;

/// Tolerance for transform equivalence checks
constexpr double s_transformEquivalentTolerance = 1e-9;

/// @brief Build a @ref Acts::Transform3 from a rotation and a translation
///
/// Maps a point as `p -> rotation * p + translation`, i.e. it rotates about
/// the origin first and translates afterwards. The translation is **not**
/// rotated - it is given in the frame that is mapped *into*, and ends up
/// verbatim in @c transform.translation().
///
/// Read as a frame rather than as a sequence, which avoids the question
/// entirely: the columns of @p rotation are the axes of the local frame and
/// @p translation is its origin, both expressed in the target frame.
///
/// This is the Eigen product `Translation3(translation) * rotation`, which
/// cannot be assigned to a @ref Acts::Transform3 because a plain matrix
/// carries no orthogonality guarantee. This function makes that guarantee
/// explicit and asserts it.
///
/// To translate along the *rotated* axes instead - the Eigen product
/// `rotation * Translation3(translation)` - rotate the translation yourself:
/// `makeTransform3(rotation, rotation * translation)`.
///
/// @param rotation The orthogonal linear part, i.e. the local frame axes
/// @param translation The local frame origin, in the target frame
/// @return The combined rigid transformation
inline Transform3 makeTransform3(const RotationMatrix3& rotation,
                                 const Vector3& translation) {
  assert((rotation * rotation.transpose())
             .isApprox(RotationMatrix3::Identity(),
                       s_transformEquivalentTolerance) &&
         "Transform3 requires an orthogonal rotation part");
  Transform3 transform = Transform3::Identity();
  transform.linear() = rotation;
  transform.translation() = translation;
  return transform;
}

/// @}

}  // namespace Acts
