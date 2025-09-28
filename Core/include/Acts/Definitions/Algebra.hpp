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

namespace Acts {

/// @defgroup acts-algebra-types Vector/matrix types with a common scalar type
///
/// These are the default vector/matrix types that should be used throughout the
/// codebase. They all use the common Acts scalar type but support variable size
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
/// `ActsDynamicMatrix` instead.

/// @brief Fixed-size vector type for N-dimensional vectors
/// @tparam kSize The dimension of the vector
template <unsigned int kSize>
using ActsVector = Eigen::Matrix<double, kSize, 1>;

/// @brief Fixed-size matrix type for NxM matrices
/// @tparam kRows Number of rows
/// @tparam kCols Number of columns
template <unsigned int kRows, unsigned int kCols>
using ActsMatrix = Eigen::Matrix<double, kRows, kCols>;

/// @brief Fixed-size square matrix type for NxN matrices
/// @tparam kSize The dimension of the square matrix
template <unsigned int kSize>
using ActsSquareMatrix = Eigen::Matrix<double, kSize, kSize>;

/// @brief Dynamic-sized vector type
using ActsDynamicVector = Eigen::Matrix<double, Eigen::Dynamic, 1>;

/// @brief Dynamic-sized matrix type
using ActsDynamicMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

/// @defgroup coordinates-types Fixed-size vector/matrix types for coordinates
///
/// These predefined types should always be used when handling coordinate
/// vectors in different coordinate systems, i.e. on surfaces (2d), spatial
/// position (3d), or space-time (4d).
///

// coordinate vectors
/// @brief 2-dimensional vector type for 2D coordinates
using Vector2 = ActsVector<2>;
/// @brief 3-dimensional vector type for e.g. spatial coordinates and momenta
using Vector3 = ActsVector<3>;
/// @brief 4-dimensional vector type for space-time coordinates
using Vector4 = ActsVector<4>;

// square matrices e.g. for coordinate covariance matrices
/// @brief 2x2 square matrix type, typically used for 2D coordinate covariance
using SquareMatrix2 = ActsSquareMatrix<2>;
/// @brief 3x3 square matrix type, typically used for 3D coordinate covariance
using SquareMatrix3 = ActsSquareMatrix<3>;
/// @brief 4x4 square matrix type, typically used for 4D coordinate covariance
using SquareMatrix4 = ActsSquareMatrix<4>;

// pure translation transformations
/// @brief 2D translation transformation
using Translation2 = Eigen::Translation<double, 2>;
/// @brief 3D translation transformation
using Translation3 = Eigen::Translation<double, 3>;

/// @brief 2D rotation matrix
using RotationMatrix2 = SquareMatrix2;
/// @brief 3D rotation matrix
using RotationMatrix3 = SquareMatrix3;

// pure rotation defined by a rotation angle around a rotation axis
/// @brief Rotation defined by an angle around a rotation axis in 3D
using AngleAxis3 = Eigen::AngleAxis<double>;

/// @brief 2D affine transformation stored as a compact 2x3 matrix
using Transform2 = Eigen::Transform<double, 2, Eigen::AffineCompact>;
/// @brief 3D affine transformation stored as a 4x4 matrix
using Transform3 = Eigen::Transform<double, 3, Eigen::Affine>;

/// Tolerance for transform equivalence checks
constexpr double s_transformEquivalentTolerance = 1e-9;

}  // namespace Acts
