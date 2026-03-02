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

/// @brief Fixed-size vector type for N-dimensional vectors
/// @tparam kSize The dimension of the vector
/// @deprecated Use Vector instead
template <unsigned int kSize>
using ActsVector [[deprecated("Use Vector instead")]] = Vector<kSize>;
/// @brief Fixed-size matrix type for NxM matrices
/// @tparam kRows Number of rows
/// @tparam kCols Number of columns
/// @deprecated Use Matrix instead
template <unsigned int kRows, unsigned int kCols>
using ActsMatrix [[deprecated("Use Matrix instead")]] = Matrix<kRows, kCols>;
/// @brief Fixed-size square matrix type for NxN matrices
/// @tparam kSize The dimension of the square matrix
/// @deprecated Use SquareMatrix instead
template <unsigned int kSize>
using ActsSquareMatrix [[deprecated("Use SquareMatrix instead")]] =
    SquareMatrix<kSize>;
/// @brief Dynamic-sized vector type
/// @deprecated Use DynamicVector instead
using ActsDynamicVector [[deprecated("Use DynamicVector instead")]] =
    DynamicVector;
/// @brief Dynamic-sized matrix type
/// @deprecated Use DynamicMatrix instead
using ActsDynamicMatrix [[deprecated("Use DynamicMatrix instead")]] =
    DynamicMatrix;

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

/// @brief 2D affine transformation stored as a compact 2x3 matrix
using Transform2 = Eigen::Transform<double, 2, Eigen::AffineCompact>;
/// @brief 3D affine transformation stored as a 4x4 matrix
using Transform3 = Eigen::Transform<double, 3, Eigen::Affine>;

/// Tolerance for transform equivalence checks
constexpr double s_transformEquivalentTolerance = 1e-9;

/// @}

}  // namespace Acts
