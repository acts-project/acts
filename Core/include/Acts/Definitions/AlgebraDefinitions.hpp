// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// for GNU: ignore this specific warning, otherwise just include Eigen/Dense
#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#include <Eigen/Core>
#include <Eigen/Geometry>
#pragma GCC diagnostic pop
#else
#include <Eigen/Core>
#include <Eigen/Geometry>
#endif

namespace Acts {

/// This is the Scalar used for all EDM and Geometry based object
using AlgebraScalar = float;

template <typename T, unsigned int rows, unsigned int cols>
using ActsMatrix = Eigen::Matrix<T, rows, cols>;

template <unsigned int rows, unsigned int cols>
using ActsMatrixD = ActsMatrix<AlgebraScalar, rows, cols>;

template <typename T, unsigned int rows>
using ActsSymMatrix = Eigen::Matrix<T, rows, rows>;

template <unsigned int rows>
using ActsSymMatrixD = ActsSymMatrix<AlgebraScalar, rows>;

template <typename T, unsigned int rows>
using ActsVector = Eigen::Matrix<T, rows, 1>;

template <unsigned int rows>
using ActsVectorD = ActsVector<AlgebraScalar, rows>;

template <typename T, unsigned int cols>
using ActsRowVector = Eigen::Matrix<T, 1, cols>;

template <unsigned int cols>
using ActsRowVectorD = ActsRowVector<AlgebraScalar, cols>;

template <typename T>
using ActsMatrixX = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

using ActsMatrixXd = ActsMatrixX<AlgebraScalar>;

template <typename T>
using ActsVectorX = Eigen::Matrix<T, Eigen::Dynamic, 1>;

/// @defgroup fixed-algebra-types Fixed-size vector/matrix e.g. for coordinates
///
/// These predefined types should always be used when handling the coordinate
/// vectors in different coordinate systems, i.e. local (2d), global spatial
/// (3d), or global space-time (4d).
///
/// @{

// coordinate vectors
using Vector2D = ActsVector<AlgebraScalar, 2>;
using Vector3D = ActsVector<AlgebraScalar, 3>;
using Vector4D = ActsVector<AlgebraScalar, 4>;

// symmetric matrices e.g. for coordinate covariance matrices
using SymMatrix2D = ActsSymMatrix<AlgebraScalar, 2>;
using SymMatrix3D = ActsSymMatrix<AlgebraScalar, 3>;
using SymMatrix4D = ActsSymMatrix<AlgebraScalar, 4>;

// pure translation transformations
using Translation2D = Eigen::Translation<AlgebraScalar, 2>;
using Translation3D = Eigen::Translation<AlgebraScalar, 3>;

// linear (rotation) matrices
using RotationMatrix2D = Eigen::Matrix<AlgebraScalar, 2, 2>;
using RotationMatrix3D = Eigen::Matrix<AlgebraScalar, 3, 3>;

// pure rotation transformations. only available in 2d and 3d
using Rotation2D = Eigen::Rotation2D<AlgebraScalar>;
using Rotation3D = Eigen::Quaternion<AlgebraScalar>;
using AngleAxis3D = Eigen::AngleAxis<AlgebraScalar>;

// combined affine transformations. types are chosen for better data alignment:
// - 2d affine compact stored as 2x3 matrix
// - 3d affine stored as 4x4 matrix
using Transform2D = Eigen::Transform<AlgebraScalar, 2, Eigen::AffineCompact>;
using Transform3D = Eigen::Transform<AlgebraScalar, 3, Eigen::Affine>;

}  // namespace Acts
