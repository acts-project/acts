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

/// This is the Scalar used for all EDM and Geometry based object,
/// can be customised.
#ifdef ACTS_CUSTOM_SCALARTYPE
using ActsScalar = ACTS_CUSTOM_SCALARTYPE;
#else
using ActsScalar = double;
#endif

// Specify the algebra types

template <typename T, unsigned int rows, unsigned int cols>
using Matrix = Eigen::Matrix<T, rows, cols>;

template <typename T, unsigned int rows>
using SymMatrix = Eigen::Matrix<T, rows, rows>;

template <typename T, unsigned int rows>
using Vector = Eigen::Matrix<T, rows, 1>;

template <typename T, unsigned int cols>
using RowVector = Eigen::Matrix<T, 1, cols>;

template <typename T>
using DynamicMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

/// Specify the Acts types with the ActsScalar

template <unsigned int rows, unsigned int cols>
using ActsMatrix = Matrix<ActsScalar, rows, cols>;

template <unsigned int rows>
using ActsSymMatrix = SymMatrix<ActsScalar, rows>;

template <unsigned int rows>
using ActsVector = Vector<ActsScalar, rows>;

template <unsigned int cols>
using ActsRowVector = RowVector<ActsScalar, cols>;

using ActsDynamicMatrix = DynamicMatrix<ActsScalar>;

template <typename T>
using ActsDynamicVector = Eigen::Matrix<T, Eigen::Dynamic, 1>;

/// @defgroup fixed-algebra-types Fixed-size vector/matrix e.g. for coordinates
///
/// These predefined types should always be used when handling the coordinate
/// vectors in different coordinate systems, i.e. local (2d), global spatial
/// (3d), or global space-time (4d).
///
/// @{

// coordinate vectors
using Vector2 = ActsVector<2>;
using Vector3 = ActsVector<3>;
using Vector4 = ActsVector<4>;

// symmetric matrices e.g. for coordinate covariance matrices
using SymMatrix2 = ActsSymMatrix<2>;
using SymMatrix3 = ActsSymMatrix<3>;
using SymMatrix4 = ActsSymMatrix<4>;

// pure translation transformations
using Translation2 = Eigen::Translation<ActsScalar, 2>;
using Translation3 = Eigen::Translation<ActsScalar, 3>;

// linear (rotation) matrices
using RotationMatrix2 = ActsMatrix<2, 2>;
using RotationMatrix3 = ActsMatrix<3, 3>;

// pure rotation transformations. only available in 2d and 3d
using Rotation2 = Eigen::Rotation2D<ActsScalar>;
using Rotation3 = Eigen::Quaternion<ActsScalar>;
using AngleAxis3 = Eigen::AngleAxis<ActsScalar>;

// combined affine transformations. types are chosen for better data alignment:
// - 2d affine compact stored as 2x3 matrix
// - 3d affine stored as 4x4 matrix
using Transform2 = Eigen::Transform<ActsScalar, 2, Eigen::AffineCompact>;
using Transform3 = Eigen::Transform<ActsScalar, 3, Eigen::Affine>;

}  // namespace Acts
