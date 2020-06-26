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
#include <Eigen/Dense>
#pragma GCC diagnostic pop
#else
#include <Eigen/Dense>
#endif

namespace Acts {

/// Tolerance for bein numerical equal for geometry building
static constexpr double s_epsilon = 3 * std::numeric_limits<double>::epsilon();

/// Tolerance for being on Surface
///
/// @note This is intentionally given w/o an explicit unit to avoid having
///       to include the units header unneccessarily. With the native length
///       unit of mm this corresponds to 0.1um.
static constexpr double s_onSurfaceTolerance = 1e-4;

/// Tolerance for not being within curvilinear projection
/// this allows using the same curvilinear frame to eta = 6,
/// validity tested with IntegrationTests/PropagationTest
static constexpr double s_curvilinearProjTolerance = 0.999995;

/// @enum NavigationDirection
/// The navigation direciton is always with
/// respect to a given momentum or direction
enum NavigationDirection : int { backward = -1, forward = 1 };

///  This is a steering enum to tell which material update stage:
/// - preUpdate  : update on approach of a surface
/// - fullUpdate : update when passing a surface
/// - postUpdate : update when leaving a surface
enum MaterialUpdateStage : int {
  preUpdate = -1,
  fullUpdate = 0,
  postUpdate = 1
};

/// @enum NoiseUpdateMode to tell how to deal with noise term in covariance
/// transport
/// - removeNoise: subtract noise term
/// - addNoise: add noise term
enum NoiseUpdateMode : int { removeNoise = -1, addNoise = 1 };

// Eigen definitions
template <typename T, unsigned int rows, unsigned int cols>
using ActsMatrix = Eigen::Matrix<T, rows, cols>;

template <unsigned int rows, unsigned int cols>
using ActsMatrixD = ActsMatrix<double, rows, cols>;

template <unsigned int rows, unsigned int cols>
using ActsMatrixF = ActsMatrix<float, rows, cols>;

template <typename T, unsigned int rows>
using ActsSymMatrix = Eigen::Matrix<T, rows, rows>;

template <unsigned int rows>
using ActsSymMatrixD = ActsSymMatrix<double, rows>;

template <unsigned int rows>
using ActsSymMatrixF = ActsSymMatrix<float, rows>;

template <typename T, unsigned int rows>
using ActsVector = Eigen::Matrix<T, rows, 1>;

template <unsigned int rows>
using ActsVectorD = ActsVector<double, rows>;

template <unsigned int rows>
using ActsVectorF = ActsVector<float, rows>;

template <typename T, unsigned int cols>
using ActsRowVector = Eigen::Matrix<T, 1, cols>;

template <unsigned int cols>
using ActsRowVectorD = ActsRowVector<double, cols>;

template <unsigned int cols>
using ActsRowVectorF = ActsRowVector<float, cols>;

template <typename T>
using ActsMatrixX = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

using ActsMatrixXd = ActsMatrixX<double>;
using ActsMatrixXf = ActsMatrixX<float>;

template <typename T>
using ActsVectorX = Eigen::Matrix<T, Eigen::Dynamic, 1>;

using ActsVectorXd = ActsVectorX<double>;
using ActsVectorXf = ActsVectorX<float>;

template <typename T>
using ActsRowVectorX = Eigen::Matrix<T, 1, Eigen::Dynamic>;

using ActsRowVectorXd = ActsRowVectorX<double>;
using ActsRowVectorXf = ActsRowVectorX<float>;

/// @defgroup fixed-algebra-types Fixed-size vector/matrix e.g. for coordinates
///
/// These predefined types should always be used when handling the coordinate
/// vectors in different coordinate systems, i.e. local (2d), global spatial
/// (3d), or global space-time (4d).
///
/// @{

// coordinate vectors
using Vector2F = ActsVector<float, 2>;
using Vector3F = ActsVector<float, 3>;
using Vector4F = ActsVector<float, 4>;
using Vector2D = ActsVector<double, 2>;
using Vector3D = ActsVector<double, 3>;
using Vector4D = ActsVector<double, 4>;
// symmetric matrices e.g. for coordinate covariance matrices
using SymMatrix2F = ActsSymMatrix<float, 2>;
using SymMatrix3F = ActsSymMatrix<float, 3>;
using SymMatrix4F = ActsSymMatrix<float, 4>;
using SymMatrix2D = ActsSymMatrix<double, 2>;
using SymMatrix3D = ActsSymMatrix<double, 3>;
using SymMatrix4D = ActsSymMatrix<double, 4>;

/// Components of coordinate vectors.
///
/// To be used to access coordinate components by named indices instead of magic
/// numbers. This must be a regular `enum` and not a scoped `enum class` to
/// allow implicit conversion to an integer. The enum value are thus visible
/// directly in `namespace Acts`.
///
/// This index enum is not user-configurable (in contrast e.g. to the track
/// parameter index enums) since it must be compatible with varying
/// dimensionality (2d-4d) and other access methods (`.{x,y,z}()` accessors).
enum CoordinateIndices : unsigned int {
  // generic position-like access
  ePos0 = 0,
  ePos1 = 1,
  ePos2 = 2,
  eTime = 3,
  // generic momentum-like access
  eMom0 = ePos0,
  eMom1 = ePos1,
  eMom2 = ePos2,
  eEnergy = eTime,
  // Cartesian spatial coordinates
  eX = ePos0,
  eY = ePos1,
  eZ = ePos2,
};

// pure translation transformations
using Translation2F = Eigen::Translation<float, 2>;
using Translation3F = Eigen::Translation<float, 3>;
using Translation4F = Eigen::Translation<float, 4>;
using Translation2D = Eigen::Translation<double, 2>;
using Translation3D = Eigen::Translation<double, 3>;
using Translation4D = Eigen::Translation<double, 4>;
// linear (rotation) matrices
using RotationMatrix2F = Eigen::Matrix<float, 2, 2>;
using RotationMatrix3F = Eigen::Matrix<float, 3, 3>;
using RotationMatrix4F = Eigen::Matrix<float, 4, 4>;
using RotationMatrix2D = Eigen::Matrix<double, 2, 2>;
using RotationMatrix3D = Eigen::Matrix<double, 3, 3>;
using RotationMatrix4D = Eigen::Matrix<double, 4, 4>;
// pure rotation transformations. only available in 2d and 3d
using Rotation2F = Eigen::Rotation2D<float>;
using Rotation3F = Eigen::Quaternion<float>;
using AngleAxis3F = Eigen::AngleAxis<float>;
using Rotation2D = Eigen::Rotation2D<double>;
using Rotation3D = Eigen::Quaternion<double>;
using AngleAxis3D = Eigen::AngleAxis<double>;
// combined affine transformations. types are chosen for better data alignment:
// - 2d affine compact stored as 2x3 matrix
// - 3d affine stored as 4x4 matrix
// - 4d affine compact stored as 4x5 matrix
using Transform2F = Eigen::Transform<float, 2, Eigen::AffineCompact>;
using Transform3F = Eigen::Transform<float, 3, Eigen::Affine>;
using Transform4F = Eigen::Transform<float, 4, Eigen::AffineCompact>;
using Transform2D = Eigen::Transform<double, 2, Eigen::AffineCompact>;
using Transform3D = Eigen::Transform<double, 3, Eigen::Affine>;
using Transform4D = Eigen::Transform<double, 4, Eigen::AffineCompact>;

/// @}

}  // namespace Acts
