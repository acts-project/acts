// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
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

#include "Acts/Utilities/Units.hpp"

#ifdef TRKDETDESCR_USEFLOATPRECISON
using TDD_real_t = float;
#else
using TDD_real_t = double;
#endif

#define TDD_max_bound_value 10e10

namespace Acts {

///  Tolerance for being on Surface
static const double s_onSurfaceTolerance = 10e-5 * units::_mm;

/// Tolerance for not being within curvilinear projection
/// this allows using the same curvilinear frame to eta = 6,
/// validity tested with IntegrationTests/PropagationTest
static const double s_curvilinearProjTolerance = 0.999995;

/// @enum NavigationDirection
/// The navigation direciton is always with
/// respect to a given momentum or direction
enum NavigationDirection : int { backward = -1, anyDirection = 0, forward = 1 };

///  This is a steering enum to tell which material update stage:
/// - preUpdate  : update on approach of a surface
/// - fullUpdate : update when passing a surface
/// - postUpdate : update when leaving a surface
enum MaterialUpdateStage : int {
  preUpdate  = -1,
  fullUpdate = 0,
  postUpdate = 1
};

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

using Rotation3D       = Eigen::Quaternion<double>;
using Translation3D    = Eigen::Translation<double, 3>;
using AngleAxis3D      = Eigen::AngleAxisd;
using Transform3D      = Eigen::Affine3d;
using Vector3D         = Eigen::Matrix<double, 3, 1>;
using Vector2D         = Eigen::Matrix<double, 2, 1>;
using RotationMatrix3D = Eigen::Matrix<double, 3, 3>;

/// axis defintion element for code readability
/// - please use these for access to the member variables if needed, e.g.
///     double z  = position[Acts::eZ];
///     double px = momentum[Acts::ePX];
///
enum AxisDefs : int {
  // position access
  eX = 0,
  eY = 1,
  eZ = 2,
  // momentum access
  ePX = 0,
  ePY = 1,
  ePZ = 2
};

}  // namespace Acts
