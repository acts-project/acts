// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_Definitions_H
#define ACTS_Definitions_H

// API addons
#define EIGEN_MATRIXBASE_PLUGIN "ACTS/Utilities/detail/MatrixBasePlugin.hpp"
#define EIGEN_MATRIX_PLUGIN "ACTS/Utilities/detail/MatrixPlugin.hpp"
#define EIGEN_TRANSFORM_PLUGIN "ACTS/Utilities/detail/TransformPlugin.hpp"

// external include(s)
#pragma GCC diagnostic push
#if defined(__USE_GNU)
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#endif
#include <Eigen/Dense>
#pragma GCC diagnostic pop

#ifdef TRKDETDESCR_USEFLOATPRECISON
typedef float TDD_real_t;
#else
typedef double TDD_real_t;
#endif

#define TDD_max_bound_value 10e10

namespace Acts {
/** Tolerance for being on Surface */
static const double s_onSurfaceTolerance = 10e-5;

/** @enum PropDirection
  PropDirection, enum for direction of the propagation.

  */
enum PropDirection {
  alongMomentum    = 1,
  oppositeMomentum = -1,
  anyDirection     = 0,
  mappingMode      = 2
};

/// @brief simple enum for searching surfaces
enum SearchDirection {
  outside            = 1,
  inside             = -1,
  bothway            = 0,
  undefinedDirection = 0
};

/** This is a steering enum to tell which material update stage:
   - preUpdate  : when reaching a layer before layer is resolved
   - fullUpdate : just pass through the layer
   - postUpdate : when leaving the layer
 */
enum MaterialUpdateStage { preUpdate = -1, fullUpdate = 0, postUpdate = 1 };

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

/** elment for code readability
    - please use these for access to the member variables if needed, e.g.
        double z  = position[Acts::eZ];
        double px = momentum[Acts::ePX];
*/
enum AxisDefs {
  // position access
  eX = 0,
  eY = 1,
  eZ = 2,
  // momentum access
  ePX = 0,
  ePY = 1,
  ePZ = 2
};

typedef Eigen::Quaternion<double> Rotation3D;
typedef Eigen::Translation<double, 3> Translation3D;
typedef Eigen::AngleAxisd AngleAxis3D;
typedef Eigen::Affine3d   Transform3D;
typedef Eigen::Matrix<double, 3, 1> Vector3D;
typedef Eigen::Matrix<double, 2, 1> Vector2D;
typedef Eigen::Matrix<double, 3, 3> RotationMatrix3D;

}  // end of namespace Acts

#endif
