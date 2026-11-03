// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>

namespace Acts::detail {

/// Compute the quaternion that rotates vector @p a onto vector @p b.
///
/// This wraps @c Eigen::Quaternion::setFromTwoVectors, whose implementation
/// instantiates a (fixed-size) @c Eigen::JacobiSVD. That SVD instantiation is
/// very expensive in compiler memory (~0.5 GB per translation unit). Declaring
/// this helper here and providing the definition out-of-line (explicitly
/// instantiated for @c float and @c double in QuaternionFromTwoVectors.cpp)
/// means the SVD is instantiated once in the Acts core library instead of in
/// every translation unit that needs it.
///
/// @param a the source vector
/// @param b the target vector
/// @return the quaternion rotating @p a onto @p b
Eigen::Quaternion<float> quaternionFromTwoVectors(
    const Eigen::Matrix<float, 3, 1>& a, const Eigen::Matrix<float, 3, 1>& b);

/// Compute the quaternion that rotates vector @p a onto vector @p b.
///
/// This wraps @c Eigen::Quaternion::setFromTwoVectors, whose implementation
/// instantiates a (fixed-size) @c Eigen::JacobiSVD. That SVD instantiation is
/// very expensive in compiler memory (~0.5 GB per translation unit). Declaring
/// this helper here and providing the definition out-of-line (explicitly
/// instantiated for @c float and @c double in QuaternionFromTwoVectors.cpp)
/// means the SVD is instantiated once in the Acts core library instead of in
/// every translation unit that needs it.
///
/// @param a the source vector
/// @param b the target vector
/// @return the quaternion rotating @p a onto @p b
Eigen::Quaternion<double> quaternionFromTwoVectors(
    const Eigen::Matrix<double, 3, 1>& a, const Eigen::Matrix<double, 3, 1>& b);

}  // namespace Acts::detail
