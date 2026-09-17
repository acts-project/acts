// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

namespace Acts {
/// Invert a placement, i.e. a transform that only rotates and translates.
/// `Transform3` is an `Eigen::Affine`, so its own `inverse()` is a general one.
///
/// @param transform The placement to invert, whose linear part must be a rotation
/// @return The inverse transform
inline Transform3 inverseTransform(const Transform3& transform) {
  assert(transform.linear().isUnitary() &&
         "a placement transform must not scale or shear");
  return transform.inverse(Eigen::Isometry);
}
/// Rotate the coordinate system by an angle around the x-axis
/// @param angle: The rotation angle around x.
/// @return Transform that rotates the coordinate system by the parsed angle around the x-axis
inline Transform3 getRotateX3D(double angle) {
  return Transform3{AngleAxis3{angle, Vector3::UnitX()}};
}
/// Rotate the coordinate system by an angle around the z-axis
/// @param angle: The rotation angle around y.
/// @return Transform that rotates the coordinate system by the parsed angle around the y-axis
inline Transform3 getRotateY3D(double angle) {
  return Transform3{AngleAxis3{angle, Vector3::UnitY()}};
}
/// Rotate the coordinate system by an angle around the z-axis
/// @param angle: The rotation angle around z.
/// @return Transform that rotates the coordinate system by the parsed angle around the z-axis
inline Transform3 getRotateZ3D(double angle) {
  return Transform3{AngleAxis3{angle, Vector3::UnitZ()}};
}
/// Returns a shift transformation along the x-axis
/// @param X: The value by which the coordinate system is shifted
/// @return Transform that shifts the coordinate system by the parsed amount in X
inline Transform3 getTranslateX3D(const double X) {
  return Transform3{Translation3{X * Vector3::UnitX()}};
}
/// Returns a shift transformation along the y-axis
/// @param Y: The value by which the coordinate system is shifted
/// @return Transform that shifts the coordinate system by the parsed amount in Y
inline Transform3 getTranslateY3D(const double Y) {
  return Transform3{Translation3{Y * Vector3::UnitY()}};
}
/// Returns a shift transformation along the z-axis
/// @param Z: The value by which the coordinate system is shifted
/// @return Transform that shifts the coordinate system by the parsed amount in Z
inline Transform3 getTranslateZ3D(const double Z) {
  return Transform3{Translation3{Z * Vector3::UnitZ()}};
}

/// Returns a shift transformation for an arbitrary position tuple
/// @param X: Shift along the x-axis
/// @param Y: Shift along the y-axis
/// @param Z: Shift along the z-axis
/// @return Transform that shifts the coordinate system by an arbitrary point
inline Transform3 getTranslate3D(const double X, const double Y,
                                 const double Z) {
  return getTranslateX3D(X) * getTranslateY3D(Y) * getTranslateZ3D(Z);
}
/// Returns a shift transformation according a vector
/// @param v: The vector by which the coordinate system is shifted
/// @return Transform that shifts the coordinate system by an arbitrary point
inline Transform3 getTranslate3D(const Vector3& v) {
  return Transform3{Translation3{v}};
}
}  // namespace Acts
