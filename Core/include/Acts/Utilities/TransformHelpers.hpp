// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <cassert>

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

/// Create a rotation around the x-axis
/// @param angle The rotation angle
/// @return The rotation transform
inline Transform3 getRotateX3D(double angle) {
  return Transform3{AngleAxis3{angle, Vector3::UnitX()}};
}

/// Create a rotation around the y-axis
/// @param angle The rotation angle
/// @return The rotation transform
inline Transform3 getRotateY3D(double angle) {
  return Transform3{AngleAxis3{angle, Vector3::UnitY()}};
}

/// Create a rotation around the z-axis
/// @param angle The rotation angle
/// @return The rotation transform
inline Transform3 getRotateZ3D(double angle) {
  return Transform3{AngleAxis3{angle, Vector3::UnitZ()}};
}

/// Create a translation along the x-axis
/// @param x The shift along the x-axis
/// @return The translation transform
inline Transform3 getTranslateX3D(double x) {
  return Transform3{Translation3{x * Vector3::UnitX()}};
}

/// Create a translation along the y-axis
/// @param y The shift along the y-axis
/// @return The translation transform
inline Transform3 getTranslateY3D(double y) {
  return Transform3{Translation3{y * Vector3::UnitY()}};
}

/// Create a translation along the z-axis
/// @param z The shift along the z-axis
/// @return The translation transform
inline Transform3 getTranslateZ3D(double z) {
  return Transform3{Translation3{z * Vector3::UnitZ()}};
}

/// Create a translation from its components
/// @param x The shift along the x-axis
/// @param y The shift along the y-axis
/// @param z The shift along the z-axis
/// @return The translation transform
inline Transform3 getTranslate3D(double x, double y, double z) {
  return Transform3{Translation3{x, y, z}};
}

/// Create a translation from a vector
/// @param v The shift
/// @return The translation transform
inline Transform3 getTranslate3D(const Vector3& v) {
  return Transform3{Translation3{v}};
}

}  // namespace Acts
