// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/Intersection.hpp"

/// @brief Helpers for planar surfaces that share the same maths
namespace Acts::PlanarHelper {

/// @brief Intersect a line in 3D space with a plane represented by the Hesse-Normal form
/// @param linePos: Arbitrary point on the line to intersect
/// @param lineDir: Direction of the line to intersect
/// @param planeNorm: Normal vector of the plane
/// @param offSet: Offset to move the plane along the normal vector
inline Intersection3D intersectPlane(const Vector3& linePos,
                                     const Vector3& lineDir,
                                     const Vector3& planeNorm,
                                     const double offset) {
  /// Use the formula: <P, N> - C = 0
  ///  --> insert line equation: <A + lambda * B, N> - C = 0
  ///  --> lambda = (C - <A,N>)/ <N, B> */
  const double normDot = planeNorm.dot(lineDir);
  if (std::abs(normDot) < std::numeric_limits<double>::epsilon()) {
    return Intersection3D::Invalid();
  }
  const double path = (offset - linePos.dot(planeNorm)) / normDot;
  return Intersection3D{linePos + path * lineDir, path,
                        IntersectionStatus::onSurface};
}
/// @brief Intersect a line in 3D space with a plane represented by the Hesse-Normal form
/// @param linePos: Arbitrary point on the line to intersect
/// @param lineDir: Direction of the line to intersect
/// @param planeNorm: Normal vector of the plane
/// @param planePoint: Point on the plane
inline Intersection3D intersectPlane(const Vector3& linePos,
                                     const Vector3& lineDir,
                                     const Vector3& planeNorm,
                                     const Vector3& planePoint) {
  return intersectPlane(linePos, lineDir, planeNorm, planePoint.dot(planeNorm));
}

/// Intersection with a planar surface
///
/// @param transform The 3D affine transform that places the surface
/// @param position The starting position for the intersection
/// @param direction The starting direction for the intersection
///
/// @return The intersection
inline Intersection3D intersect(const Transform3& transform,
                                const Vector3& position,
                                const Vector3& direction, double tolerance) {
  // Get the matrix from the transform (faster access)
  const auto& tMatrix = transform.matrix();
  const Vector3 pnormal = tMatrix.block<3, 1>(0, 2).transpose();
  const Vector3 pcenter = tMatrix.block<3, 1>(0, 3).transpose();
  // It is solvable, so go on
  double denom = direction.dot(pnormal);
  if (denom == 0) {
    // The line is parallel to the plane, hence no intersection
    return Intersection3D::Invalid();
  }
  // Translate that into a path
  double path = (pnormal.dot((pcenter - position))) / (denom);
  // Is valid hence either on surface or reachable
  IntersectionStatus status = std::abs(path) < std::abs(tolerance)
                                  ? IntersectionStatus::onSurface
                                  : IntersectionStatus::reachable;
  // Return the intersection
  return Intersection3D{(position + path * direction), path, status};
}

}  // namespace Acts::PlanarHelper
