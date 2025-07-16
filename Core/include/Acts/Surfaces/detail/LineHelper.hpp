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
#include "Acts/Utilities/MathHelpers.hpp"

namespace Acts::detail::LineHelper {
/// @brief Intersect two straight N-dimensional lines with each other or more generally
///        calculate the point of closest approach of the second line to the
///        first line.
/// @param linePosA: Arbitrary point on the first line
/// @param lineDirA: Direction of the first line (Unit-length)
/// @param linePosB: Arbitrary point on the second line
/// @param lineDirB: Direction of the second line (Unit-length)
template <int N>
inline Intersection<N> lineIntersect(const ActsVector<N>& linePosA,
                                     const ActsVector<N>& lineDirA,
                                     const ActsVector<N>& linePosB,
                                     const ActsVector<N>& lineDirB) {
  static_assert(N >= 2, "One dimensional intersect not sensible");
  /// Use the formula
  ///   A + lambda dirA  = B + mu dirB
  ///   (A-B) + lambda dirA = mu dirB
  ///   <A-B, dirB> + lambda <dirA,dirB> = mu
  ///    A + lambda dirA = B + (<A-B, dirB> + lambda <dirA,dirB>)dirB
  ///    <A-B,dirA> + lambda <dirA, dirA> = <A-B, dirB><dirA,dirB> +
  ///    lambda<dirA,dirB><dirA,dirB>
  ///  -> lambda = -(<A-B, dirA> - <A-B, dirB> * <dirA, dirB>) / (1-
  ///  <dirA,dirB>^2)
  ///  --> mu    =  (<A-B, dirB> - <A-B, dirA> * <dirA, dirB>) / (1-
  ///  <dirA,dirB>^2)
  const double dirDots = lineDirA.dot(lineDirB);
  const double divisor = (1. - square(dirDots));
  /// If the two directions are parallel to each other there's no way of
  /// intersection
  if (std::abs(divisor) < std::numeric_limits<double>::epsilon()) {
    return Intersection<N>::invalid();
  }
  const ActsVector<N> aMinusB = linePosA - linePosB;
  const double pathLength =
      (aMinusB.dot(lineDirB) - aMinusB.dot(lineDirA) * dirDots) / divisor;

  return Intersection<N>{linePosB + pathLength * lineDirB, pathLength,
                         IntersectionStatus::onSurface};
}
/// @brief Intersect the lines of two line surfaces using their respective transforms.
/// @param lineSurfTrf1: local -> global transform of the first surface
/// @param lineSurfTrf2: local -> global transform of the second surface
inline Intersection3D lineSurfaceIntersect(
    const Acts::Transform3& lineSurfTrf1,
    const Acts::Transform3& lineSurfTrf2) {
  return lineIntersect<3>(
      lineSurfTrf1.translation(), lineSurfTrf1.linear().col(2),
      lineSurfTrf2.translation(), lineSurfTrf2.linear().col(2));
}
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
    return Intersection3D::invalid();
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
/// @brief Calculates the signed distance between two lines in 3D space
/// @param linePosA: Arbitrary point on the first line
/// @param lineDirA: Direction of the first line (Unit-length)
/// @param linePosB: Arbitrary point on the second line
/// @param lineDirB: Direction of the second line (Unit-length)
inline double signedDistance(const Vector3& linePosA, const Vector3& lineDirA,
                             const Vector3& linePosB, const Vector3& lineDirB) {
  /// Project the first direction onto the second & renormalize to a unit vector
  const double dirDots = lineDirA.dot(lineDirB);
  const Vector3 aMinusB = linePosA - linePosB;
  if (std::abs(dirDots - 1.) < std::numeric_limits<double>::epsilon()) {
    return (aMinusB - lineDirA.dot(aMinusB) * lineDirA).norm();
  }
  const Vector3 projDir = (lineDirA - dirDots * lineDirB).normalized();
  return aMinusB.cross(lineDirB).dot(projDir);
}
}  // namespace Acts::LineHelper
