// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/PointSurface.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Geometry/GeometryObject.hpp"
#include "Acts/Surfaces/InfiniteBounds.hpp"
#include "Acts/Surfaces/PointBounds.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Surfaces/SurfaceError.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"
#include "Acts/Utilities/detail/OstreamStateGuard.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>

namespace Acts {

namespace {
/// Build the curvilinear measurement frame from the momentum @p direction.
/// Columns are (U, V, T) with T the (normalized) direction, U built from global
/// Z (falling back to global X when the direction is (anti-)parallel to Z).
RotationMatrix3 pointReferenceFrame(const Vector3& direction) {
  const Vector3 T = direction.normalized();
  const bool standard =
      std::abs(T.dot(Vector3::UnitZ())) < s_curvilinearProjTolerance;
  Vector3 U = (standard ? Vector3::UnitZ() : Vector3::UnitX()).cross(T);
  U.normalize();
  const Vector3 V = T.cross(U);

  RotationMatrix3 rframe;
  rframe << U, V, T;
  return rframe;
}
}  // namespace

PointSurface::PointSurface(const Vector3& center)
    : Surface(Transform3(Translation3(center))), m_bounds(nullptr) {}

PointSurface::PointSurface(const Vector3& center, double maxDistance)
    : Surface(Transform3(Translation3(center))),
      m_bounds(std::make_shared<const PointBounds>(maxDistance)) {}

PointSurface::PointSurface(const Transform3& transform,
                           std::shared_ptr<const PointBounds> pbounds)
    : Surface(transform), m_bounds(std::move(pbounds)) {}

PointSurface::PointSurface(const GeometryContext& gctx,
                           const PointSurface& other, const Transform3& shift)
    : GeometryObject{}, Surface(gctx, other, shift), m_bounds(other.m_bounds) {}

Vector3 PointSurface::normal(const GeometryContext& /*gctx*/,
                             const Vector3& /*pos*/,
                             const Vector3& direction) const {
  return direction.normalized();
}

Vector3 PointSurface::referencePosition(const GeometryContext& gctx,
                                        AxisDirection /*aDir*/) const {
  return center(gctx);
}

RotationMatrix3 PointSurface::referenceFrame(const GeometryContext& /*gctx*/,
                                             const Vector3& /*position*/,
                                             const Vector3& direction) const {
  return pointReferenceFrame(direction);
}

Vector3 PointSurface::localToGlobal(const GeometryContext& gctx,
                                    const Vector2& lposition,
                                    const Vector3& direction) const {
  const RotationMatrix3 rframe = pointReferenceFrame(direction);
  return center(gctx) + lposition[0] * rframe.col(0) +
         lposition[1] * rframe.col(1);
}

Result<Vector2> PointSurface::globalToLocal(const GeometryContext& gctx,
                                            const Vector3& position,
                                            const Vector3& direction,
                                            double tolerance) const {
  // Bring the global position into the measurement frame relative to the point.
  const Vector3 localPosition =
      referenceFrame(gctx, position, direction).inverse() *
      (position - center(gctx));

  // `localPosition.z()` is the distance along the direction. It must vanish for
  // `position` to be the point of closest approach (i.e. on the measurement
  // plane through the center).
  if (std::abs(localPosition.z()) > std::abs(tolerance)) {
    return Result<Vector2>::failure(SurfaceError::GlobalPositionNotOnSurface);
  }

  return Result<Vector2>::success(localPosition.head<2>());
}

MultiIntersection3D PointSurface::intersect(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction, const BoundaryTolerance& boundaryTolerance,
    double tolerance) const {
  // The point (center) and the track ray (position, direction)
  const Vector3 pc = center(gctx);

  // Path length along the track to the point of closest approach. Since the
  // measurement plane normal equals the track direction there is always a
  // unique crossing, hence no degeneracy check is needed.
  const double u = (pc - position).dot(direction);

  IntersectionStatus status = std::abs(u) > std::abs(tolerance)
                                  ? IntersectionStatus::reachable
                                  : IntersectionStatus::onSurface;
  const Vector3 result = position + u * direction;

  // Evaluate the boundary check if requested. m_bounds == nullptr means an
  // unbounded point surface, which skips the check.
  if (m_bounds != nullptr && !boundaryTolerance.isInfinite()) {
    const RotationMatrix3 rframe = pointReferenceFrame(direction);
    const Vector3 vecLocal = result - pc;
    const Vector2 local(vecLocal.dot(rframe.col(0)),
                        vecLocal.dot(rframe.col(1)));
    if (!m_bounds->inside(local, boundaryTolerance)) {
      status = IntersectionStatus::unreachable;
    }
  }

  return MultiIntersection3D(Intersection3D(result, u, status));
}

BoundToFreeMatrix PointSurface::boundToFreeJacobian(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction) const {
  assert(isOnSurface(gctx, position, direction, BoundaryTolerance::Infinite()));

  const RotationMatrix3 rframe = referenceFrame(gctx, position, direction);
  const Vector3 U = rframe.col(0);
  const Vector3 V = rframe.col(1);
  const Vector3 T = rframe.col(2);

  // The local (x, y) position on the measurement plane
  const Vector2 local = *globalToLocal(gctx, position, direction,
                                       std::numeric_limits<double>::max());
  const double loc0 = local.x();
  const double loc1 = local.y();

  // Start from the generic jacobian (correct local, direction, time, q/p
  // blocks) and add the position/angle coupling that arises because the
  // measurement frame rotates with the direction.
  BoundToFreeMatrix jacToGlobal =
      Surface::boundToFreeJacobian(gctx, position, direction);

  // For the global-Z frame branch the frame axes are
  //   U = (-sinPhi, cosPhi, 0)
  //   V = (-cosTheta cosPhi, -cosTheta sinPhi, sinTheta)
  // with derivatives
  //   dU/dPhi = cosTheta V - sinTheta T,   dU/dTheta = 0
  //   dV/dPhi = -cosTheta U,               dV/dTheta = T
  // and the free position on surface is P = c + loc0 U + loc1 V, so
  //   dP/dPhi   = loc0 (cosTheta V - sinTheta T) - loc1 cosTheta U
  //   dP/dTheta = loc1 T
  // This generalizes the LineSurface patch, which keeps only the loc0 term
  // (its second axis is fixed to the line). The closed form assumes the
  // standard (Z-based) branch; near direction parallel to Z the frame uses the
  // X-based fall-back and this correction is only approximate (as with the
  // documented forward-eta limitation of the LineSurface).
  const double cosTheta = direction.z();
  const double sinTheta = VectorHelpers::perp(direction);

  jacToGlobal.block<3, 1>(eFreePos0, eBoundPhi) =
      loc0 * (cosTheta * V - sinTheta * T) - loc1 * cosTheta * U;
  jacToGlobal.block<3, 1>(eFreePos0, eBoundTheta) = loc1 * T;

  return jacToGlobal;
}

FreeToPathMatrix PointSurface::freeToPathDerivative(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction) const {
  assert(isOnSurface(gctx, position, direction, BoundaryTolerance::Infinite()));

  // Path to the PCA: u = (c - position) . direction, with the measurement plane
  // normal equal to the direction (so the usual 1/cos term is unity).
  FreeToPathMatrix freeToPath = FreeToPathMatrix::Zero();
  // d u / d position = -direction
  freeToPath.segment<3>(eFreePos0) = -direction.transpose();
  // d u / d direction = (center - position); nonzero at the PCA since the
  // residual is generally nonzero (unlike the curvilinear surface).
  freeToPath.segment<3>(eFreeDir0) = (center(gctx) - position).transpose();

  return freeToPath;
}

AlignmentToPathMatrix PointSurface::alignmentToPathDerivative(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction) const {
  assert(isOnSurface(gctx, position, direction, BoundaryTolerance::Infinite()));
  static_cast<void>(position);

  // u = (c - position) . direction, so d u / d center = direction. A point has
  // no orientation, hence the rotation block is zero.
  AlignmentToPathMatrix alignToPath = AlignmentToPathMatrix::Zero();
  alignToPath.segment<3>(eAlignmentCenter0) = direction.transpose();

  return alignToPath;
}

Matrix<2, 3> PointSurface::localCartesianToBoundLocalDerivative(
    const GeometryContext& /*gctx*/, const Vector3& /*position*/) const {
  // The bound local coordinates are the (x, y) components of the local 3D
  // Cartesian position, so the derivative selects the first two rows. Note the
  // local frame here is the transform frame; the direction-dependent
  // measurement frame is handled in the jacobians above.
  Matrix<2, 3> loc3DToLocBound = Matrix<2, 3>::Zero();
  loc3DToLocBound << 1, 0, 0, 0, 1, 0;
  return loc3DToLocBound;
}

double PointSurface::pathCorrection(const GeometryContext& /*gctx*/,
                                    const Vector3& /*pos*/,
                                    const Vector3& /*mom*/) const {
  return 1.;
}

const SurfaceBounds& PointSurface::bounds() const {
  if (m_bounds != nullptr) {
    return *m_bounds;
  }
  return s_noBounds;
}

const std::shared_ptr<const PointBounds>& PointSurface::boundsPtr() const {
  return m_bounds;
}

void PointSurface::assignSurfaceBounds(
    std::shared_ptr<const PointBounds> newBounds) {
  m_bounds = std::move(newBounds);
}

std::string PointSurface::name() const {
  return "Acts::PointSurface";
}

Polyhedron PointSurface::polyhedronRepresentation(
    const GeometryContext& gctx, unsigned int /*quarterSegments*/) const {
  // Non-physical marker: a small cross centered on the point in the transform
  // frame (a point has no extent, this is for display/debugging only).
  std::vector<Vector3> vertices;
  std::vector<Polyhedron::FaceType> faces;
  std::vector<Polyhedron::FaceType> triangularMesh;

  const Transform3& ctransform = localToGlobalTransform(gctx);
  const double d = 1.;
  vertices.push_back(ctransform * Vector3(-d, 0., 0.));
  vertices.push_back(ctransform * Vector3(d, 0., 0.));
  vertices.push_back(ctransform * Vector3(0., -d, 0.));
  vertices.push_back(ctransform * Vector3(0., d, 0.));
  faces.push_back({0, 1});
  faces.push_back({2, 3});
  triangularMesh.push_back({0, 1, 2});

  return Polyhedron(vertices, faces, triangularMesh);
}

std::ostream& PointSurface::toStreamImpl(const GeometryContext& gctx,
                                         std::ostream& sl) const {
  detail::OstreamStateGuard guard{sl};
  sl << std::fixed << std::setprecision(7);
  sl << "Acts::PointSurface:" << std::endl;
  const Vector3& sfCenter = center(gctx);
  sl << "     Center position  (x, y, z) = (" << sfCenter.x() << ", "
     << sfCenter.y() << ", " << sfCenter.z() << ")";
  return sl;
}

}  // namespace Acts
