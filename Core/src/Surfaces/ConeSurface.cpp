// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/ConeSurface.hpp"

#include "Acts/Geometry/GeometryObject.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/SurfaceError.hpp"
#include "Acts/Surfaces/detail/AlignmentHelper.hpp"
#include "Acts/Surfaces/detail/FacesHelper.hpp"
#include "Acts/Surfaces/detail/VerticesHelper.hpp"
#include "Acts/Utilities/AlgebraHelpers.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/ThrowAssert.hpp"
#include "Acts/Utilities/detail/RealQuadraticEquation.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Acts {

using VectorHelpers::perp;
using VectorHelpers::phi;

ConeSurface::ConeSurface(const ConeSurface& other)
    : GeometryObject{}, RegularSurface(other), m_bounds(other.m_bounds) {}

ConeSurface::ConeSurface(const GeometryContext& gctx, const ConeSurface& other,
                         const Transform3& shift)
    : RegularSurface(gctx, other, shift), m_bounds(other.m_bounds) {}

ConeSurface::ConeSurface(const Transform3& transform, double alpha,
                         bool symmetric)
    : RegularSurface(transform),
      m_bounds(std::make_shared<const ConeBounds>(alpha, symmetric)) {}

ConeSurface::ConeSurface(const Transform3& transform, double alpha, double zmin,
                         double zmax, double halfPhi)
    : RegularSurface(transform),
      m_bounds(std::make_shared<const ConeBounds>(alpha, zmin, zmax, halfPhi)) {
}

ConeSurface::ConeSurface(const Transform3& transform,
                         std::shared_ptr<const ConeBounds> cbounds)
    : RegularSurface(transform), m_bounds(std::move(cbounds)) {
  throw_assert(m_bounds, "ConeBounds must not be nullptr");
}

Vector3 ConeSurface::referencePosition(const GeometryContext& gctx,
                                       AxisDirection aDir) const {
  const Vector3& sfCenter = center(gctx);

  // special binning type for R-type methods
  if (aDir == AxisDirection::AxisR || aDir == AxisDirection::AxisRPhi) {
    return Vector3(sfCenter.x() + bounds().r(sfCenter.z()), sfCenter.y(),
                   sfCenter.z());
  }
  // give the center as default for all of these binning types
  // AxisDirection::AxisX, AxisDirection::AxisY, AxisDirection::AxisZ,
  // AxisDirection::AxisR, AxisDirection::AxisPhi, AxisDirection::AxisRPhi,
  // AxisDirection::AxisTheta, AxisDirection::AxisEta
  return sfCenter;
}

Surface::SurfaceType ConeSurface::type() const {
  return Surface::Cone;
}

ConeSurface& ConeSurface::operator=(const ConeSurface& other) {
  if (this != &other) {
    Surface::operator=(other);
    m_bounds = other.m_bounds;
  }
  return *this;
}

Vector3 ConeSurface::rotSymmetryAxis(const GeometryContext& gctx) const {
  return localToGlobalTransform(gctx).matrix().block<3, 1>(0, 2);
}

RotationMatrix3 ConeSurface::referenceFrame(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& /*direction*/) const {
  RotationMatrix3 mFrame;
  // construct the measurement frame
  // measured Y is the local z axis
  Vector3 measY = rotSymmetryAxis(gctx);
  // measured z is the position transverse normalized
  Vector3 measDepth = Vector3(position.x(), position.y(), 0.).normalized();
  // measured X is what comoes out of it
  Vector3 measX(measY.cross(measDepth).normalized());
  // the columnes
  mFrame.col(0) = measX;
  mFrame.col(1) = measY;
  mFrame.col(2) = measDepth;
  // return the rotation matrix
  //!< @todo fold in alpha
  // return it
  return mFrame;
}

FreeToBoundMatrix ConeSurface::freeToBoundJacobian(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction) const {
  // t, φ, θ and q/p come from the base. Its position block follows
  // referenceFrame, whose transverse axis is the global xy direction and does
  // not stay along eφ once the cone is offset. Bound (rφ, z) is obtained by
  // applying (δu, δz) -> ( (R/r) δu + φ tanα δz, δz ) to (eφ, ez).
  FreeToBoundMatrix jacToLocal =
      Surface::freeToBoundJacobian(gctx, position, direction);

  const Transform3& toGlobal = localToGlobalTransform(gctx);
  const Vector3 localPosition = toGlobal.inverse() * position;
  const double localPhi = phi(localPosition);
  const auto& rotation = toGlobal.rotation();
  const Vector3 phiAxis = rotation.col(2).cross(
      rotation * Vector3(std::cos(localPhi), std::sin(localPhi), 0.));
  const double radiusOverPerp =
      localPosition.z() * bounds().tanAlpha() / perp(localPosition);
  Matrix<2, 2> toBound;
  toBound << radiusOverPerp, localPhi * bounds().tanAlpha(), 0., 1.;
  Matrix<2, 3> axisFrame;
  axisFrame.row(0) = phiAxis.transpose();
  axisFrame.row(1) = rotation.col(2).transpose();
  jacToLocal.block<2, 3>(eBoundLoc0, eFreePos0) = toBound * axisFrame;
  return jacToLocal;
}

Vector3 ConeSurface::localToGlobal(const GeometryContext& gctx,
                                   const Vector2& lposition) const {
  // create the position in the local 3d frame
  double r = lposition[1] * bounds().tanAlpha();
  double phi = lposition[0] / r;
  Vector3 loc3Dframe(r * std::cos(phi), r * std::sin(phi), lposition[1]);
  return localToGlobalTransform(gctx) * loc3Dframe;
}

Result<Vector2> ConeSurface::globalToLocal(const GeometryContext& gctx,
                                           const Vector3& position,
                                           double tolerance) const {
  Vector3 loc3Dframe = localToGlobalTransform(gctx).inverse() * position;
  double r = loc3Dframe.z() * bounds().tanAlpha();
  if (std::abs(perp(loc3Dframe) - r) > tolerance) {
    return Result<Vector2>::failure(SurfaceError::GlobalPositionNotOnSurface);
  }
  return Result<Vector2>::success(
      Vector2(r * std::atan2(loc3Dframe.y(), loc3Dframe.x()), loc3Dframe.z()));
}

double ConeSurface::pathCorrection(const GeometryContext& gctx,
                                   const Vector3& position,
                                   const Vector3& direction) const {
  // (cos phi cos alpha, sin phi cos alpha, sgn z sin alpha)
  Vector3 posLocal = localToGlobalTransform(gctx).inverse() * position;
  double phi = VectorHelpers::phi(posLocal);
  double sgn = -std::copysign(1., posLocal.z());
  double cosAlpha = std::cos(bounds().get(ConeBounds::eAlpha));
  double sinAlpha = std::sin(bounds().get(ConeBounds::eAlpha));
  Vector3 normalC(std::cos(phi) * cosAlpha, std::sin(phi) * cosAlpha,
                  sgn * sinAlpha);
  normalC = localToGlobalTransform(gctx).linear() * normalC;
  // Back to the global frame
  double cAlpha = normalC.dot(direction);
  return std::abs(1. / cAlpha);
}

std::string ConeSurface::name() const {
  return "Acts::ConeSurface";
}

Vector3 ConeSurface::normal(const GeometryContext& gctx,
                            const Vector2& lposition) const {
  // (cos phi cos alpha, sin phi cos alpha, sgn z sin alpha)
  double phi = lposition[0] / (bounds().r(lposition[1])),
         sgn = -std::copysign(1., lposition[1]);
  double cosAlpha = std::cos(bounds().get(ConeBounds::eAlpha));
  double sinAlpha = std::sin(bounds().get(ConeBounds::eAlpha));
  Vector3 localNormal(std::cos(phi) * cosAlpha, std::sin(phi) * cosAlpha,
                      sgn * sinAlpha);
  return Vector3(localToGlobalTransform(gctx).linear() * localNormal);
}

Vector3 ConeSurface::normal(const GeometryContext& gctx,
                            const Vector3& position) const {
  // get it into the cylinder frame if needed
  // @todo respect opening angle
  Vector3 pos3D = localToGlobalTransform(gctx).inverse() * position;
  pos3D.z() = 0;
  return pos3D.normalized();
}

const ConeBounds& ConeSurface::bounds() const {
  // is safe because no constructor w/o bounds exists
  return *m_bounds;
}

Polyhedron ConeSurface::polyhedronRepresentation(
    const GeometryContext& gctx, unsigned int quarterSegments) const {
  // Prepare vertices and faces
  std::vector<Vector3> vertices;
  std::vector<Polyhedron::FaceType> faces;
  std::vector<Polyhedron::FaceType> triangularMesh;
  double minZ = bounds().get(ConeBounds::eMinZ);
  double maxZ = bounds().get(ConeBounds::eMaxZ);

  if (minZ == -std::numeric_limits<double>::infinity() ||
      maxZ == std::numeric_limits<double>::infinity()) {
    throw std::domain_error(
        "Polyhedron representation of boundless surface is not possible");
  }

  auto ctransform = localToGlobalTransform(gctx);

  // The tip - created only once and only, if it is not a cut-off cone
  bool tipExists = false;
  if (minZ * maxZ <= s_onSurfaceTolerance) {
    vertices.push_back(ctransform * Vector3(0., 0., 0.));
    tipExists = true;
  }

  // Cone parameters
  double hPhiSec = bounds().get(ConeBounds::eHalfPhiSector);
  double avgPhi = bounds().get(ConeBounds::eAveragePhi);
  std::vector<double> refPhi = {};
  if (bool fullCone =
          std::abs(hPhiSec - std::numbers::pi) < s_fullAzimuthTolerance;
      !fullCone) {
    refPhi = {avgPhi};
  }

  // Add the cone sizes
  std::vector<double> coneSides;
  if (std::abs(minZ) > s_onSurfaceTolerance) {
    coneSides.push_back(minZ);
  }
  if (std::abs(maxZ) > s_onSurfaceTolerance) {
    coneSides.push_back(maxZ);
  }

  for (auto& z : coneSides) {
    std::size_t firstIv = vertices.size();
    // Radius and z offset
    double r = std::abs(z) * bounds().tanAlpha();
    Vector3 zoffset(0., 0., z);
    auto svertices = detail::VerticesHelper::segmentVertices(
        {r, r}, avgPhi - hPhiSec, avgPhi + hPhiSec, refPhi, quarterSegments,
        zoffset, ctransform);
    vertices.insert(vertices.end(), svertices.begin(), svertices.end());
    // If the tip exists, the faces need to be triangular
    if (tipExists) {
      for (std::size_t iv = firstIv + 1; iv < svertices.size() + firstIv;
           ++iv) {
        std::size_t one = 0, two = iv, three = iv - 1;
        if (z < 0.) {
          std::swap(two, three);
        }
        faces.push_back({one, two, three});
      }
    }
  }

  // if no tip exists, connect the two bows
  if (tipExists) {
    triangularMesh = faces;
  } else {
    auto facesMesh = detail::FacesHelper::cylindricalFaceMesh(vertices);
    faces = facesMesh.first;
    triangularMesh = facesMesh.second;
  }

  return Polyhedron(vertices, faces, triangularMesh, false);
}

detail::RealQuadraticEquation ConeSurface::intersectionSolver(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction) const {
  // Transform into the local frame
  Transform3 invTrans = localToGlobalTransform(gctx).inverse();
  Vector3 point1 = invTrans * position;
  Vector3 dir1 = invTrans.linear() * direction;

  // See file header for the formula derivation
  double tan2Alpha = bounds().tanAlpha() * bounds().tanAlpha(),
         A = dir1.x() * dir1.x() + dir1.y() * dir1.y() -
             tan2Alpha * dir1.z() * dir1.z(),
         B = 2 * (dir1.x() * point1.x() + dir1.y() * point1.y() -
                  tan2Alpha * dir1.z() * point1.z()),
         C = point1.x() * point1.x() + point1.y() * point1.y() -
             tan2Alpha * point1.z() * point1.z();
  if (A == 0.) {
    A += 1e-16;  // avoid division by zero
  }

  return detail::RealQuadraticEquation(A, B, C);
}

MultiIntersection3D ConeSurface::intersect(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction, const BoundaryTolerance& boundaryTolerance,
    double tolerance) const {
  // Solve the quadratic equation
  auto qe = intersectionSolver(gctx, position, direction);

  // If no valid solution return a non-valid surfaceIntersection
  if (qe.solutions == 0) {
    return MultiIntersection3D(Intersection3D::Invalid(),
                               Intersection3D::Invalid());
  }

  // Check the validity of the first solution
  Vector3 solution1 = position + qe.first * direction;
  IntersectionStatus status1 = std::abs(qe.first) < std::abs(tolerance)
                                   ? IntersectionStatus::onSurface
                                   : IntersectionStatus::reachable;

  if (!boundaryTolerance.isInfinite() &&
      !isOnSurface(gctx, solution1, direction, boundaryTolerance)) {
    status1 = IntersectionStatus::unreachable;
  }

  // Check the validity of the second solution
  Vector3 solution2 = position + qe.first * direction;
  IntersectionStatus status2 = std::abs(qe.second) < std::abs(tolerance)
                                   ? IntersectionStatus::onSurface
                                   : IntersectionStatus::reachable;
  if (!boundaryTolerance.isInfinite() &&
      !isOnSurface(gctx, solution2, direction, boundaryTolerance)) {
    status2 = IntersectionStatus::unreachable;
  }

  const auto& tf = localToGlobalTransform(gctx);
  // Set the intersection
  Intersection3D first(tf * solution1, qe.first, status1);
  Intersection3D second(tf * solution2, qe.second, status2);
  // Order based on path length
  if (first.pathLength() <= second.pathLength()) {
    return MultiIntersection3D(first, second);
  }
  return MultiIntersection3D(second, first);
}

AlignmentToPathMatrix ConeSurface::alignmentToPathDerivative(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction) const {
  assert(isOnSurface(gctx, position, direction, BoundaryTolerance::Infinite()));

  // Constraint F(x, y, z) = x^2 + y^2 - z^2 tan^2(alpha) = 0.
  // Implicit differentiation gives the alignment-to-path derivative at s = 0.
  const auto& rotation = localToGlobalTransform(gctx).rotation();
  const Vector3 localPos = rotation.transpose() * (position - center(gctx));
  const Vector3 localDir = rotation.transpose() * direction;
  const double tanAlpha2 = bounds().tanAlpha() * bounds().tanAlpha();
  // Vanishes when the track is tangent to the cone.
  const double denom = localPos.x() * localDir.x() +
                       localPos.y() * localDir.y() -
                       tanAlpha2 * localPos.z() * localDir.z();

  AlignmentToPathMatrix alignToPath = AlignmentToPathMatrix::Zero();
  const Vector3 globalCenterDerivative =
      (localPos.x() * rotation.col(0) + localPos.y() * rotation.col(1) -
       tanAlpha2 * localPos.z() * rotation.col(2)) /
      denom;
  detail::setAlignToPathLocalCenterDerivative(alignToPath,
                                              globalCenterDerivative, rotation);
  // A rotation about local x changes the axis by -local y, about local y by
  // +local x, and about local z leaves the axis unchanged.
  const double rotationScale = (1. + tanAlpha2) * localPos.z() / denom;
  alignToPath.segment<3>(eAlignmentRotation0) =
      Vector3(-rotationScale * localPos.y(), rotationScale * localPos.x(), 0.)
          .transpose();

  return alignToPath;
}

Matrix<2, 3> ConeSurface::localCartesianToBoundLocalDerivative(
    const GeometryContext& gctx, const Vector3& position) const {
  using VectorHelpers::perp;
  using VectorHelpers::phi;
  // The local frame transform
  const auto& sTransform = localToGlobalTransform(gctx);
  // calculate the transformation to local coordinates
  const Vector3 localPos = sTransform.inverse() * position;
  const double lr = perp(localPos);
  const double lphi = phi(localPos);
  const double lcphi = std::cos(lphi);
  const double lsphi = std::sin(lphi);
  // Solve for radius R
  const double R = localPos.z() * bounds().tanAlpha();
  Matrix<2, 3> loc3DToLocBound = Matrix<2, 3>::Zero();
  loc3DToLocBound << -R * lsphi / lr, R * lcphi / lr,
      lphi * bounds().tanAlpha(), 0, 0, 1;

  return loc3DToLocBound;
}

const std::shared_ptr<const ConeBounds>& ConeSurface::boundsPtr() const {
  return m_bounds;
}

void ConeSurface::assignSurfaceBounds(
    std::shared_ptr<const ConeBounds> newBounds) {
  m_bounds = std::move(newBounds);
}

}  // namespace Acts
