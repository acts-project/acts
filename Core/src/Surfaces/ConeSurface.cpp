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

namespace {
/// Local azimuth and the bound axes at a cone position, in the global frame.
/// radial = (cos phi, sin phi, 0), phiAxis = (-sin phi, cos phi, 0) and
/// axis = ez. referenceFrame is (phiAxis, n x phiAxis, n). The bound z
/// derivative tan(alpha)*radial - phi*tan(alpha)*phiAxis + axis is not that
/// second column, so boundToFreeJacobian replaces it and freeToBoundJacobian
/// replaces both loc rows.
struct ConeBoundAxes {
  double localPhi;
  Vector3 radial;
  Vector3 phiAxis;
  Vector3 axis;
};

ConeBoundAxes coneBoundAxes(const Transform3& toGlobal,
                            const Vector3& position) {
  const Vector3 localPosition = toGlobal.inverse() * position;
  const double localPhi = phi(localPosition);
  const double cosPhi = std::cos(localPhi);
  const double sinPhi = std::sin(localPhi);
  const auto& rotation = toGlobal.rotation();
  return {localPhi, rotation * Vector3(cosPhi, sinPhi, 0.),
          rotation * Vector3(-sinPhi, cosPhi, 0.), rotation.col(2)};
}
}  // namespace

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
  // Columns are (e_phi, n x e_phi, n). The normal is what the base
  // freeToPathDerivative projects onto.
  const Vector3 n = normal(gctx, position);
  const Transform3& toGlobal = localToGlobalTransform(gctx);
  const double localPhi = phi(toGlobal.inverse() * position);
  const Vector3 phiAxis = toGlobal.rotation() *
                          Vector3(-std::sin(localPhi), std::cos(localPhi), 0.);
  RotationMatrix3 frame;
  frame.col(0) = phiAxis;
  frame.col(1) = n.cross(phiAxis);
  frame.col(2) = n;
  return frame;
}

BoundToFreeMatrix ConeSurface::boundToFreeJacobian(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction) const {
  assert(isOnSurface(gctx, position, direction, BoundaryTolerance::Infinite()));

  // referenceFrame supplies e_phi as the loc0 column. Replace loc1.
  BoundToFreeMatrix jacToGlobal =
      Surface::boundToFreeJacobian(gctx, position, direction);

  const ConeBoundAxes axes =
      coneBoundAxes(localToGlobalTransform(gctx), position);
  const double tanAlpha = bounds().tanAlpha();
  jacToGlobal.block<3, 1>(eFreePos0, eBoundLoc1) =
      tanAlpha * axes.radial - axes.localPhi * tanAlpha * axes.phiAxis +
      axes.axis;
  return jacToGlobal;
}

FreeToBoundMatrix ConeSurface::freeToBoundJacobian(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction) const {
  assert(isOnSurface(gctx, position, direction, BoundaryTolerance::Infinite()));

  FreeToBoundMatrix jacToLocal =
      Surface::freeToBoundJacobian(gctx, position, direction);

  // On the surface R/r = 1, so the (rphi, z) rows are e_phi + phi*tan(alpha)*ez
  // and ez.
  const ConeBoundAxes axes =
      coneBoundAxes(localToGlobalTransform(gctx), position);
  const double tanAlpha = bounds().tanAlpha();
  jacToLocal.block<1, 3>(eBoundLoc0, eFreePos0) =
      (axes.phiAxis + axes.localPhi * tanAlpha * axes.axis).transpose();
  jacToLocal.block<1, 3>(eBoundLoc1, eFreePos0) = axes.axis.transpose();
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
  return std::abs(1. / normal(gctx, position).dot(direction));
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
  // (cos phi cos alpha, sin phi cos alpha, -sign(z) sin alpha)
  const Transform3& toGlobal = localToGlobalTransform(gctx);
  const Vector3 localPosition = toGlobal.inverse() * position;
  const double localPhi = phi(localPosition);
  const double cosAlpha = std::cos(bounds().get(ConeBounds::eAlpha));
  const double sinAlpha = std::sin(bounds().get(ConeBounds::eAlpha));
  const double sgn = -std::copysign(1., localPosition.z());
  const Vector3 localNormal(std::cos(localPhi) * cosAlpha,
                            std::sin(localPhi) * cosAlpha, sgn * sinAlpha);
  return toGlobal.linear() * localNormal;
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
  Vector3 solution2 = position + qe.second * direction;
  IntersectionStatus status2 = std::abs(qe.second) < std::abs(tolerance)
                                   ? IntersectionStatus::onSurface
                                   : IntersectionStatus::reachable;
  if (!boundaryTolerance.isInfinite() &&
      !isOnSurface(gctx, solution2, direction, boundaryTolerance)) {
    status2 = IntersectionStatus::unreachable;
  }

  // solution1/2 are already in the global frame.
  Intersection3D first(solution1, qe.first, status1);
  Intersection3D second(solution2, qe.second, status2);
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

  // Local-frame constraint F(x, y, z) = x^2 + y^2 - z^2 tan^2(alpha) = 0.
  // Implicit differentiation at the current position gives the
  // alignment-to-path derivative.
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
  // A rotation w about the local axes moves the surface point at localPos by
  // cross(w, localPos), so ds/dw = cross(localPos, grad(F) / 2) / denom
  // = (1 + tan^2(alpha)) z (-y, x, 0) / denom.
  const double rotationScale = (1. + tanAlpha2) * localPos.z() / denom;
  alignToPath.segment<3>(eAlignmentRotation0) =
      Vector3(-rotationScale * localPos.y(), rotationScale * localPos.x(), 0.)
          .transpose();

  return alignToPath;
}

Matrix<2, 3> ConeSurface::localCartesianToBoundLocalDerivative(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& /*direction*/) const {
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
