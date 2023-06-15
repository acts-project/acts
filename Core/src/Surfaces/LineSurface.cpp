// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/LineSurface.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/detail/TransformationBoundToFree.hpp"
#include "Acts/Geometry/GeometryObject.hpp"
#include "Acts/Surfaces/InfiniteBounds.hpp"
#include "Acts/Surfaces/LineBounds.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Surfaces/SurfaceError.hpp"
#include "Acts/Surfaces/detail/AlignmentHelper.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/ThrowAssert.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace Acts {
class DetectorElementBase;
}  // namespace Acts

Acts::LineSurface::LineSurface(const Transform3& transform, double radius,
                               double halez)
    : GeometryObject(),
      Surface(transform),
      m_bounds(std::make_shared<const LineBounds>(radius, halez)) {}

Acts::LineSurface::LineSurface(const Transform3& transform,
                               std::shared_ptr<const LineBounds> lbounds)
    : GeometryObject(), Surface(transform), m_bounds(std::move(lbounds)) {}

Acts::LineSurface::LineSurface(std::shared_ptr<const LineBounds> lbounds,
                               const DetectorElementBase& detelement)
    : GeometryObject(), Surface(detelement), m_bounds(std::move(lbounds)) {
  throw_assert(m_bounds, "LineBounds must not be nullptr");
}

Acts::LineSurface::LineSurface(const LineSurface& other)
    : GeometryObject(), Surface(other), m_bounds(other.m_bounds) {}

Acts::LineSurface::LineSurface(const GeometryContext& gctx,
                               const LineSurface& other,
                               const Transform3& shift)
    : GeometryObject(), Surface(gctx, other, shift), m_bounds(other.m_bounds) {}

Acts::LineSurface& Acts::LineSurface::operator=(const LineSurface& other) {
  if (this != &other) {
    Surface::operator=(other);
    m_bounds = other.m_bounds;
  }
  return *this;
}

Acts::Vector3 Acts::LineSurface::localToGlobal(const GeometryContext& gctx,
                                               const Vector2& lposition,
                                               const Vector3& momentum) const {
  const auto& sTransform = transform(gctx);
  const auto& tMatrix = sTransform.matrix();
  Vector3 lineDirection(tMatrix(0, 2), tMatrix(1, 2), tMatrix(2, 2));

  // get the vector perpendicular to the momentum and the straw axis
  Vector3 radiusAxisGlobal(lineDirection.cross(momentum));
  Vector3 locZinGlobal = sTransform * Vector3(0., 0., lposition[eBoundLoc1]);
  // add eBoundLoc0 * radiusAxis
  return Vector3(locZinGlobal +
                 lposition[eBoundLoc0] * radiusAxisGlobal.normalized());
}

Acts::Result<Acts::Vector2> Acts::LineSurface::globalToLocal(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& momentum, double tolerance) const {
  using VectorHelpers::perp;
  const auto& sTransform = transform(gctx);
  const auto& tMatrix = sTransform.matrix();
  Vector3 lineDirection(tMatrix(0, 2), tMatrix(1, 2), tMatrix(2, 2));
  // Bring the global position into the local frame
  Vector3 loc3Dframe = sTransform.inverse() * position;
  // construct localPosition with sign*perp(candidate) and z.()
  Vector2 lposition(perp(loc3Dframe), loc3Dframe.z());
  Vector3 sCenter(tMatrix(0, 3), tMatrix(1, 3), tMatrix(2, 3));
  Vector3 decVec(position - sCenter);
  // assign the right sign
  double sign = ((lineDirection.cross(momentum)).dot(decVec) < 0.) ? -1. : 1.;
  lposition[eBoundLoc0] *= sign;

  if ((localToGlobal(gctx, lposition, momentum) - position).norm() >
      tolerance) {
    return Result<Vector2>::failure(SurfaceError::GlobalPositionNotOnSurface);
  }

  return Result<Vector2>::success(lposition);
}

std::string Acts::LineSurface::name() const {
  return "Acts::LineSurface";
}

Acts::RotationMatrix3 Acts::LineSurface::referenceFrame(
    const GeometryContext& gctx, const Vector3& /*position*/,
    const Vector3& momentum) const {
  RotationMatrix3 mFrame;
  const auto& tMatrix = transform(gctx).matrix();
  Vector3 measY(tMatrix(0, 2), tMatrix(1, 2), tMatrix(2, 2));
  Vector3 measX(measY.cross(momentum).normalized());
  Vector3 measDepth(measX.cross(measY));
  // assign the columnes
  mFrame.col(0) = measX;
  mFrame.col(1) = measY;
  mFrame.col(2) = measDepth;
  // return the rotation matrix
  return mFrame;
}

double Acts::LineSurface::pathCorrection(const GeometryContext& /*gctx*/,
                                         const Vector3& /*pos*/,
                                         const Vector3& /*mom*/) const {
  return 1.;
}

Acts::Vector3 Acts::LineSurface::binningPosition(
    const GeometryContext& gctx, BinningValue /*bValue*/) const {
  return center(gctx);
}

Acts::Vector3 Acts::LineSurface::normal(const GeometryContext& gctx,
                                        const Vector2& /*lpos*/) const {
  const auto& tMatrix = transform(gctx).matrix();
  return Vector3(tMatrix(0, 2), tMatrix(1, 2), tMatrix(2, 2));
}

const Acts::SurfaceBounds& Acts::LineSurface::bounds() const {
  if (m_bounds) {
    return (*m_bounds.get());
  }
  return s_noBounds;
}

Acts::SurfaceIntersection Acts::LineSurface::intersect(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction, const BoundaryCheck& bcheck,
    ActsScalar tolerance) const {
  // following nominclature found in header file and doxygen documentation
  // line one is the straight track
  const Vector3& ma = position;
  const Vector3& ea = direction;
  // line two is the line surface
  const auto& tMatrix = transform(gctx).matrix();
  Vector3 mb = tMatrix.block<3, 1>(0, 3).transpose();
  Vector3 eb = tMatrix.block<3, 1>(0, 2).transpose();
  // now go ahead and solve for the closest approach
  Vector3 mab(mb - ma);
  double eaTeb = ea.dot(eb);
  double denom = 1 - eaTeb * eaTeb;
  // validity parameter
  Intersection3D::Status status = Intersection3D::Status::unreachable;
  if (std::abs(denom) > std::abs(tolerance)) {
    double u = (mab.dot(ea) - mab.dot(eb) * eaTeb) / denom;
    // Check if we are on the surface already
    status = std::abs(u) < std::abs(tolerance)
                 ? Intersection3D::Status::onSurface
                 : Intersection3D::Status::reachable;
    Vector3 result = (ma + u * ea);
    // Evaluate the boundary check if requested
    // m_bounds == nullptr prevents unecessary calulations for PerigeeSurface
    if (bcheck and m_bounds) {
      // At closest approach: check inside R or and inside Z
      const Vector3 vecLocal(result - mb);
      double cZ = vecLocal.dot(eb);
      double hZ = m_bounds->get(LineBounds::eHalfLengthZ) + tolerance;
      if ((std::abs(cZ) > std::abs(hZ)) or
          ((vecLocal - cZ * eb).norm() >
           m_bounds->get(LineBounds::eR) + tolerance)) {
        status = Intersection3D::Status::missed;
      }
    }
    return {Intersection3D(result, u, status), this};
  }
  // return a false intersection
  return {Intersection3D(position, std::numeric_limits<double>::max(), status),
          this};
}

Acts::BoundToFreeMatrix Acts::LineSurface::boundToFreeJacobian(
    const GeometryContext& gctx, const BoundVector& boundParams) const {
  // Transform from bound to free parameters
  FreeVector freeParams =
      detail::transformBoundToFreeParameters(*this, gctx, boundParams);
  // The global position
  const Vector3 position = freeParams.segment<3>(eFreePos0);
  // The direction
  const Vector3 direction = freeParams.segment<3>(eFreeDir0);
  // Get the sines and cosines directly
  const double cos_theta = std::cos(boundParams[eBoundTheta]);
  const double sin_theta = std::sin(boundParams[eBoundTheta]);
  const double cos_phi = std::cos(boundParams[eBoundPhi]);
  const double sin_phi = std::sin(boundParams[eBoundPhi]);
  // retrieve the reference frame
  const auto rframe = referenceFrame(gctx, position, direction);
  // Initialize the jacobian from local to global
  BoundToFreeMatrix jacToGlobal = BoundToFreeMatrix::Zero();
  // the local error components - given by the reference frame
  jacToGlobal.topLeftCorner<3, 2>() = rframe.topLeftCorner<3, 2>();
  // the time component
  jacToGlobal(eFreeTime, eBoundTime) = 1;
  // the momentum components
  jacToGlobal(eFreeDir0, eBoundPhi) = (-sin_theta) * sin_phi;
  jacToGlobal(eFreeDir0, eBoundTheta) = cos_theta * cos_phi;
  jacToGlobal(eFreeDir1, eBoundPhi) = sin_theta * cos_phi;
  jacToGlobal(eFreeDir1, eBoundTheta) = cos_theta * sin_phi;
  jacToGlobal(eFreeDir2, eBoundTheta) = (-sin_theta);
  jacToGlobal(eFreeQOverP, eBoundQOverP) = 1;

  // the projection of direction onto ref frame normal
  double ipdn = 1. / direction.dot(rframe.col(2));
  // build the cross product of d(D)/d(eBoundPhi) components with y axis
  auto dDPhiY = rframe.block<3, 1>(0, 1).cross(
      jacToGlobal.block<3, 1>(eFreeDir0, eBoundPhi));
  // and the same for the d(D)/d(eTheta) components
  auto dDThetaY = rframe.block<3, 1>(0, 1).cross(
      jacToGlobal.block<3, 1>(eFreeDir0, eBoundTheta));
  // and correct for the x axis components
  dDPhiY -= rframe.block<3, 1>(0, 0) * (rframe.block<3, 1>(0, 0).dot(dDPhiY));
  dDThetaY -=
      rframe.block<3, 1>(0, 0) * (rframe.block<3, 1>(0, 0).dot(dDThetaY));
  // set the jacobian components for global d/ phi/Theta
  jacToGlobal.block<3, 1>(eFreePos0, eBoundPhi) =
      dDPhiY * boundParams[eBoundLoc0] * ipdn;
  jacToGlobal.block<3, 1>(eFreePos0, eBoundTheta) =
      dDThetaY * boundParams[eBoundLoc0] * ipdn;
  return jacToGlobal;
}

Acts::FreeToPathMatrix Acts::LineSurface::freeToPathDerivative(
    const GeometryContext& gctx, const FreeVector& parameters) const {
  // The global posiiton
  const auto position = parameters.segment<3>(eFreePos0);
  // The direction
  const auto direction = parameters.segment<3>(eFreeDir0);
  // The vector between position and center
  const auto pcRowVec = (position - center(gctx)).transpose().eval();
  // The rotation
  const auto& rotation = transform(gctx).rotation();
  // The local frame z axis
  const auto& localZAxis = rotation.col(2);
  // The local z coordinate
  const auto pz = pcRowVec * localZAxis;
  // Cosine of angle between momentum direction and local frame z axis
  const auto dz = localZAxis.dot(direction);
  const auto norm = 1 / (1 - dz * dz);
  // Initialize the derivative of propagation path w.r.t. free parameter
  FreeToPathMatrix freeToPath = FreeToPathMatrix::Zero();
  // The derivative of path w.r.t. position
  freeToPath.segment<3>(eFreePos0) =
      norm * (dz * localZAxis.transpose() - direction.transpose());
  // The derivative of path w.r.t. direction
  freeToPath.segment<3>(eFreeDir0) =
      norm * (pz * localZAxis.transpose() - pcRowVec);

  return freeToPath;
}

Acts::AlignmentToPathMatrix Acts::LineSurface::alignmentToPathDerivative(
    const GeometryContext& gctx, const FreeVector& parameters) const {
  // The global posiiton
  const auto position = parameters.segment<3>(eFreePos0);
  // The direction
  const auto direction = parameters.segment<3>(eFreeDir0);
  // The vector between position and center
  const auto pcRowVec = (position - center(gctx)).transpose().eval();
  // The rotation
  const auto& rotation = transform(gctx).rotation();
  // The local frame z axis
  const Vector3 localZAxis = rotation.col(2);
  // The local z coordinate
  const auto pz = pcRowVec * localZAxis;
  // Cosine of angle between momentum direction and local frame z axis
  const auto dz = localZAxis.dot(direction);
  const auto norm = 1 / (1 - dz * dz);
  // Calculate the derivative of local frame axes w.r.t its rotation
  const auto [rotToLocalXAxis, rotToLocalYAxis, rotToLocalZAxis] =
      detail::rotationToLocalAxesDerivative(rotation);
  // Initialize the derivative of propagation path w.r.t. local frame
  // translation (origin) and rotation
  AlignmentToPathMatrix alignToPath = AlignmentToPathMatrix::Zero();
  alignToPath.segment<3>(eAlignmentCenter0) =
      norm * (direction.transpose() - dz * localZAxis.transpose());
  alignToPath.segment<3>(eAlignmentRotation0) =
      norm * (dz * pcRowVec + pz * direction.transpose()) * rotToLocalZAxis;

  return alignToPath;
}

Acts::ActsMatrix<2, 3> Acts::LineSurface::localCartesianToBoundLocalDerivative(
    const GeometryContext& gctx, const Vector3& position) const {
  using VectorHelpers::phi;
  // The local frame transform
  const auto& sTransform = transform(gctx);
  // calculate the transformation to local coorinates
  const Vector3 localPos = sTransform.inverse() * position;
  const double lphi = phi(localPos);
  const double lcphi = std::cos(lphi);
  const double lsphi = std::sin(lphi);
  ActsMatrix<2, 3> loc3DToLocBound = ActsMatrix<2, 3>::Zero();
  loc3DToLocBound << lcphi, lsphi, 0, 0, 0, 1;

  return loc3DToLocBound;
}