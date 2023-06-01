// This file is part of the Acts project.
//
// Copyright (C) 2016-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/DiscSurface.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/detail/TransformationBoundToFree.hpp"
#include "Acts/Surfaces/DiscTrapezoidBounds.hpp"
#include "Acts/Surfaces/InfiniteBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/SurfaceError.hpp"
#include "Acts/Surfaces/detail/FacesHelper.hpp"
#include "Acts/Surfaces/detail/PlanarHelper.hpp"
#include "Acts/Utilities/ThrowAssert.hpp"

#include <cmath>
#include <system_error>
#include <vector>

using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;

Acts::DiscSurface::DiscSurface(const DiscSurface& other)
    : GeometryObject(), Surface(other), m_bounds(other.m_bounds) {}

Acts::DiscSurface::DiscSurface(const GeometryContext& gctx,
                               const DiscSurface& other,
                               const Transform3& shift)
    : GeometryObject(), Surface(gctx, other, shift), m_bounds(other.m_bounds) {}

Acts::DiscSurface::DiscSurface(const Transform3& transform, double rmin,
                               double rmax, double hphisec)
    : GeometryObject(),
      Surface(transform),
      m_bounds(std::make_shared<const RadialBounds>(rmin, rmax, hphisec)) {}

Acts::DiscSurface::DiscSurface(const Transform3& transform, double minhalfx,
                               double maxhalfx, double minR, double maxR,
                               double avephi, double stereo)
    : GeometryObject(),
      Surface(transform),
      m_bounds(std::make_shared<const DiscTrapezoidBounds>(
          minhalfx, maxhalfx, minR, maxR, avephi, stereo)) {}

Acts::DiscSurface::DiscSurface(const Transform3& transform,
                               std::shared_ptr<const DiscBounds> dbounds)
    : GeometryObject(), Surface(transform), m_bounds(std::move(dbounds)) {}

Acts::DiscSurface::DiscSurface(const std::shared_ptr<const DiscBounds>& dbounds,
                               const DetectorElementBase& detelement)
    : GeometryObject(), Surface(detelement), m_bounds(dbounds) {
  throw_assert(dbounds, "nullptr as DiscBounds");
}

Acts::DiscSurface& Acts::DiscSurface::operator=(const DiscSurface& other) {
  if (this != &other) {
    Acts::Surface::operator=(other);
    m_bounds = other.m_bounds;
  }
  return *this;
}

Acts::Surface::SurfaceType Acts::DiscSurface::type() const {
  return Surface::Disc;
}

Acts::Vector3 Acts::DiscSurface::localToGlobal(const GeometryContext& gctx,
                                               const Vector2& lposition,
                                               const Vector3& /*gmom*/) const {
  // create the position in the local 3d frame
  Vector3 loc3Dframe(
      lposition[Acts::eBoundLoc0] * cos(lposition[Acts::eBoundLoc1]),
      lposition[Acts::eBoundLoc0] * sin(lposition[Acts::eBoundLoc1]), 0.);
  // transform to globalframe
  return transform(gctx) * loc3Dframe;
}

Acts::Result<Acts::Vector2> Acts::DiscSurface::globalToLocal(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& /*gmom*/, double tolerance) const {
  // transport it to the globalframe
  Vector3 loc3Dframe = (transform(gctx).inverse()) * position;
  if (std::abs(loc3Dframe.z()) > std::abs(tolerance)) {
    return Result<Vector2>::failure(SurfaceError::GlobalPositionNotOnSurface);
  }
  return Result<Acts::Vector2>::success({perp(loc3Dframe), phi(loc3Dframe)});
}

Acts::Vector2 Acts::DiscSurface::localPolarToLocalCartesian(
    const Vector2& locpol) const {
  const DiscTrapezoidBounds* dtbo =
      dynamic_cast<const Acts::DiscTrapezoidBounds*>(&(bounds()));
  if (dtbo != nullptr) {
    double rMedium = dtbo->rCenter();
    double phi = dtbo->get(DiscTrapezoidBounds::eAveragePhi);

    Vector2 polarCenter(rMedium, phi);
    Vector2 cartCenter = localPolarToCartesian(polarCenter);
    Vector2 cartPos = localPolarToCartesian(locpol);
    Vector2 Pos = cartPos - cartCenter;

    Acts::Vector2 locPos(
        Pos[Acts::eBoundLoc0] * sin(phi) - Pos[Acts::eBoundLoc1] * cos(phi),
        Pos[Acts::eBoundLoc1] * sin(phi) + Pos[Acts::eBoundLoc0] * cos(phi));
    return Vector2(locPos[Acts::eBoundLoc0], locPos[Acts::eBoundLoc1]);
  }
  return Vector2(locpol[Acts::eBoundLoc0] * cos(locpol[Acts::eBoundLoc1]),
                 locpol[Acts::eBoundLoc0] * sin(locpol[Acts::eBoundLoc1]));
}

Acts::Vector3 Acts::DiscSurface::localCartesianToGlobal(
    const GeometryContext& gctx, const Vector2& lposition) const {
  Vector3 loc3Dframe(lposition[Acts::eBoundLoc0], lposition[Acts::eBoundLoc1],
                     0.);
  return Vector3(transform(gctx) * loc3Dframe);
}

Acts::Vector2 Acts::DiscSurface::globalToLocalCartesian(
    const GeometryContext& gctx, const Vector3& position,
    double /*momentum*/) const {
  Vector3 loc3Dframe = (transform(gctx).inverse()) * position;
  return Vector2(loc3Dframe.x(), loc3Dframe.y());
}

std::string Acts::DiscSurface::name() const {
  return "Acts::DiscSurface";
}

const Acts::SurfaceBounds& Acts::DiscSurface::bounds() const {
  if (m_bounds) {
    return (*(m_bounds.get()));
  }
  return s_noBounds;
}

Acts::Polyhedron Acts::DiscSurface::polyhedronRepresentation(
    const GeometryContext& gctx, size_t lseg) const {
  // Prepare vertices and faces
  std::vector<Vector3> vertices;
  std::vector<Polyhedron::FaceType> faces;
  std::vector<Polyhedron::FaceType> triangularMesh;

  // Understand the disc
  bool fullDisc = m_bounds->coversFullAzimuth();
  bool toCenter = m_bounds->rMin() < s_onSurfaceTolerance;
  // If you have bounds you can create a polyhedron representation
  bool exactPolyhedron = (m_bounds->type() == SurfaceBounds::eDiscTrapezoid);
  if (m_bounds) {
    auto vertices2D = m_bounds->vertices(lseg);
    vertices.reserve(vertices2D.size() + 1);
    Vector3 wCenter(0., 0., 0);
    for (const auto& v2D : vertices2D) {
      vertices.push_back(transform(gctx) * Vector3(v2D.x(), v2D.y(), 0.));
      wCenter += (*vertices.rbegin());
    }
    // These are convex shapes, use the helper method
    // For rings there's a sweet spot when this stops working
    if (m_bounds->type() == SurfaceBounds::eDiscTrapezoid or toCenter or
        not fullDisc) {
      // Transform them into the vertex frame
      wCenter *= 1. / vertices.size();
      vertices.push_back(wCenter);
      auto facesMesh = detail::FacesHelper::convexFaceMesh(vertices, true);
      faces = facesMesh.first;
      triangularMesh = facesMesh.second;
    } else {
      // Two concentric rings, we use the pure concentric method momentarily,
      // but that creates too  many unneccesarry faces, when only two
      // are needed to descibe the mesh, @todo investigate merging flag
      auto facesMesh = detail::FacesHelper::cylindricalFaceMesh(vertices, true);
      faces = facesMesh.first;
      triangularMesh = facesMesh.second;
    }
  } else {
    throw std::domain_error(
        "Polyhedron repr of boundless surface not possible.");
  }
  return Polyhedron(vertices, faces, triangularMesh, exactPolyhedron);
}

Acts::Vector2 Acts::DiscSurface::localPolarToCartesian(
    const Vector2& lpolar) const {
  return Vector2(lpolar[eBoundLoc0] * cos(lpolar[eBoundLoc1]),
                 lpolar[eBoundLoc0] * sin(lpolar[eBoundLoc1]));
}

Acts::Vector2 Acts::DiscSurface::localCartesianToPolar(
    const Vector2& lcart) const {
  return Vector2(std::hypot(lcart[eBoundLoc0], lcart[eBoundLoc1]),
                 std::atan2(lcart[eBoundLoc1], lcart[eBoundLoc0]));
}

Acts::BoundToFreeMatrix Acts::DiscSurface::boundToFreeJacobian(
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
  // special polar coordinates for the Disc
  double lrad = boundParams[eBoundLoc0];
  double lphi = boundParams[eBoundLoc1];
  double lcos_phi = cos(lphi);
  double lsin_phi = sin(lphi);
  // retrieve the reference frame
  const auto rframe = referenceFrame(gctx, position, direction);
  // Initialize the jacobian from local to global
  BoundToFreeMatrix jacToGlobal = BoundToFreeMatrix::Zero();
  // the local error components - rotated from reference frame
  jacToGlobal.block<3, 1>(eFreePos0, eBoundLoc0) =
      lcos_phi * rframe.block<3, 1>(0, 0) + lsin_phi * rframe.block<3, 1>(0, 1);
  jacToGlobal.block<3, 1>(eFreePos0, eBoundLoc1) =
      lrad * (lcos_phi * rframe.block<3, 1>(0, 1) -
              lsin_phi * rframe.block<3, 1>(0, 0));
  // the time component
  jacToGlobal(eFreeTime, eBoundTime) = 1;
  // the momentum components
  jacToGlobal(eFreeDir0, eBoundPhi) = (-sin_theta) * sin_phi;
  jacToGlobal(eFreeDir0, eBoundTheta) = cos_theta * cos_phi;
  jacToGlobal(eFreeDir1, eBoundPhi) = sin_theta * cos_phi;
  jacToGlobal(eFreeDir1, eBoundTheta) = cos_theta * sin_phi;
  jacToGlobal(eFreeDir2, eBoundTheta) = (-sin_theta);
  jacToGlobal(eFreeQOverP, eBoundQOverP) = 1;
  return jacToGlobal;
}

Acts::FreeToBoundMatrix Acts::DiscSurface::freeToBoundJacobian(
    const GeometryContext& gctx, const FreeVector& parameters) const {
  using VectorHelpers::perp;
  using VectorHelpers::phi;
  // The global position
  const auto position = parameters.segment<3>(eFreePos0);
  // The direction
  const auto direction = parameters.segment<3>(eFreeDir0);
  // Optimized trigonometry on the propagation direction
  const double x = direction(0);  // == cos(phi) * sin(theta)
  const double y = direction(1);  // == sin(phi) * sin(theta)
  const double z = direction(2);  // == cos(theta)
  // can be turned into cosine/sine
  const double cosTheta = z;
  const double sinTheta = std::hypot(x, y);
  const double invSinTheta = 1. / sinTheta;
  const double cosPhi = x * invSinTheta;
  const double sinPhi = y * invSinTheta;
  // The measurement frame of the surface
  RotationMatrix3 rframeT =
      referenceFrame(gctx, position, direction).transpose();
  // calculate the transformation to local coorinates
  const Vector3 pos_loc = transform(gctx).inverse() * position;
  const double lr = perp(pos_loc);
  const double lphi = phi(pos_loc);
  const double lcphi = cos(lphi);
  const double lsphi = sin(lphi);
  // rotate into the polar coorindates
  auto lx = rframeT.block<1, 3>(0, 0);
  auto ly = rframeT.block<1, 3>(1, 0);
  // Initalize the jacobian from global to local
  FreeToBoundMatrix jacToLocal = FreeToBoundMatrix::Zero();
  // Local position component
  jacToLocal.block<1, 3>(eBoundLoc0, eFreePos0) = lcphi * lx + lsphi * ly;
  jacToLocal.block<1, 3>(eBoundLoc1, eFreePos0) =
      (lcphi * ly - lsphi * lx) / lr;
  // Time element
  jacToLocal(eBoundTime, eFreeTime) = 1;
  // Directional and momentum elements for reference frame surface
  jacToLocal(eBoundPhi, eFreeDir0) = -sinPhi * invSinTheta;
  jacToLocal(eBoundPhi, eFreeDir1) = cosPhi * invSinTheta;
  jacToLocal(eBoundTheta, eFreeDir0) = cosPhi * cosTheta;
  jacToLocal(eBoundTheta, eFreeDir1) = sinPhi * cosTheta;
  jacToLocal(eBoundTheta, eFreeDir2) = -sinTheta;
  jacToLocal(eBoundQOverP, eFreeQOverP) = 1;
  return jacToLocal;
}

Acts::SurfaceIntersection Acts::DiscSurface::intersect(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction, const BoundaryCheck& bcheck,
    ActsScalar tolerance) const {
  // Get the contextual transform
  auto gctxTransform = transform(gctx);
  // Use the intersection helper for planar surfaces
  auto intersection =
      PlanarHelper::intersect(gctxTransform, position, direction, tolerance);
  // Evaluate boundary check if requested (and reachable)
  if (intersection.status != Intersection3D::Status::unreachable and bcheck and
      m_bounds != nullptr) {
    // Built-in local to global for speed reasons
    const auto& tMatrix = gctxTransform.matrix();
    const Vector3 vecLocal(intersection.position - tMatrix.block<3, 1>(0, 3));
    const Vector2 lcartesian = tMatrix.block<3, 2>(0, 0).transpose() * vecLocal;
    if (bcheck.type() == BoundaryCheck::Type::eAbsolute and
        m_bounds->coversFullAzimuth()) {
      double modifiedTolerance = tolerance + bcheck.tolerance()[eBoundLoc0];
      if (not m_bounds->insideRadialBounds(VectorHelpers::perp(lcartesian),
                                           modifiedTolerance)) {
        intersection.status = Intersection3D::Status::missed;
      }
    } else if (not insideBounds(localCartesianToPolar(lcartesian), bcheck)) {
      intersection.status = Intersection3D::Status::missed;
    }
  }
  return {intersection, this};
}

Acts::ActsMatrix<2, 3> Acts::DiscSurface::localCartesianToBoundLocalDerivative(
    const GeometryContext& gctx, const Vector3& position) const {
  using VectorHelpers::perp;
  using VectorHelpers::phi;
  // The local frame transform
  const auto& sTransform = transform(gctx);
  // calculate the transformation to local coorinates
  const Vector3 localPos = sTransform.inverse() * position;
  const double lr = perp(localPos);
  const double lphi = phi(localPos);
  const double lcphi = std::cos(lphi);
  const double lsphi = std::sin(lphi);
  ActsMatrix<2, 3> loc3DToLocBound = ActsMatrix<2, 3>::Zero();
  loc3DToLocBound << lcphi, lsphi, 0, -lsphi / lr, lcphi / lr, 0;

  return loc3DToLocBound;
}

Acts::Vector3 Acts::DiscSurface::normal(const GeometryContext& gctx,
                                        const Vector2& /*lposition*/) const {
  // fast access via tranform matrix (and not rotation())
  const auto& tMatrix = transform(gctx).matrix();
  return Vector3(tMatrix(0, 2), tMatrix(1, 2), tMatrix(2, 2));
}

Acts::Vector3 Acts::DiscSurface::binningPosition(const GeometryContext& gctx,
                                                 BinningValue bValue) const {
  if (bValue == binR || bValue == binPhi) {
    double r = m_bounds->binningValueR();
    double phi = m_bounds->binningValuePhi();
    return localToGlobal(gctx, Vector2{r, phi}, Vector3{});
  }
  return center(gctx);
}

double Acts::DiscSurface::binningPositionValue(const GeometryContext& gctx,
                                               BinningValue bValue) const {
  if (bValue == binR) {
    return VectorHelpers::perp(binningPosition(gctx, bValue));
  }
  if (bValue == binPhi) {
    return VectorHelpers::phi(binningPosition(gctx, bValue));
  }

  return GeometryObject::binningPositionValue(gctx, bValue);
}

double Acts::DiscSurface::pathCorrection(const GeometryContext& gctx,
                                         const Vector3& position,
                                         const Vector3& direction) const {
  /// we can ignore the global position here
  return 1. / std::abs(Surface::normal(gctx, position).dot(direction));
}
