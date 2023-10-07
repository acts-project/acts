// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/ConeSurface.hpp"

#include "Acts/EventData/detail/TransformationBoundToFree.hpp"
#include "Acts/Surfaces/SurfaceError.hpp"
#include "Acts/Surfaces/detail/FacesHelper.hpp"
#include "Acts/Surfaces/detail/VerticesHelper.hpp"
#include "Acts/Utilities/ThrowAssert.hpp"
#include "Acts/Utilities/detail/RealQuadraticEquation.hpp"

#include <cassert>
#include <cmath>

using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;

Acts::ConeSurface::ConeSurface(const ConeSurface& other)
    : GeometryObject(), Surface(other), m_bounds(other.m_bounds) {}

Acts::ConeSurface::ConeSurface(const GeometryContext& gctx,
                               const ConeSurface& other,
                               const Transform3& shift)
    : GeometryObject(), Surface(gctx, other, shift), m_bounds(other.m_bounds) {}

Acts::ConeSurface::ConeSurface(const Transform3& transform, double alpha,
                               bool symmetric)
    : GeometryObject(),
      Surface(transform),
      m_bounds(std::make_shared<const ConeBounds>(alpha, symmetric)) {}

Acts::ConeSurface::ConeSurface(const Transform3& transform, double alpha,
                               double zmin, double zmax, double halfPhi)
    : GeometryObject(),
      Surface(transform),
      m_bounds(std::make_shared<const ConeBounds>(alpha, zmin, zmax, halfPhi)) {
}

Acts::ConeSurface::ConeSurface(const Transform3& transform,
                               std::shared_ptr<const ConeBounds> cbounds)
    : GeometryObject(), Surface(transform), m_bounds(std::move(cbounds)) {
  throw_assert(m_bounds, "ConeBounds must not be nullptr");
}

Acts::Vector3 Acts::ConeSurface::binningPosition(
    const GeometryContext& gctx, Acts::BinningValue bValue) const {
  const Vector3& sfCenter = center(gctx);

  // special binning type for R-type methods
  if (bValue == Acts::binR || bValue == Acts::binRPhi) {
    return Vector3(sfCenter.x() + bounds().r(sfCenter.z()), sfCenter.y(),
                   sfCenter.z());
  }
  // give the center as default for all of these binning types
  // binX, binY, binZ, binR, binPhi, binRPhi, binH, binEta
  return sfCenter;
}

Acts::Surface::SurfaceType Acts::ConeSurface::type() const {
  return Surface::Cone;
}

Acts::ConeSurface& Acts::ConeSurface::operator=(const ConeSurface& other) {
  if (this != &other) {
    Surface::operator=(other);
    m_bounds = other.m_bounds;
  }
  return *this;
}

Acts::Vector3 Acts::ConeSurface::rotSymmetryAxis(
    const GeometryContext& gctx) const {
  return transform(gctx).matrix().block<3, 1>(0, 2);
}

Acts::RotationMatrix3 Acts::ConeSurface::referenceFrame(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& /*momentum*/) const {
  RotationMatrix3 mFrame;
  // construct the measurement frame
  // measured Y is the local z axis
  Vector3 measY = rotSymmetryAxis(gctx);
  // measured z is the position transverse normalized
  Vector3 measDepth = Vector3(position.x(), position.y(), 0.).normalized();
  // measured X is what comoes out of it
  Acts::Vector3 measX(measY.cross(measDepth).normalized());
  // the columnes
  mFrame.col(0) = measX;
  mFrame.col(1) = measY;
  mFrame.col(2) = measDepth;
  // return the rotation matrix
  //!< @todo fold in alpha
  // return it
  return mFrame;
}

Acts::Vector3 Acts::ConeSurface::localToGlobal(
    const GeometryContext& gctx, const Vector2& lposition,
    const Vector3& /*momentum*/) const {
  // create the position in the local 3d frame
  double r = lposition[Acts::eBoundLoc1] * bounds().tanAlpha();
  double phi = lposition[Acts::eBoundLoc0] / r;
  Vector3 loc3Dframe(r * cos(phi), r * sin(phi), lposition[Acts::eBoundLoc1]);
  return transform(gctx) * loc3Dframe;
}

Acts::Result<Acts::Vector2> Acts::ConeSurface::globalToLocal(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& /*momentum*/, double tolerance) const {
  Vector3 loc3Dframe = transform(gctx).inverse() * position;
  double r = loc3Dframe.z() * bounds().tanAlpha();
  if (std::abs(perp(loc3Dframe) - r) > tolerance) {
    return Result<Vector2>::failure(SurfaceError::GlobalPositionNotOnSurface);
  }
  return Result<Acts::Vector2>::success(
      Vector2(r * atan2(loc3Dframe.y(), loc3Dframe.x()), loc3Dframe.z()));
}

double Acts::ConeSurface::pathCorrection(const GeometryContext& gctx,
                                         const Vector3& position,
                                         const Vector3& direction) const {
  // (cos phi cos alpha, sin phi cos alpha, sgn z sin alpha)
  Vector3 posLocal = transform(gctx).inverse() * position;
  double phi = VectorHelpers::phi(posLocal);
  double sgn = posLocal.z() > 0. ? -1. : +1.;
  double cosAlpha = std::cos(bounds().get(ConeBounds::eAlpha));
  double sinAlpha = std::sin(bounds().get(ConeBounds::eAlpha));
  Vector3 normalC(cos(phi) * cosAlpha, sin(phi) * cosAlpha, sgn * sinAlpha);
  normalC = transform(gctx) * normalC;
  // Back to the global frame
  double cAlpha = normalC.dot(direction);
  return std::abs(1. / cAlpha);
}

std::string Acts::ConeSurface::name() const {
  return "Acts::ConeSurface";
}

Acts::Vector3 Acts::ConeSurface::normal(const GeometryContext& gctx,
                                        const Acts::Vector2& lposition) const {
  // (cos phi cos alpha, sin phi cos alpha, sgn z sin alpha)
  double phi = lposition[Acts::eBoundLoc0] /
               (bounds().r(lposition[Acts::eBoundLoc1])),
         sgn = lposition[Acts::eBoundLoc1] > 0 ? -1. : +1.;
  double cosAlpha = std::cos(bounds().get(ConeBounds::eAlpha));
  double sinAlpha = std::sin(bounds().get(ConeBounds::eAlpha));
  Vector3 localNormal(cos(phi) * cosAlpha, sin(phi) * cosAlpha, sgn * sinAlpha);
  return Vector3(transform(gctx).linear() * localNormal);
}

Acts::Vector3 Acts::ConeSurface::normal(const GeometryContext& gctx,
                                        const Acts::Vector3& position) const {
  // get it into the cylinder frame if needed
  // @todo respect opening angle
  Vector3 pos3D = transform(gctx).inverse() * position;
  pos3D.z() = 0;
  return pos3D.normalized();
}

const Acts::ConeBounds& Acts::ConeSurface::bounds() const {
  // is safe because no constructor w/o bounds exists
  return (*m_bounds.get());
}

Acts::Polyhedron Acts::ConeSurface::polyhedronRepresentation(
    const GeometryContext& gctx, size_t lseg) const {
  // Prepare vertices and faces
  std::vector<Vector3> vertices;
  std::vector<Polyhedron::FaceType> faces;
  std::vector<Polyhedron::FaceType> triangularMesh;

  double minZ = bounds().get(ConeBounds::eMinZ);
  double maxZ = bounds().get(ConeBounds::eMaxZ);

  if (minZ == -std::numeric_limits<double>::infinity() or
      maxZ == std::numeric_limits<double>::infinity()) {
    throw std::domain_error(
        "Polyhedron repr of boundless surface not possible");
  }

  auto ctransform = transform(gctx);

  // The tip - created only once and only, if the it's not a cut-off cone
  bool tipExists = false;
  if (minZ * maxZ <= s_onSurfaceTolerance) {
    vertices.push_back(ctransform * Vector3(0., 0., 0.));
    tipExists = true;
  }

  // Cone parameters
  double hPhiSec = bounds().get(ConeBounds::eHalfPhiSector);
  double avgPhi = bounds().get(ConeBounds::eAveragePhi);
  bool fullCone = (hPhiSec == M_PI);

  // Get the phi segments from the helper
  auto phiSegs = fullCone ? detail::VerticesHelper::phiSegments()
                          : detail::VerticesHelper::phiSegments(
                                avgPhi - hPhiSec, avgPhi + hPhiSec,
                                {static_cast<ActsScalar>(avgPhi)});

  // Negative cone if exists
  std::vector<double> coneSides;
  if (std::abs(minZ) > s_onSurfaceTolerance) {
    coneSides.push_back(minZ);
  }
  if (std::abs(maxZ) > s_onSurfaceTolerance) {
    coneSides.push_back(maxZ);
  }
  for (auto& z : coneSides) {
    // Remember the first vertex
    size_t firstIv = vertices.size();
    // Radius and z offset
    double r = std::abs(z) * bounds().tanAlpha();
    Vector3 zoffset(0., 0., z);
    for (unsigned int iseg = 0; iseg < phiSegs.size() - 1; ++iseg) {
      int addon = (iseg == phiSegs.size() - 2 and not fullCone) ? 1 : 0;
      detail::VerticesHelper::createSegment(vertices, {r, r}, phiSegs[iseg],
                                            phiSegs[iseg + 1], lseg, addon,
                                            zoffset, ctransform);
    }
    // Create the faces
    if (tipExists) {
      for (size_t iv = firstIv + 2; iv < vertices.size() + 1; ++iv) {
        size_t one = 0, two = iv - 1, three = iv - 2;
        if (z < 0.) {
          std::swap(two, three);
        }
        faces.push_back({one, two, three});
      }
      // Complete cone if necessary
      if (fullCone) {
        if (z > 0.) {
          faces.push_back({0, firstIv, vertices.size() - 1});
        } else {
          faces.push_back({0, vertices.size() - 1, firstIv});
        }
      }
    }
  }
  // if no tip exists, connect the two bows
  if (tipExists) {
    triangularMesh = faces;
  } else {
    auto facesMesh =
        detail::FacesHelper::cylindricalFaceMesh(vertices, fullCone);
    faces = facesMesh.first;
    triangularMesh = facesMesh.second;
  }
  return Polyhedron(vertices, faces, triangularMesh, false);
}

Acts::detail::RealQuadraticEquation Acts::ConeSurface::intersectionSolver(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction) const {
  // Transform into the local frame
  Transform3 invTrans = transform(gctx).inverse();
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

Acts::SurfaceIntersection Acts::ConeSurface::intersect(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction, const BoundaryCheck& bcheck) const {
  // Solve the quadratic equation
  auto qe = intersectionSolver(gctx, position, direction);

  // If no valid solution return a non-valid surfaceIntersection
  if (qe.solutions == 0) {
    return SurfaceIntersection();
  }

  // Check the validity of the first solution
  Vector3 solution1 = position + qe.first * direction;
  Intersection3D::Status status1 =
      std::abs(qe.first) < std::abs(s_onSurfaceTolerance)
          ? Intersection3D::Status::onSurface
          : Intersection3D::Status::reachable;

  if (bcheck and not isOnSurface(gctx, solution1, direction, bcheck)) {
    status1 = Intersection3D::Status::missed;
  }

  // Check the validity of the second solution
  Vector3 solution2 = position + qe.first * direction;
  Intersection3D::Status status2 =
      std::abs(qe.second) < std::abs(s_onSurfaceTolerance)
          ? Intersection3D::Status::onSurface
          : Intersection3D::Status::reachable;
  if (bcheck and not isOnSurface(gctx, solution2, direction, bcheck)) {
    status2 = Intersection3D::Status::missed;
  }

  const auto& tf = transform(gctx);
  // Set the intersection
  Intersection3D first(tf * solution1, qe.first, status1);
  Intersection3D second(tf * solution2, qe.second, status2);
  SurfaceIntersection cIntersection(first, this);
  // Check one if its valid or neither is valid
  bool check1 = status1 != Intersection3D::Status::missed or
                (status1 == Intersection3D::Status::missed and
                 status2 == Intersection3D::Status::missed);
  // Check and (eventually) go with the first solution
  if ((check1 and (std::abs(qe.first) < std::abs(qe.second))) or
      status2 == Intersection3D::Status::missed) {
    // And add the alternative
    if (qe.solutions > 1) {
      cIntersection.alternative = second;
    }
  } else {
    // And add the alternative
    if (qe.solutions > 1) {
      cIntersection.alternative = first;
    }
    cIntersection.intersection = second;
  }
  return cIntersection;
}

Acts::AlignmentToPathMatrix Acts::ConeSurface::alignmentToPathDerivative(
    const GeometryContext& gctx, const FreeVector& parameters) const {
  // The global position
  const auto position = parameters.segment<3>(eFreePos0);
  // The direction
  const auto direction = parameters.segment<3>(eFreeDir0);
  // The vector between position and center
  const auto pcRowVec = (position - center(gctx)).transpose().eval();
  // The rotation
  const auto& rotation = transform(gctx).rotation();
  // The local frame x/y/z axis
  const auto& localXAxis = rotation.col(0);
  const auto& localYAxis = rotation.col(1);
  const auto& localZAxis = rotation.col(2);
  // The local coordinates
  const auto localPos = (rotation.transpose() * position).eval();
  const auto dx = direction.dot(localXAxis);
  const auto dy = direction.dot(localYAxis);
  const auto dz = direction.dot(localZAxis);
  // The normalization factor
  const auto tanAlpha2 = bounds().tanAlpha() * bounds().tanAlpha();
  const auto norm = 1 / (1 - dz * dz * (1 + tanAlpha2));
  // The direction transpose
  const auto& dirRowVec = direction.transpose();
  // The derivative of path w.r.t. the local axes
  // @note The following calculations assume that the intersection of the track
  // with the cone always satisfy: localPos.z()*tanAlpha =perp(localPos)
  const auto localXAxisToPath =
      (-2 * norm * (dx * pcRowVec + localPos.x() * dirRowVec)).eval();
  const auto localYAxisToPath =
      (-2 * norm * (dy * pcRowVec + localPos.y() * dirRowVec)).eval();
  const auto localZAxisToPath =
      (2 * norm * tanAlpha2 * (dz * pcRowVec + localPos.z() * dirRowVec) -
       4 * norm * norm * (1 + tanAlpha2) *
           (dx * localPos.x() + dy * localPos.y() -
            dz * localPos.z() * tanAlpha2) *
           dz * dirRowVec)
          .eval();
  // Calculate the derivative of local frame axes w.r.t its rotation
  const auto [rotToLocalXAxis, rotToLocalYAxis, rotToLocalZAxis] =
      detail::rotationToLocalAxesDerivative(rotation);
  // Initialize the derivative of propagation path w.r.t. local frame
  // translation (origin) and rotation
  AlignmentToPathMatrix alignToPath = AlignmentToPathMatrix::Zero();
  alignToPath.segment<3>(eAlignmentCenter0) =
      2 * norm * (dx * localXAxis.transpose() + dy * localYAxis.transpose());
  alignToPath.segment<3>(eAlignmentRotation0) =
      localXAxisToPath * rotToLocalXAxis + localYAxisToPath * rotToLocalYAxis +
      localZAxisToPath * rotToLocalZAxis;

  return alignToPath;
}

Acts::ActsMatrix<2, 3> Acts::ConeSurface::localCartesianToBoundLocalDerivative(
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
  // Solve for radius R
  const double R = localPos.z() * bounds().tanAlpha();
  ActsMatrix<2, 3> loc3DToLocBound = ActsMatrix<2, 3>::Zero();
  loc3DToLocBound << -R * lsphi / lr, R * lcphi / lr,
      lphi * bounds().tanAlpha(), 0, 0, 1;

  return loc3DToLocBound;
}
