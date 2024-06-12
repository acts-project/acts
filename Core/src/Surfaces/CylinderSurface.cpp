// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/CylinderSurface.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/GeometryObject.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/SurfaceError.hpp"
#include "Acts/Surfaces/detail/AlignmentHelper.hpp"
#include "Acts/Surfaces/detail/FacesHelper.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/ThrowAssert.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <utility>
#include <vector>

namespace Acts {
class DetectorElementBase;
}  // namespace Acts

using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;

Acts::CylinderSurface::CylinderSurface(const CylinderSurface& other)
    : GeometryObject(), RegularSurface(other), m_bounds(other.m_bounds) {}

Acts::CylinderSurface::CylinderSurface(const GeometryContext& gctx,
                                       const CylinderSurface& other,
                                       const Transform3& shift)
    : GeometryObject(),
      RegularSurface(gctx, other, shift),
      m_bounds(other.m_bounds) {}

Acts::CylinderSurface::CylinderSurface(const Transform3& transform,
                                       double radius, double halfz,
                                       double halfphi, double avphi,
                                       double bevelMinZ, double bevelMaxZ)
    : RegularSurface(transform),
      m_bounds(std::make_shared<const CylinderBounds>(
          radius, halfz, halfphi, avphi, bevelMinZ, bevelMaxZ)) {}

Acts::CylinderSurface::CylinderSurface(
    std::shared_ptr<const CylinderBounds> cbounds,
    const DetectorElementBase& detelement)
    : RegularSurface(detelement), m_bounds(std::move(cbounds)) {
  /// surfaces representing a detector element must have bounds
  throw_assert(m_bounds, "CylinderBounds must not be nullptr");
}

Acts::CylinderSurface::CylinderSurface(
    const Transform3& transform, std::shared_ptr<const CylinderBounds> cbounds)
    : RegularSurface(transform), m_bounds(std::move(cbounds)) {
  throw_assert(m_bounds, "CylinderBounds must not be nullptr");
}

Acts::CylinderSurface& Acts::CylinderSurface::operator=(
    const CylinderSurface& other) {
  if (this != &other) {
    Surface::operator=(other);
    m_bounds = other.m_bounds;
  }
  return *this;
}

// return the binning position for ordering in the BinnedArray
Acts::Vector3 Acts::CylinderSurface::binningPosition(
    const GeometryContext& gctx, BinningValue bValue) const {
  // special binning type for R-type methods
  if (bValue == Acts::binR || bValue == Acts::binRPhi) {
    double R = bounds().get(CylinderBounds::eR);
    double phi = bounds().get(CylinderBounds::eAveragePhi);
    return localToGlobal(gctx, Vector2{phi * R, 0}, Vector3{});
  }
  // give the center as default for all of these binning types
  // binX, binY, binZ, binR, binPhi, binRPhi, binH, binEta
  return center(gctx);
}

// return the measurement frame: it's the tangential plane
Acts::RotationMatrix3 Acts::CylinderSurface::referenceFrame(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& /*direction*/) const {
  RotationMatrix3 mFrame;
  // construct the measurement frame
  // measured Y is the z axis
  Vector3 measY = rotSymmetryAxis(gctx);
  // measured z is the position normalized transverse (in local)
  Vector3 measDepth = normal(gctx, position);
  // measured X is what comoes out of it
  Vector3 measX(measY.cross(measDepth).normalized());
  // assign the columnes
  mFrame.col(0) = measX;
  mFrame.col(1) = measY;
  mFrame.col(2) = measDepth;
  // return the rotation matrix
  return mFrame;
}

Acts::Surface::SurfaceType Acts::CylinderSurface::type() const {
  return Surface::Cylinder;
}

Acts::Vector3 Acts::CylinderSurface::localToGlobal(
    const GeometryContext& gctx, const Vector2& lposition) const {
  // create the position in the local 3d frame
  double r = bounds().get(CylinderBounds::eR);
  double phi = lposition[Acts::eBoundLoc0] / r;
  Vector3 position(r * cos(phi), r * sin(phi), lposition[Acts::eBoundLoc1]);
  return transform(gctx) * position;
}

Acts::Result<Acts::Vector2> Acts::CylinderSurface::globalToLocal(
    const GeometryContext& gctx, const Vector3& position,
    double tolerance) const {
  double inttol = tolerance;
  if (tolerance == s_onSurfaceTolerance) {
    // transform default value!
    // @TODO: check if s_onSurfaceTolerance would do here
    inttol = bounds().get(CylinderBounds::eR) * 0.0001;
  }
  if (inttol < 0.01) {
    inttol = 0.01;
  }
  const Transform3& sfTransform = transform(gctx);
  Transform3 inverseTrans(sfTransform.inverse());
  Vector3 loc3Dframe(inverseTrans * position);
  if (std::abs(perp(loc3Dframe) - bounds().get(CylinderBounds::eR)) > inttol) {
    return Result<Vector2>::failure(SurfaceError::GlobalPositionNotOnSurface);
  }
  return Result<Vector2>::success(
      {bounds().get(CylinderBounds::eR) * phi(loc3Dframe), loc3Dframe.z()});
}

std::string Acts::CylinderSurface::name() const {
  return "Acts::CylinderSurface";
}

Acts::Vector3 Acts::CylinderSurface::normal(
    const GeometryContext& gctx, const Acts::Vector2& lposition) const {
  double phi = lposition[Acts::eBoundLoc0] / m_bounds->get(CylinderBounds::eR);
  Vector3 localNormal(cos(phi), sin(phi), 0.);
  return transform(gctx).linear() * localNormal;
}

Acts::Vector3 Acts::CylinderSurface::normal(
    const GeometryContext& gctx, const Acts::Vector3& position) const {
  const Transform3& sfTransform = transform(gctx);
  // get it into the cylinder frame
  Vector3 pos3D = sfTransform.inverse() * position;
  // set the z coordinate to 0
  pos3D.z() = 0.;
  // normalize and rotate back into global
  return sfTransform.linear() * pos3D.normalized();
}

double Acts::CylinderSurface::pathCorrection(
    const GeometryContext& gctx, const Acts::Vector3& position,
    const Acts::Vector3& direction) const {
  Vector3 normalT = normal(gctx, position);
  double cosAlpha = normalT.dot(direction);
  return std::fabs(1. / cosAlpha);
}

const Acts::CylinderBounds& Acts::CylinderSurface::bounds() const {
  return (*m_bounds.get());
}

Acts::Polyhedron Acts::CylinderSurface::polyhedronRepresentation(
    const GeometryContext& gctx, std::size_t lseg) const {
  auto ctrans = transform(gctx);

  // Prepare vertices and faces
  std::vector<Vector3> vertices = bounds().createCircles(ctrans, lseg);
  std::vector<Polyhedron::FaceType> faces;
  std::vector<Polyhedron::FaceType> triangularMesh;

  bool fullCylinder = bounds().coversFullAzimuth();

  auto facesMesh =
      detail::FacesHelper::cylindricalFaceMesh(vertices, fullCylinder);
  return Polyhedron(vertices, facesMesh.first, facesMesh.second, false);
}

Acts::Vector3 Acts::CylinderSurface::rotSymmetryAxis(
    const GeometryContext& gctx) const {
  // fast access via transform matrix (and not rotation())
  return transform(gctx).matrix().block<3, 1>(0, 2);
}

Acts::detail::RealQuadraticEquation Acts::CylinderSurface::intersectionSolver(
    const Transform3& transform, const Vector3& position,
    const Vector3& direction) const {
  // Solve for radius R
  double R = bounds().get(CylinderBounds::eR);

  // Get the transformation matrtix
  const auto& tMatrix = transform.matrix();
  Vector3 caxis = tMatrix.block<3, 1>(0, 2).transpose();
  Vector3 ccenter = tMatrix.block<3, 1>(0, 3).transpose();

  // Check documentation for explanation
  Vector3 pc = position - ccenter;
  Vector3 pcXcd = pc.cross(caxis);
  Vector3 ldXcd = direction.cross(caxis);
  double a = ldXcd.dot(ldXcd);
  double b = 2. * (ldXcd.dot(pcXcd));
  double c = pcXcd.dot(pcXcd) - (R * R);
  // And solve the qaudratic equation
  return detail::RealQuadraticEquation(a, b, c);
}

Acts::SurfaceMultiIntersection Acts::CylinderSurface::intersect(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction, const BoundaryCheck& bcheck,
    ActsScalar tolerance) const {
  const auto& gctxTransform = transform(gctx);

  // Solve the quadratic equation
  auto qe = intersectionSolver(gctxTransform, position, direction);

  // If no valid solution return a non-valid surfaceIntersection
  if (qe.solutions == 0) {
    return {{Intersection3D::invalid(), Intersection3D::invalid()}, this};
  }

  // Check the validity of the first solution
  Vector3 solution1 = position + qe.first * direction;
  Intersection3D::Status status1 = std::abs(qe.first) < std::abs(tolerance)
                                       ? Intersection3D::Status::onSurface
                                       : Intersection3D::Status::reachable;

  // Helper method for boundary check
  auto boundaryCheck =
      [&](const Vector3& solution,
          Intersection3D::Status status) -> Intersection3D::Status {
    // No check to be done, return current status
    if (!bcheck.isEnabled()) {
      return status;
    }
    const auto& cBounds = bounds();
    if (cBounds.coversFullAzimuth() &&
        bcheck.type() == BoundaryCheck::Type::eAbsolute) {
      // Project out the current Z value via local z axis
      // Built-in local to global for speed reasons
      const auto& tMatrix = gctxTransform.matrix();
      // Create the reference vector in local
      const Vector3 vecLocal(solution - tMatrix.block<3, 1>(0, 3));
      double cZ = vecLocal.dot(tMatrix.block<3, 1>(0, 2));
      double modifiedTolerance = tolerance + bcheck.tolerance()[eBoundLoc1];
      double hZ = cBounds.get(CylinderBounds::eHalfLengthZ) + modifiedTolerance;
      return std::abs(cZ) < std::abs(hZ) ? status
                                         : Intersection3D::Status::missed;
    }
    return (isOnSurface(gctx, solution, direction, bcheck)
                ? status
                : Intersection3D::Status::missed);
  };
  // Check first solution for boundary compatibility
  status1 = boundaryCheck(solution1, status1);
  // Set the intersection
  Intersection3D first(solution1, qe.first, status1);
  if (qe.solutions == 1) {
    return {{first, first}, this};
  }
  // Check the validity of the second solution
  Vector3 solution2 = position + qe.second * direction;
  Intersection3D::Status status2 = std::abs(qe.second) < std::abs(tolerance)
                                       ? Intersection3D::Status::onSurface
                                       : Intersection3D::Status::reachable;
  // Check first solution for boundary compatibility
  status2 = boundaryCheck(solution2, status2);
  Intersection3D second(solution2, qe.second, status2);
  // Order based on path length
  if (first.pathLength() <= second.pathLength()) {
    return {{first, second}, this};
  }
  return {{second, first}, this};
}

Acts::AlignmentToPathMatrix Acts::CylinderSurface::alignmentToPathDerivative(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction) const {
  assert(isOnSurface(gctx, position, direction, BoundaryCheck(false)));

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
  const auto norm = 1 / (1 - dz * dz);
  // The direction transpose
  const auto& dirRowVec = direction.transpose();
  // The derivative of path w.r.t. the local axes
  // @note The following calculations assume that the intersection of the track
  // with the cylinder always satisfy: perp(localPos) = R
  const auto localXAxisToPath =
      (-2 * norm * (dx * pcRowVec + localPos.x() * dirRowVec)).eval();
  const auto localYAxisToPath =
      (-2 * norm * (dy * pcRowVec + localPos.y() * dirRowVec)).eval();
  const auto localZAxisToPath =
      (-4 * norm * norm * (dx * localPos.x() + dy * localPos.y()) * dz *
       dirRowVec)
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

Acts::ActsMatrix<2, 3>
Acts::CylinderSurface::localCartesianToBoundLocalDerivative(
    const GeometryContext& gctx, const Vector3& position) const {
  using VectorHelpers::perp;
  using VectorHelpers::phi;
  // The local frame transform
  const auto& sTransform = transform(gctx);
  // calculate the transformation to local coordinates
  const Vector3 localPos = sTransform.inverse() * position;
  const double lr = perp(localPos);
  const double lphi = phi(localPos);
  const double lcphi = std::cos(lphi);
  const double lsphi = std::sin(lphi);
  // Solve for radius R
  double R = bounds().get(CylinderBounds::eR);
  ActsMatrix<2, 3> loc3DToLocBound = ActsMatrix<2, 3>::Zero();
  loc3DToLocBound << -R * lsphi / lr, R * lcphi / lr, 0, 0, 0, 1;

  return loc3DToLocBound;
}

std::shared_ptr<Acts::CylinderSurface> Acts::CylinderSurface::merge(
    const GeometryContext& gctx, const CylinderSurface& other,
    BinningValue direction) const {
  Transform3 otherLocal = transform(gctx).inverse() * other.transform(gctx);

  constexpr auto tolerance = s_onSurfaceTolerance;

  // surface cannot have any relative rotation
  if (!otherLocal.linear().isApprox(RotationMatrix3::Identity())) {
    throw std::invalid_argument(
        "CylinderSurface::merge: surfaces have relative rotation");
  }

  auto checkNoBevel = [](const auto& bounds) {
    if (bounds.get(CylinderBounds::eBevelMinZ) != 0.0) {
      throw std::invalid_argument(
          "CylinderVolumeStack requires all volumes to have a bevel angle of "
          "0");
    }

    if (bounds.get(CylinderBounds::eBevelMaxZ) != 0.0) {
      throw std::invalid_argument(
          "CylinderVolumeStack requires all volumes to have a bevel angle of "
          "0");
    }
  };

  checkNoBevel(bounds());
  checkNoBevel(other.bounds());

  // radii need to be identical
  if (std::abs(bounds().get(CylinderBounds::eR) -
               other.bounds().get(CylinderBounds::eR)) > tolerance) {
    throw std::invalid_argument(
        "CylinderSurface::merge: surfaces have different radii");
  }

  ActsScalar r = bounds().get(CylinderBounds::eR);

  // no translation in x/z is allowed
  Vector3 translation = otherLocal.translation();

  if (std::abs(translation[0]) > tolerance ||
      std::abs(translation[1]) > tolerance) {
    throw std::invalid_argument(
        "CylinderSurface::merge: surfaces have relative translation in x/y");
  }

  if (direction == Acts::binZ) {
    // z shift must match the bounds

    ActsScalar hlZ = bounds().get(CylinderBounds::eHalfLengthZ);
    ActsScalar minZ = -hlZ;
    ActsScalar maxZ = hlZ;

    ActsScalar zShift = translation[2];
    ActsScalar otherHlZ = other.bounds().get(CylinderBounds::eHalfLengthZ);
    ActsScalar otherMinZ = -otherHlZ + zShift;
    ActsScalar otherMaxZ = otherHlZ + zShift;

    std::cout << "hlZ: " << hlZ << std::endl;
    std::cout << "minZ: " << minZ << std::endl;
    std::cout << "maxZ: " << maxZ << std::endl;
    std::cout << "zShift: " << zShift << std::endl;
    std::cout << "otherHlZ: " << otherHlZ << std::endl;
    std::cout << "otherMinZ: " << otherMinZ << std::endl;
    std::cout << "otherMaxZ: " << otherMaxZ << std::endl;

    if (std::abs(maxZ - otherMinZ) > tolerance &&
        std::abs(minZ - otherMaxZ) > tolerance) {
      throw std::invalid_argument(
          "CylinderSurface::merge: surfaces have incompatible z bounds");
    }

    ActsScalar newMaxZ = std::max(maxZ, otherMaxZ);
    ActsScalar newMinZ = std::min(minZ, otherMinZ);
    std::cout << " -> newMaxZ " << newMaxZ << std::endl;
    std::cout << " -> newMinZ " << newMinZ << std::endl;
    ActsScalar newHlZ = (newMaxZ - newMinZ) / 2.0;
    ActsScalar newMidZ = (newMaxZ + newMinZ) / 2.0;
    std::cout << " -> " << newMidZ << std::endl;
    std::cout << " -> +-" << newHlZ << std::endl;

    auto newBounds = std::make_shared<CylinderBounds>(r, newHlZ);

    Transform3 newTransform =
        transform(gctx) * Translation3{Vector3::UnitZ() * newMidZ};

    return Surface::makeShared<CylinderSurface>(newTransform, newBounds);

  } else if (direction == Acts::binRPhi) {
    // no z shift is allowed
    if (std::abs(translation[2]) > tolerance) {
      throw std::invalid_argument(
          "CylinderSurface::merge: surfaces have relative translation in z for "
          "rPhi merging");
    }

    ActsScalar hlPhi = bounds().get(CylinderBounds::eHalfPhiSector);
    ActsScalar avgPhi = bounds().get(CylinderBounds::eAveragePhi);
    ActsScalar minPhi = -hlPhi + avgPhi;
    ActsScalar maxPhi = hlPhi + avgPhi;

    ActsScalar otherHlPhi = other.bounds().get(CylinderBounds::eHalfPhiSector);
    ActsScalar otherAvgPhi = other.bounds().get(CylinderBounds::eAveragePhi);
    ActsScalar otherMinPhi = -otherHlPhi + otherAvgPhi;
    ActsScalar otherMaxPhi = otherHlPhi + otherAvgPhi;

    using namespace Acts::UnitLiterals;

    std::cout << "hlPhi: " << hlPhi / 1_degree << std::endl;
    std::cout << "avgPhi: " << avgPhi / 1_degree << std::endl;
    std::cout << "minPhi: " << minPhi / 1_degree << std::endl;
    std::cout << "maxPhi: " << maxPhi / 1_degree << std::endl;
    std::cout << "otherHlPhi: " << otherHlPhi / 1_degree << std::endl;
    std::cout << "otherAvgPhi: " << otherAvgPhi / 1_degree << std::endl;
    std::cout << "otherMinPhi: " << otherMinPhi / 1_degree << std::endl;
    std::cout << "otherMaxPhi: " << otherMaxPhi / 1_degree << std::endl;

    ActsScalar newMaxPhi = std::max(maxPhi, otherMaxPhi);
    ActsScalar newMinPhi = std::min(minPhi, otherMinPhi);
    ActsScalar newHlPhi = (newMaxPhi - newMinPhi) / 2.0;
    ActsScalar newAvgPhi = (newMaxPhi + newMinPhi) / 2.0;

    std::cout << " -> newMaxPhi " << newMaxPhi / 1_degree << std::endl;
    std::cout << " -> newMinPhi " << newMinPhi / 1_degree << std::endl;
    std::cout << " -> " << newAvgPhi / 1_degree << std::endl;
    std::cout << " -> +-" << newHlPhi / 1_degree << std::endl;

    auto newBounds = std::make_shared<CylinderBounds>(
        r, bounds().get(CylinderBounds::eHalfLengthZ), newHlPhi, newAvgPhi);

    return Surface::makeShared<CylinderSurface>(transform(gctx), newBounds);

  } else {
    throw std::invalid_argument("CylinderSurface::merge: invalid direction " +
                                binningValueNames()[direction]);
  }
  return nullptr;
}
