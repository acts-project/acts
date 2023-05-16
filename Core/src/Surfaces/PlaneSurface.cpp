// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/PlaneSurface.hpp"

#include "Acts/Surfaces/EllipseBounds.hpp"
#include "Acts/Surfaces/InfiniteBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/SurfaceError.hpp"
#include "Acts/Surfaces/detail/FacesHelper.hpp"
#include "Acts/Surfaces/detail/PlanarHelper.hpp"
#include "Acts/Utilities/ThrowAssert.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>

Acts::PlaneSurface::PlaneSurface(const PlaneSurface& other)
    : GeometryObject(), Surface(other), m_bounds(other.m_bounds) {}

Acts::PlaneSurface::PlaneSurface(const GeometryContext& gctx,
                                 const PlaneSurface& other,
                                 const Transform3& transform)
    : GeometryObject(),
      Surface(gctx, other, transform),
      m_bounds(other.m_bounds) {}

Acts::PlaneSurface::PlaneSurface(const Vector3& center, const Vector3& normal)
    : Surface(), m_bounds(nullptr) {
  /// the right-handed coordinate system is defined as
  /// T = normal
  /// U = Z x T if T not parallel to Z otherwise U = X x T
  /// V = T x U
  Vector3 T = normal.normalized();
  Vector3 U = std::abs(T.dot(Vector3::UnitZ())) < s_curvilinearProjTolerance
                  ? Vector3::UnitZ().cross(T).normalized()
                  : Vector3::UnitX().cross(T).normalized();
  Vector3 V = T.cross(U);
  RotationMatrix3 curvilinearRotation;
  curvilinearRotation.col(0) = U;
  curvilinearRotation.col(1) = V;
  curvilinearRotation.col(2) = T;

  // curvilinear surfaces are boundless
  m_transform = Transform3{curvilinearRotation};
  m_transform.pretranslate(center);
}

Acts::PlaneSurface::PlaneSurface(
    const std::shared_ptr<const PlanarBounds>& pbounds,
    const Acts::DetectorElementBase& detelement)
    : Surface(detelement), m_bounds(pbounds) {
  /// surfaces representing a detector element must have bounds
  throw_assert(pbounds, "PlaneBounds must not be nullptr");
}

Acts::PlaneSurface::PlaneSurface(const Transform3& transform,
                                 std::shared_ptr<const PlanarBounds> pbounds)
    : Surface(transform), m_bounds(std::move(pbounds)) {}

Acts::PlaneSurface& Acts::PlaneSurface::operator=(const PlaneSurface& other) {
  if (this != &other) {
    Surface::operator=(other);
    m_bounds = other.m_bounds;
  }
  return *this;
}

Acts::Surface::SurfaceType Acts::PlaneSurface::type() const {
  return Surface::Plane;
}

Acts::Vector3 Acts::PlaneSurface::localToGlobal(
    const GeometryContext& gctx, const Vector2& lposition,
    const Vector3& /*momentum*/) const {
  return transform(gctx) *
         Vector3(lposition[Acts::eBoundLoc0], lposition[Acts::eBoundLoc1], 0.);
}

Acts::Result<Acts::Vector2> Acts::PlaneSurface::globalToLocal(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& /*momentum*/, double tolerance) const {
  Vector3 loc3Dframe = transform(gctx).inverse() * position;
  if (std::abs(loc3Dframe.z()) > std::abs(tolerance)) {
    return Result<Vector2>::failure(SurfaceError::GlobalPositionNotOnSurface);
  }
  return Result<Vector2>::success({loc3Dframe.x(), loc3Dframe.y()});
}

std::string Acts::PlaneSurface::name() const {
  return "Acts::PlaneSurface";
}

const Acts::SurfaceBounds& Acts::PlaneSurface::bounds() const {
  if (m_bounds) {
    return (*m_bounds.get());
  }
  return s_noBounds;
}

Acts::Polyhedron Acts::PlaneSurface::polyhedronRepresentation(
    const GeometryContext& gctx, size_t lseg) const {
  // Prepare vertices and faces
  std::vector<Vector3> vertices;
  std::vector<Polyhedron::FaceType> faces;
  std::vector<Polyhedron::FaceType> triangularMesh;
  bool exactPolyhedron = true;

  // If you have bounds you can create a polyhedron representation
  if (m_bounds) {
    auto vertices2D = m_bounds->vertices(lseg);
    vertices.reserve(vertices2D.size() + 1);
    for (const auto& v2D : vertices2D) {
      vertices.push_back(transform(gctx) * Vector3(v2D.x(), v2D.y(), 0.));
    }
    bool isEllipse = bounds().type() == SurfaceBounds::eEllipse;
    bool innerExists = false, coversFull = false;
    if (isEllipse) {
      exactPolyhedron = false;
      auto vStore = bounds().values();
      innerExists = vStore[EllipseBounds::eInnerRx] > s_epsilon and
                    vStore[EllipseBounds::eInnerRy] > s_epsilon;
      coversFull =
          std::abs(vStore[EllipseBounds::eHalfPhiSector] - M_PI) < s_epsilon;
    }
    // All of those can be described as convex
    // @todo same as for Discs: coversFull is not the right criterium
    // for triangulation
    if (not isEllipse or not innerExists or not coversFull) {
      auto facesMesh = detail::FacesHelper::convexFaceMesh(vertices);
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

Acts::Vector3 Acts::PlaneSurface::normal(const GeometryContext& gctx,
                                         const Vector2& /*lpos*/) const {
  // fast access via tranform matrix (and not rotation())
  const auto& tMatrix = transform(gctx).matrix();
  return Vector3(tMatrix(0, 2), tMatrix(1, 2), tMatrix(2, 2));
}

Acts::Vector3 Acts::PlaneSurface::binningPosition(
    const GeometryContext& gctx, BinningValue /*bValue*/) const {
  return center(gctx);
}

double Acts::PlaneSurface::pathCorrection(const GeometryContext& gctx,
                                          const Vector3& position,
                                          const Vector3& direction) const {
  // We can ignore the global position here
  return 1. / std::abs(Surface::normal(gctx, position).dot(direction));
}

Acts::SurfaceIntersection Acts::PlaneSurface::intersect(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction, const BoundaryCheck& bcheck,
    ActsScalar tolerance) const {
  // Get the contextual transform
  const auto& gctxTransform = transform(gctx);
  // Use the intersection helper for planar surfaces
  auto intersection =
      PlanarHelper::intersect(gctxTransform, position, direction, tolerance);
  // Evaluate boundary check if requested (and reachable)
  if (intersection.status != Intersection3D::Status::unreachable and bcheck) {
    // Built-in local to global for speed reasons
    const auto& tMatrix = gctxTransform.matrix();
    // Create the reference vector in local
    const Vector3 vecLocal(intersection.position - tMatrix.block<3, 1>(0, 3));
    if (not insideBounds(tMatrix.block<3, 2>(0, 0).transpose() * vecLocal,
                         bcheck)) {
      intersection.status = Intersection3D::Status::missed;
    }
  }
  return {intersection, this};
}

Acts::ActsMatrix<2, 3> Acts::PlaneSurface::localCartesianToBoundLocalDerivative(
    const GeometryContext& /*gctx*/, const Vector3& /*position*/) const {
  const ActsMatrix<2, 3> loc3DToLocBound = ActsMatrix<2, 3>::Identity();
  return loc3DToLocBound;
}
