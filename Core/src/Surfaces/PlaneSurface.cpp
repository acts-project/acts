// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/PlaneSurface.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>

#include "Acts/Surfaces/EllipseBounds.hpp"
#include "Acts/Surfaces/InfiniteBounds.hpp"
#include "Acts/Surfaces/PlanarBooleanBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/detail/FacesHelper.hpp"
#include "Acts/Utilities/ThrowAssert.hpp"

Acts::PlaneSurface::PlaneSurface(const PlaneSurface& other)
    : GeometryObject(), Surface(other), m_bounds(other.m_bounds) {}

Acts::PlaneSurface::PlaneSurface(const GeometryContext& gctx,
                                 const PlaneSurface& other,
                                 const Transform3D& transf)
    : GeometryObject(),
      Surface(gctx, other, transf),
      m_bounds(other.m_bounds) {}

Acts::PlaneSurface::PlaneSurface(const Vector3D& center, const Vector3D& normal)
    : Surface(), m_bounds(nullptr) {
  /// the right-handed coordinate system is defined as
  /// T = normal
  /// U = Z x T if T not parallel to Z otherwise U = X x T
  /// V = T x U
  Vector3D T = normal.normalized();
  Vector3D U = std::abs(T.dot(Vector3D::UnitZ())) < s_curvilinearProjTolerance
                   ? Vector3D::UnitZ().cross(T).normalized()
                   : Vector3D::UnitX().cross(T).normalized();
  Vector3D V = T.cross(U);
  RotationMatrix3D curvilinearRotation;
  curvilinearRotation.col(0) = U;
  curvilinearRotation.col(1) = V;
  curvilinearRotation.col(2) = T;

  // curvilinear surfaces are boundless
  Transform3D transform{curvilinearRotation};
  transform.pretranslate(center);
  Surface::m_transform = std::make_shared<const Transform3D>(transform);
}

Acts::PlaneSurface::PlaneSurface(
    const std::shared_ptr<const PlanarBounds>& pbounds,
    const Acts::DetectorElementBase& detelement)
    : Surface(detelement), m_bounds(pbounds) {
  /// surfaces representing a detector element must have bounds
  throw_assert(pbounds, "PlaneBounds must not be nullptr");
}

Acts::PlaneSurface::PlaneSurface(std::shared_ptr<const Transform3D> htrans,
                                 std::shared_ptr<const PlanarBounds> pbounds)
    : Surface(std::move(htrans)), m_bounds(std::move(pbounds)) {}

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

void Acts::PlaneSurface::localToGlobal(const GeometryContext& gctx,
                                       const Vector2D& lposition,
                                       const Vector3D& /*gmom*/,
                                       Vector3D& position) const {
  Vector3D loc3Dframe(lposition[Acts::eLOC_X], lposition[Acts::eLOC_Y], 0.);
  /// the chance that there is no transform is almost 0, let's apply it
  position = transform(gctx) * loc3Dframe;
}

bool Acts::PlaneSurface::globalToLocal(const GeometryContext& gctx,
                                       const Vector3D& position,
                                       const Vector3D& /*gmom*/,
                                       Acts::Vector2D& lposition) const {
  /// the chance that there is no transform is almost 0, let's apply it
  Vector3D loc3Dframe = (transform(gctx).inverse()) * position;
  lposition = Vector2D(loc3Dframe.x(), loc3Dframe.y());
  return ((loc3Dframe.z() * loc3Dframe.z() >
           s_onSurfaceTolerance * s_onSurfaceTolerance)
              ? false
              : true);
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
  // If you have bounds you can create a polyhedron representation
  if (m_bounds == nullptr) {
    throw std::domain_error(
        "Polyhedron repr of boundless surface not possible.");
  }

  auto baseTransform = transform(gctx);
  if (m_bounds->type() == SurfaceBounds::ePlanarBoolean) {
    auto booleanBounds =
        dynamic_cast<const PlanarBooleanBounds*>(m_bounds.get());
    if (booleanBounds->bOperation() == eUnion) {
      auto phedron =
          createPolyhedron(booleanBounds->leftBounds(), baseTransform, lseg,
                           booleanBounds->leftShift());
      phedron.merge(createPolyhedron(booleanBounds->rightBounds(),
                                     baseTransform, lseg,
                                     booleanBounds->rightShift()));
      return phedron;
    }
  }

  return createPolyhedron(*m_bounds, baseTransform, lseg);
}

Acts::Polyhedron Acts::PlaneSurface::createPolyhedron(
    const PlanarBounds& pBounds, const Transform3D& bTransform, size_t lseg,
    const Vector2D& shift) const {
  // Prepare vertices and faces
  std::vector<Vector3D> vertices;
  std::vector<Polyhedron::FaceType> faces;
  std::vector<Polyhedron::FaceType> triangularMesh;

  bool exactPolyhedron = true;

  auto vertices2D = pBounds.vertices(lseg);
  vertices.reserve(vertices2D.size() + 1);
  for (const auto& v2D : vertices2D) {
    vertices.push_back(bTransform *
                       Vector3D(v2D.x() + shift.x(), v2D.y() + shift.y(), 0.));
  }
  bool isEllipse = pBounds.type() == SurfaceBounds::eEllipse;
  bool innerExists = false, coversFull = false;
  if (isEllipse) {
    exactPolyhedron = false;
    auto vStore = pBounds.values();
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

  return Polyhedron(vertices, faces, triangularMesh, exactPolyhedron);
}