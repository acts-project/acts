// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/DiscSurface.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/DiscTrapezoidBounds.hpp"
#include "Acts/Surfaces/InfiniteBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/SurfaceError.hpp"
#include "Acts/Surfaces/detail/FacesHelper.hpp"
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
      Surface(std::move(transform)),
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
  if (loc3Dframe.z() * loc3Dframe.z() > tolerance * tolerance) {
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
    double /*unused*/) const {
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
