// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/DiscSurface.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <utility>
#include <vector>

#include "Acts/Surfaces/DiscTrapezoidBounds.hpp"
#include "Acts/Surfaces/InfiniteBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/detail/FacesHelper.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/ThrowAssert.hpp"

using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;

Acts::DiscSurface::DiscSurface(const DiscSurface& other)
    : GeometryObject(), Surface(other), m_bounds(other.m_bounds) {}

Acts::DiscSurface::DiscSurface(const GeometryContext& gctx,
                               const DiscSurface& other,
                               const Transform3D& transf)
    : GeometryObject(),
      Surface(gctx, other, transf),
      m_bounds(other.m_bounds) {}

Acts::DiscSurface::DiscSurface(std::shared_ptr<const Transform3D> htrans,
                               double rmin, double rmax, double hphisec)
    : GeometryObject(),
      Surface(std::move(htrans)),
      m_bounds(std::make_shared<const RadialBounds>(rmin, rmax, hphisec)) {}

Acts::DiscSurface::DiscSurface(std::shared_ptr<const Transform3D> htrans,
                               double minhalfx, double maxhalfx, double maxR,
                               double minR, double avephi, double stereo)
    : GeometryObject(),
      Surface(std::move(htrans)),
      m_bounds(std::make_shared<const DiscTrapezoidBounds>(
          minhalfx, maxhalfx, maxR, minR, avephi, stereo)) {}

Acts::DiscSurface::DiscSurface(std::shared_ptr<const Transform3D> htrans,
                               std::shared_ptr<const DiscBounds> dbounds)
    : GeometryObject(),
      Surface(std::move(htrans)),
      m_bounds(std::move(dbounds)) {}

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

void Acts::DiscSurface::localToGlobal(const GeometryContext& gctx,
                                      const Vector2D& lposition,
                                      const Vector3D& /*gmom*/,
                                      Vector3D& position) const {
  // create the position in the local 3d frame
  Vector3D loc3Dframe(lposition[Acts::eLOC_R] * cos(lposition[Acts::eLOC_PHI]),
                      lposition[Acts::eLOC_R] * sin(lposition[Acts::eLOC_PHI]),
                      0.);
  // transport it to the globalframe (very unlikely that this is not needed)
  position = transform(gctx) * loc3Dframe;
}

bool Acts::DiscSurface::globalToLocal(const GeometryContext& gctx,
                                      const Vector3D& position,
                                      const Vector3D& /*gmom*/,
                                      Vector2D& lposition) const {
  // transport it to the globalframe (very unlikely that this is not needed)
  Vector3D loc3Dframe = (transform(gctx).inverse()) * position;
  lposition = Acts::Vector2D(perp(loc3Dframe), phi(loc3Dframe));
  return ((std::abs(loc3Dframe.z()) > s_onSurfaceTolerance) ? false : true);
}

const Acts::Vector2D Acts::DiscSurface::localPolarToLocalCartesian(
    const Vector2D& locpol) const {
  const DiscTrapezoidBounds* dtbo =
      dynamic_cast<const Acts::DiscTrapezoidBounds*>(&(bounds()));
  if (dtbo != nullptr) {
    double rMedium = dtbo->rCenter();
    double phi = dtbo->averagePhi();

    Vector2D polarCenter(rMedium, phi);
    Vector2D cartCenter = localPolarToCartesian(polarCenter);
    Vector2D cartPos = localPolarToCartesian(locpol);
    Vector2D Pos = cartPos - cartCenter;

    Acts::Vector2D locPos(
        Pos[Acts::eLOC_X] * sin(phi) - Pos[Acts::eLOC_Y] * cos(phi),
        Pos[Acts::eLOC_Y] * sin(phi) + Pos[Acts::eLOC_X] * cos(phi));
    return Vector2D(locPos[Acts::eLOC_X], locPos[Acts::eLOC_Y]);
  }
  return Vector2D(locpol[Acts::eLOC_R] * cos(locpol[Acts::eLOC_PHI]),
                  locpol[Acts::eLOC_R] * sin(locpol[Acts::eLOC_PHI]));
}

const Acts::Vector3D Acts::DiscSurface::localCartesianToGlobal(
    const GeometryContext& gctx, const Vector2D& lposition) const {
  Vector3D loc3Dframe(lposition[Acts::eLOC_X], lposition[Acts::eLOC_Y], 0.);
  return Vector3D(transform(gctx) * loc3Dframe);
}

const Acts::Vector2D Acts::DiscSurface::globalToLocalCartesian(
    const GeometryContext& gctx, const Vector3D& position,
    double /*unused*/) const {
  Vector3D loc3Dframe = (transform(gctx).inverse()) * position;
  return Vector2D(loc3Dframe.x(), loc3Dframe.y());
}

std::string Acts::DiscSurface::name() const {
  return "Acts::DiscSurface";
}

std::shared_ptr<Acts::DiscSurface> Acts::DiscSurface::clone(
    const GeometryContext& gctx, const Transform3D& shift) const {
  return std::shared_ptr<DiscSurface>(this->clone_impl(gctx, shift));
}

Acts::DiscSurface* Acts::DiscSurface::clone_impl(
    const GeometryContext& gctx, const Transform3D& shift) const {
  return new DiscSurface(gctx, *this, shift);
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
  std::vector<Vector3D> vertices;
  std::vector<Polyhedron::Face> faces;
  std::vector<Polyhedron::Face> triangularMesh;

  // Understand the disc
  bool fullDisc = m_bounds->coversFullAzimuth();
  bool toCenter = m_bounds->rMin() < s_onSurfaceTolerance;
  // If you have bounds you can create a polyhedron representation
  if (m_bounds) {
    auto vertices2D = m_bounds->vertices(lseg);
    vertices.reserve(vertices2D.size() + 1);
    Vector3D wCenter(0., 0., 0);
    for (const auto& v2D : vertices2D) {
      vertices.push_back(transform(gctx) * Vector3D(v2D.x(), v2D.y(), 0.));
      wCenter += (*vertices.rbegin());
    }
    // These are convex shapes, use the helper method
    // For rings there's a sweet spot when this stops working
    if (m_bounds->type() == SurfaceBounds::DiscTrapezoidal or toCenter or
        not fullDisc) {
      // Transformt hem into the vertex frame
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
  return Polyhedron(vertices, faces, triangularMesh);
}
