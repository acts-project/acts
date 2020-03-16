// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/PerigeeSurface.hpp"

#include <iomanip>
#include <iostream>
#include <utility>

Acts::PerigeeSurface::PerigeeSurface(const Vector3D& gp)
    : LineSurface(nullptr, nullptr) {
  Surface::m_transform = std::make_shared<const Transform3D>(
      Translation3D(gp.x(), gp.y(), gp.z()));
}

Acts::PerigeeSurface::PerigeeSurface(
    std::shared_ptr<const Transform3D> tTransform)
    : GeometryObject(), LineSurface(std::move(tTransform)) {}

Acts::PerigeeSurface::PerigeeSurface(const PerigeeSurface& other)
    : GeometryObject(), LineSurface(other) {}

Acts::PerigeeSurface::PerigeeSurface(const GeometryContext& gctx,
                                     const PerigeeSurface& other,
                                     const Transform3D& transf)
    : GeometryObject(), LineSurface(gctx, other, transf) {}

Acts::PerigeeSurface& Acts::PerigeeSurface::operator=(
    const PerigeeSurface& other) {
  if (this != &other) {
    LineSurface::operator=(other);
  }
  return *this;
}

std::shared_ptr<Acts::PerigeeSurface> Acts::PerigeeSurface::clone(
    const GeometryContext& gctx, const Transform3D& shift) const {
  return std::shared_ptr<PerigeeSurface>(this->clone_impl(gctx, shift));
}

Acts::PerigeeSurface* Acts::PerigeeSurface::clone_impl(
    const GeometryContext& gctx, const Transform3D& shift) const {
  return new PerigeeSurface(gctx, *this, shift);
}

Acts::Surface::SurfaceType Acts::PerigeeSurface::type() const {
  return Surface::Perigee;
}

std::string Acts::PerigeeSurface::name() const {
  return "Acts::PerigeeSurface";
}

std::ostream& Acts::PerigeeSurface::toStream(const GeometryContext& gctx,
                                             std::ostream& sl) const {
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::PerigeeSurface:" << std::endl;
  const Vector3D& sfCenter = center(gctx);
  sl << "     Center position  (x, y, z) = (" << sfCenter.x() << ", "
     << sfCenter.y() << ", " << sfCenter.z() << ")";
  sl << std::setprecision(-1);
  return sl;
}

Acts::Polyhedron Acts::PerigeeSurface::polyhedronRepresentation(
    const GeometryContext& gctx, size_t /*lseg*/) const {
  // Prepare vertices and faces
  std::vector<Vector3D> vertices;
  std::vector<Polyhedron::Face> faces;
  std::vector<Polyhedron::Face> triangularMesh;

  const Transform3D& ctransform = transform(gctx);
  Vector3D left(0, 0, -100.);
  Vector3D right(0, 0, 100.);

  // The central wire/straw
  vertices.push_back(ctransform * left);
  vertices.push_back(ctransform * right);
  faces.push_back({0, 1});
  vertices.push_back(ctransform * Vector3D(0., 0., 0.));
  triangularMesh.push_back({0, 2, 1});

  return Polyhedron(vertices, faces, triangularMesh);
}
