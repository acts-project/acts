// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// DiscSurface.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Surfaces/DiscSurface.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <utility>

#include "Acts/Surfaces/DiscTrapezoidalBounds.hpp"
#include "Acts/Surfaces/InfiniteBounds.hpp"
#include "Acts/Surfaces/PolyhedronRepresentation.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/ThrowAssert.hpp"

using Acts::VectorHelpers::phi;
using Acts::VectorHelpers::perp;

Acts::DiscSurface::DiscSurface(const DiscSurface& other)
  : GeometryObject(), Surface(other), m_bounds(other.m_bounds)
{
}

Acts::DiscSurface::DiscSurface(const GeometryContext& gctx,
                               const DiscSurface&     other,
                               const Transform3D&     transf)
  : GeometryObject(), Surface(gctx, other, transf), m_bounds(other.m_bounds)
{
}

Acts::DiscSurface::DiscSurface(std::shared_ptr<const Transform3D> htrans,
                               double                             rmin,
                               double                             rmax,
                               double                             hphisec)
  : GeometryObject()
  , Surface(std::move(htrans))
  , m_bounds(std::make_shared<const RadialBounds>(rmin, rmax, hphisec))
{
}

Acts::DiscSurface::DiscSurface(std::shared_ptr<const Transform3D> htrans,
                               double                             minhalfx,
                               double                             maxhalfx,
                               double                             maxR,
                               double                             minR,
                               double                             avephi,
                               double                             stereo)
  : GeometryObject()
  , Surface(std::move(htrans))
  , m_bounds(std::make_shared<const DiscTrapezoidalBounds>(minhalfx,
                                                           maxhalfx,
                                                           maxR,
                                                           minR,
                                                           avephi,
                                                           stereo))
{
}

Acts::DiscSurface::DiscSurface(std::shared_ptr<const Transform3D> htrans,
                               std::shared_ptr<const DiscBounds>  dbounds)
  : GeometryObject(), Surface(std::move(htrans)), m_bounds(std::move(dbounds))
{
}

Acts::DiscSurface::DiscSurface(const std::shared_ptr<const DiscBounds>& dbounds,
                               const DetectorElementBase& detelement)
  : GeometryObject(), Surface(detelement), m_bounds(dbounds)
{
  throw_assert(dbounds, "nullptr as DiscBounds");
}

Acts::DiscSurface&
Acts::DiscSurface::operator=(const DiscSurface& other)
{
  if (this != &other) {
    Acts::Surface::operator=(other);
    m_bounds               = other.m_bounds;
  }
  return *this;
}

Acts::Surface::SurfaceType
Acts::DiscSurface::type() const
{
  return Surface::Disc;
}

void
Acts::DiscSurface::localToGlobal(const GeometryContext& gctx,
                                 const Vector2D&        lpos,
                                 const Vector3D& /*gmom*/,
                                 Vector3D& gpos) const
{
  // create the position in the local 3d frame
  Vector3D loc3Dframe(lpos[Acts::eLOC_R] * cos(lpos[Acts::eLOC_PHI]),
                      lpos[Acts::eLOC_R] * sin(lpos[Acts::eLOC_PHI]),
                      0.);
  // transport it to the globalframe (very unlikely that this is not needed)
  gpos = transform(gctx) * loc3Dframe;
}

bool
Acts::DiscSurface::globalToLocal(const GeometryContext& gctx,
                                 const Vector3D&        gpos,
                                 const Vector3D& /*gmom*/,
                                 Vector2D& lpos) const
{
  // transport it to the globalframe (very unlikely that this is not needed)
  Vector3D loc3Dframe = (transform(gctx).inverse()) * gpos;
  lpos                = Acts::Vector2D(perp(loc3Dframe), phi(loc3Dframe));
  return ((std::abs(loc3Dframe.z()) > s_onSurfaceTolerance) ? false : true);
}

const Acts::Vector2D
Acts::DiscSurface::localPolarToLocalCartesian(const Vector2D& locpol) const
{
  const DiscTrapezoidalBounds* dtbo
      = dynamic_cast<const Acts::DiscTrapezoidalBounds*>(&(bounds()));
  if (dtbo != nullptr) {
    double rMedium = dtbo->rCenter();
    double phi     = dtbo->averagePhi();

    Vector2D polarCenter(rMedium, phi);
    Vector2D cartCenter = localPolarToCartesian(polarCenter);
    Vector2D cartPos    = localPolarToCartesian(locpol);
    Vector2D Pos        = cartPos - cartCenter;

    Acts::Vector2D locPos(
        Pos[Acts::eLOC_X] * sin(phi) - Pos[Acts::eLOC_Y] * cos(phi),
        Pos[Acts::eLOC_Y] * sin(phi) + Pos[Acts::eLOC_X] * cos(phi));
    return Vector2D(locPos[Acts::eLOC_X], locPos[Acts::eLOC_Y]);
  }
  return Vector2D(locpol[Acts::eLOC_R] * cos(locpol[Acts::eLOC_PHI]),
                  locpol[Acts::eLOC_R] * sin(locpol[Acts::eLOC_PHI]));
}

const Acts::Vector3D
Acts::DiscSurface::localCartesianToGlobal(const GeometryContext& gctx,
                                          const Vector2D&        lpos) const
{
  Vector3D loc3Dframe(lpos[Acts::eLOC_X], lpos[Acts::eLOC_Y], 0.);
  return Vector3D(transform(gctx) * loc3Dframe);
}

const Acts::Vector2D
Acts::DiscSurface::globalToLocalCartesian(const GeometryContext& gctx,
                                          const Vector3D&        gpos,
                                          double /*unused*/) const
{
  Vector3D loc3Dframe = (transform(gctx).inverse()) * gpos;
  return Vector2D(loc3Dframe.x(), loc3Dframe.y());
}

std::string
Acts::DiscSurface::name() const
{
  return "Acts::DiscSurface";
}

std::shared_ptr<Acts::DiscSurface>
Acts::DiscSurface::clone(const GeometryContext& gctx,
                         const Transform3D&     shift) const
{
  return std::shared_ptr<DiscSurface>(this->clone_impl(gctx, shift));
}

Acts::DiscSurface*
Acts::DiscSurface::clone_impl(const GeometryContext& gctx,
                              const Transform3D&     shift) const
{
  return new DiscSurface(gctx, *this, shift);
}

const Acts::SurfaceBounds&
Acts::DiscSurface::bounds() const
{
  if (m_bounds) {
    return (*(m_bounds.get()));
  }
  return s_noBounds;
}

Acts::PolyhedronRepresentation
Acts::DiscSurface::polyhedronRepresentation(const GeometryContext& gctx,
                                            size_t                 l0div,
                                            size_t /*unused*/) const
{
  std::vector<Vector3D>            vertices;
  std::vector<std::vector<size_t>> faces;

  if (l0div < 3) {
    throw std::domain_error("Polyhedron repr of disk with <3 div is undefined");
  }

  auto bounds = std::dynamic_pointer_cast<const RadialBounds>(m_bounds);
  if (!bounds) {
    throw std::domain_error(
        "Polyhedron repr of disk with non RadialBounds currently unsupported");
  }

  double phistep = 2 * M_PI / l0div;
  double rMin    = bounds->rMin();
  double rMax    = bounds->rMax();

  Vector3D inner(rMin, 0, 0);
  Vector3D outer(rMax, 0, 0);

  const Transform3D& sfTransform = transform(gctx);

  for (size_t i = 0; i < l0div; i++) {
    Transform3D rot(AngleAxis3D(i * phistep, Vector3D::UnitZ()));
    vertices.push_back(sfTransform * rot * inner);
    vertices.push_back(sfTransform * rot * outer);
  }

  for (size_t v = 0; v < vertices.size() - 2; v = v + 2) {

    faces.push_back({v, v + 1, v + 3, v + 2});
  }
  if (l0div > 2) {
    faces.push_back({vertices.size() - 2, vertices.size() - 1, 1, 0});
  }

  return PolyhedronRepresentation(vertices, faces);
}
