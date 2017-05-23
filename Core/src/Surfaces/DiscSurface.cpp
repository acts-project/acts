// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// DiscSurface.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Surfaces/DiscSurface.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>

#include "ACTS/Surfaces/DiscTrapezoidalBounds.hpp"
#include "ACTS/Surfaces/InfiniteBounds.hpp"
#include "ACTS/Surfaces/RadialBounds.hpp"
#include "ACTS/Utilities/Definitions.hpp"

Acts::DiscSurface::DiscSurface(const DiscSurface& other)
  : Surface(other), m_bounds(other.m_bounds)
{
}

Acts::DiscSurface::DiscSurface(const DiscSurface& other,
                               const Transform3D& transf)
  : Surface(other, transf), m_bounds(other.m_bounds)
{
}

Acts::DiscSurface::DiscSurface(std::shared_ptr<const Transform3D> htrans,
                               double                             rmin,
                               double                             rmax,
                               double                             hphisec)
  : Surface(htrans)
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
  : Surface(htrans)
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
  : Surface(htrans), m_bounds(dbounds)
{
}

Acts::DiscSurface::DiscSurface(std::shared_ptr<const DiscBounds> dbounds,
                               const DetectorElementBase&        detelement,
                               const Identifier&                 identifier)
  : Surface(detelement, identifier), m_bounds(nullptr)
{
  assert(dbounds);
}

Acts::DiscSurface::~DiscSurface()
{
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
Acts::DiscSurface::localToGlobal(const Vector2D& lpos,
                                 const Vector3D&,
                                 Vector3D& gpos) const
{
  // create the position in the local 3d frame
  Vector3D loc3Dframe(lpos[Acts::eLOC_R] * cos(lpos[Acts::eLOC_PHI]),
                      lpos[Acts::eLOC_R] * sin(lpos[Acts::eLOC_PHI]),
                      0.);
  // transport it to the globalframe (very unlikely that this is not needed)
  gpos = transform() * loc3Dframe;
}

bool
Acts::DiscSurface::globalToLocal(const Acts::Vector3D& gpos,
                                 const Acts::Vector3D&,
                                 Acts::Vector2D& lpos) const
{
  // transport it to the globalframe (very unlikely that this is not needed)
  Vector3D loc3Dframe = (transform().inverse()) * gpos;
  lpos                = Acts::Vector2D(loc3Dframe.perp(), loc3Dframe.phi());
  return ((std::abs(loc3Dframe.z()) > s_onSurfaceTolerance) ? false : true);
}

const Acts::Vector2D
Acts::DiscSurface::localPolarToLocalCartesian(const Vector2D& locpol) const
{
  const DiscTrapezoidalBounds* dtbo
      = dynamic_cast<const Acts::DiscTrapezoidalBounds*>(&(bounds()));
  if (dtbo) {
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
Acts::DiscSurface::localCartesianToGlobal(const Vector2D& lpos) const
{
  Vector3D loc3Dframe(lpos[Acts::eLOC_X], lpos[Acts::eLOC_Y], 0.);
  return Vector3D(transform() * loc3Dframe);
}

const Acts::Vector2D
Acts::DiscSurface::globalToLocalCartesian(const Vector3D& gpos, double) const
{
  Vector3D loc3Dframe = (transform().inverse()) * gpos;
  return Vector2D(loc3Dframe.x(), loc3Dframe.y());
}

std::string
Acts::DiscSurface::name() const
{
  return "Acts::DiscSurface";
}

bool
Acts::DiscSurface::isOnSurface(const Vector3D&      glopo,
                               const BoundaryCheck& bcheck) const
{
  Vector3D loc3Dframe = (transform().inverse()) * glopo;
  if (std::abs(loc3Dframe.z()) > (s_onSurfaceTolerance)) return false;
  return (bcheck
              ? bounds().inside(Vector2D(loc3Dframe.perp(), loc3Dframe.phi()),
                                bcheck)
              : true);
}

Acts::DiscSurface*
Acts::DiscSurface::clone(const Transform3D* shift) const
{
  if (shift) return new DiscSurface(*this, *shift);
  return new DiscSurface(*this);
}

const Acts::SurfaceBounds&
Acts::DiscSurface::bounds() const
{
  if (m_bounds) return (*(m_bounds.get()));
  return s_noBounds;
}

const Acts::Vector3D
Acts::DiscSurface::normal(const Acts::Vector2D&) const
{
  // fast access via tranform matrix (and not rotation())
  auto tMatrix = transform().matrix();
  return Vector3D(tMatrix(0, 2), tMatrix(1, 2), tMatrix(2, 2));
}

const Acts::Vector3D
    Acts::DiscSurface::binningPosition(Acts::BinningValue) const
{
  return center();
}

double
Acts::DiscSurface::pathCorrection(const Acts::Vector3D&,
                                  const Acts::Vector3D& mom) const
{
  /// we can ignore the global position here
  return 1. / std::abs(normal().dot(mom.unit()));
}

Acts::Intersection
Acts::DiscSurface::intersectionEstimate(const Acts::Vector3D&      gpos,
                                        const Acts::Vector3D&      dir,
                                        bool                       forceDir,
                                        const Acts::BoundaryCheck& bcheck) const
{
  double denom = dir.dot(normal());
  if (denom) {
    double   u = (normal().dot((center() - gpos))) / (denom);
    Vector3D intersectPoint(gpos + u * dir);
    // evaluate the intersection in terms of direction
    bool isValid = forceDir ? (u > 0.) : true;
    // evaluate (if necessary in terms of boundaries)
    isValid
        = bcheck ? (isValid && isOnSurface(intersectPoint, bcheck)) : isValid;
    // return the result
    return Intersection(intersectPoint, u, isValid);
  }
  return Intersection(gpos, 0., false);
}
