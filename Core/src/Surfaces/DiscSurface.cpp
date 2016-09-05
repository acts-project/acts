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
#include <iomanip>
#include <iostream>
#include "ACTS/Surfaces/DiscTrapezoidalBounds.hpp"
#include "ACTS/Surfaces/InfiniteBounds.hpp"
#include "ACTS/Surfaces/RadialBounds.hpp"
#include "ACTS/Utilities/Definitions.hpp"

Acts::DiscSurface::DiscSurface(const DiscSurface& dsf)
  : Surface(dsf), m_bounds(dsf.m_bounds)
{
}

Acts::DiscSurface::DiscSurface(const DiscSurface& dsf,
                               const Transform3D& transf)
  : Surface(dsf, transf), m_bounds(dsf.m_bounds)
{
}

Acts::DiscSurface::DiscSurface(std::shared_ptr<Transform3D> htrans,
                               double                       rmin,
                               double                       rmax,
                               double                       hphisec)
  : Surface(htrans)
  , m_bounds(std::make_shared<Acts::RadialBounds>(rmin, rmax, hphisec))
{
}

Acts::DiscSurface::DiscSurface(std::shared_ptr<Transform3D> htrans,
                               double                       minhalfx,
                               double                       maxhalfx,
                               double                       maxR,
                               double                       minR,
                               double                       avephi,
                               double                       stereo)
  : Surface(htrans)
  , m_bounds(std::make_shared<Acts::DiscTrapezoidalBounds>(minhalfx,
                                                           maxhalfx,
                                                           maxR,
                                                           minR,
                                                           avephi,
                                                           stereo))
{
}

Acts::DiscSurface::DiscSurface(std::shared_ptr<Transform3D>      htrans,
                               std::shared_ptr<const DiscBounds> dbounds)
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
Acts::DiscSurface::operator=(const DiscSurface& dsf)
{
  if (this != &dsf) {
    Acts::Surface::operator=(dsf);
    m_bounds               = dsf.m_bounds;
  }
  return *this;
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
  return ((fabs(loc3Dframe.z()) > s_onSurfaceTolerance) ? false : true);
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

bool
Acts::DiscSurface::isOnSurface(const Vector3D&      glopo,
                               const BoundaryCheck& bchk) const
{
  Vector3D loc3Dframe = (transform().inverse()) * glopo;
  if (fabs(loc3Dframe.z()) > (s_onSurfaceTolerance)) return false;
  return (
      bchk
          ? bounds().inside(Vector2D(loc3Dframe.perp(), loc3Dframe.phi()), bchk)
          : true);
}
