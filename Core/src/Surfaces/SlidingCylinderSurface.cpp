// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// SlidingCylinderSurface.cpp, ACTS project
///////////////////////////////////////////////////////////////////

// Geometry module
#include "ACTS/Surfaces/SlidingCylinderSurface.hpp"
#include "ACTS/Surfaces/CylinderBounds.hpp"

// default constructor
Acts::SlidingCylinderSurface::SlidingCylinderSurface()
  : Acts::CylinderSurface(), m_depth(), m_etaBin(), m_align()
{
}

// copy constructor
Acts::SlidingCylinderSurface::SlidingCylinderSurface(
    const SlidingCylinderSurface& dsf)
  : Acts::CylinderSurface(dsf)
  , m_depth(new std::vector<float>(*(dsf.m_depth)))
  , m_etaBin(dsf.m_etaBin->clone())
  , m_align(dsf.m_align ? new Acts::Transform3D(*dsf.m_align) : nullptr)
{
}

// copy constructor with shift
Acts::SlidingCylinderSurface::SlidingCylinderSurface(
    const SlidingCylinderSurface& dsf,
    const Acts::Transform3D&      transf)
  : Acts::CylinderSurface(dsf, transf)
  , m_depth(new std::vector<float>(*(dsf.m_depth)))
  , m_etaBin(dsf.m_etaBin->clone())
  , m_align(dsf.m_align ? new Acts::Transform3D(*dsf.m_align) : nullptr)
{
}

// constructor
Acts::SlidingCylinderSurface::SlidingCylinderSurface(
    const Acts::CylinderSurface& dsf,
    Acts::BinUtility*            bu,
    const std::vector<float>*    offset,
    Acts::Transform3D*           align)
  : Acts::CylinderSurface(dsf), m_depth(offset), m_etaBin(bu), m_align(align)
{
}

// destructor (will call destructor from base class which deletes objects)
Acts::SlidingCylinderSurface::~SlidingCylinderSurface()
{
  delete m_depth;
  delete m_etaBin;
  delete m_align;
}

Acts::SlidingCylinderSurface&
Acts::SlidingCylinderSurface::operator=(const SlidingCylinderSurface& dsf)
{
  if (this != &dsf) {
    Acts::CylinderSurface::operator=(dsf);
    delete m_depth;
    m_depth = new std::vector<float>(*(dsf.m_depth));
    delete m_etaBin;
    m_etaBin = dsf.m_etaBin->clone();
    delete m_align;
    m_align = (dsf.m_align ? new Acts::Transform3D(*dsf.m_align) : nullptr);
  }
  return *this;
}

bool
Acts::SlidingCylinderSurface::operator==(const Acts::Surface& sf) const
{
  // first check the type not to compare apples with oranges
  const Acts::SlidingCylinderSurface* dsf
      = dynamic_cast<const Acts::SlidingCylinderSurface*>(&sf);
  if (!dsf) return false;
  if (this == dsf) return true;
  bool transfEqual(transform().isApprox(dsf->transform(), 10e-8));
  bool centerEqual = (transfEqual) ? (center() == dsf->center()) : false;
  bool boundsEqual = (centerEqual) ? (bounds() == dsf->bounds()) : false;
  return boundsEqual;
}

void
Acts::SlidingCylinderSurface::localToGlobal(const Acts::Vector2D& locpos,
                                            const Acts::Vector3D&,
                                            Acts::Vector3D& glopos) const
{
  // create the position in the local 3d frame
  double         r0   = bounds().r();
  double         phi0 = locpos[Acts::eLOC_RPHI] / r0;
  Acts::Vector3D loc3D0(r0 * cos(phi0), r0 * sin(phi0), locpos[Acts::eLOC_Z]);
  // correct for alignment, retrieve offset correction
  Acts::Transform3D t0
      = m_align ? m_align->inverse() * transform() : transform();
  float          offset = m_depth ? (*m_depth)[m_etaBin->bin(t0 * loc3D0)] : 0.;
  double         r      = r0 + offset;
  double         phi    = locpos[Acts::eLOC_RPHI] / r;
  Acts::Vector3D loc3Dframe(r * cos(phi), r * sin(phi), locpos[Acts::eLOC_Z]);
  // transport it to the globalframe
  glopos = Acts::Surface::transform() * loc3Dframe;
}

/** local<->global transformation in case of polar local coordinates */
bool
Acts::SlidingCylinderSurface::globalToLocal(const Acts::Vector3D& glopos,
                                            const Acts::Vector3D&,
                                            Acts::Vector2D& locpos) const
{
  // get the transform & transform global position into cylinder frame
  // transform it to the globalframe: CylinderSurfaces are allowed to have 0
  // pointer transform
  double radius             = 0.;
  double inttol             = bounds().r() * 0.0001;
  if (inttol < 0.01) inttol = 0.01;
  // realign to find local eta bin
  Acts::Vector3D loc3D0 = m_align ? (m_align->inverse() * glopos) : glopos;
  float          offset = (*m_depth)[m_etaBin->bin(loc3D0)];
  // do the transformation or not
  if (Acts::Surface::m_transform) {
    const Acts::Transform3D& surfaceTrans = transform();
    Acts::Transform3D        inverseTrans(surfaceTrans.inverse());
    Acts::Vector3D           loc3Dframe(inverseTrans * glopos);
    locpos = Acts::Vector2D((bounds().r() + offset) * loc3Dframe.phi(),
                            loc3Dframe.z());
    radius = loc3Dframe.perp();
  } else {
    locpos = Acts::Vector2D((bounds().r() + offset) * glopos.phi(), glopos.z());
    radius = glopos.perp();
  }
  // return true or false
  return ((fabs(radius - bounds().r() - offset) > inttol) ? false : true);
}

bool
Acts::SlidingCylinderSurface::isOnSurface(const Acts::Vector3D&      glopo,
                                          const Acts::BoundaryCheck& bchk) const
{
  Acts::Vector3D loc3D0 = m_align ? m_align->inverse() * glopo : glopo;
  Acts::Vector3D loc3Dframe
      = m_transform ? (transform().inverse()) * glopo : glopo;
  float offset = (*m_depth)[m_etaBin->bin(loc3D0)];
  // recalculate r to match bounds
  Acts::Vector3D loc3Dbase((loc3Dframe.perp() - offset) * cos(loc3Dframe.phi()),
                           (loc3Dframe.perp() - offset) * sin(loc3Dframe.phi()),
                           loc3Dframe.z());

  return (bchk ? bounds().inside3D(loc3Dbase, bchk) : true);
}
