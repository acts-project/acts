// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// StraightLineSurface.cpp, ACTS project
///////////////////////////////////////////////////////////////////

// Geometry module
#include "ACTS/Surfaces/StraightLineSurface.hpp"
#include "ACTS/Surfaces/CylinderBounds.hpp"
// Core module
#include "ACTS/Utilities/Identifier.hpp"
// STD/STL
#include <iomanip>
#include <iostream>

Acts::NoBounds Acts::StraightLineSurface::s_boundless;

// default constructor
Acts::StraightLineSurface::StraightLineSurface()
  : Surface(), m_lineDirection(nullptr), m_bounds()
{
}

// constructors by arguments: boundless surface
Acts::StraightLineSurface::StraightLineSurface(
    std::shared_ptr<Acts::Transform3D> htrans)
  : Surface(htrans), m_lineDirection(nullptr), m_bounds()
{
}

// constructors by arguments: boundless surface
Acts::StraightLineSurface::StraightLineSurface(
    std::unique_ptr<Acts::Transform3D> htrans)
  : Surface(std::shared_ptr<Acts::Transform3D>(std::move(htrans)))
  , m_lineDirection(nullptr)
  , m_bounds()
{
}

// constructors by arguments
Acts::StraightLineSurface::StraightLineSurface(
    std::shared_ptr<Acts::Transform3D> htrans,
    double                             radius,
    double                             halez)
  : Surface(htrans)
  , m_lineDirection(nullptr)
  , m_bounds(std::make_shared<Acts::CylinderBounds>(radius, halez))
{
}

// constructors by arguments
Acts::StraightLineSurface::StraightLineSurface(
    std::shared_ptr<Acts::Transform3D>    htrans,
    std::shared_ptr<const CylinderBounds> cbounds)
  : Surface(htrans), m_lineDirection(nullptr), m_bounds(cbounds)
{
}

// constructor from detector elements
Acts::StraightLineSurface::StraightLineSurface(
    const Acts::DetectorElementBase& detelement,
    const Identifier&                id)
  : Surface(detelement, id), m_lineDirection(nullptr), m_bounds()
{
}

// copy constructor
Acts::StraightLineSurface::StraightLineSurface(
    const Acts::StraightLineSurface& slsf)
  : Surface(slsf), m_lineDirection(nullptr), m_bounds(slsf.m_bounds)
{
}

// copy constructor with shift
Acts::StraightLineSurface::StraightLineSurface(const StraightLineSurface& csf,
                                               const Acts::Transform3D& transf)
  : Surface(csf, transf), m_lineDirection(nullptr), m_bounds(csf.m_bounds)
{
}

// destructor (will call destructor from base class which deletes objects)
Acts::StraightLineSurface::~StraightLineSurface()
{
  delete m_lineDirection;
}

// assignment operator
Acts::StraightLineSurface&
Acts::StraightLineSurface::operator=(const Acts::StraightLineSurface& slsf)
{
  if (this != &slsf) {
    delete m_lineDirection;
    m_lineDirection        = 0;
    Acts::Surface::operator=(slsf);
    m_bounds               = slsf.m_bounds;
  }
  return *this;
}

bool
Acts::StraightLineSurface::operator==(const Acts::Surface& sf) const
{
  // first check the type not to compare apples with oranges
  const Acts::StraightLineSurface* slsf
      = dynamic_cast<const Acts::StraightLineSurface*>(&sf);
  if (!slsf) return false;
  bool transfEqual(transform().isApprox(slsf->transform(), 10e-8));
  bool centerEqual = (transfEqual) ? (center() == slsf->center()) : false;
  bool boundsEqual = (centerEqual) ? (bounds() == slsf->bounds()) : false;
  return boundsEqual;
}

// true local to global method - fully defined
void
Acts::StraightLineSurface::localToGlobal(const Acts::Vector2D& locpos,
                                         const Acts::Vector3D& glomom,
                                         Acts::Vector3D&       glopos) const
{
  // get the vector perpenticular to the momentum and the straw axis
  Acts::Vector3D radiusAxisGlobal(lineDirection().cross(glomom));
  Acts::Vector3D locZinGlobal
      = transform() * Acts::Vector3D(0., 0., locpos[Acts::eLOC_Z]);
  // transform zPosition into global coordinates and add Acts::eLOC_R *
  // radiusAxis
  glopos = Acts::Vector3D(
      locZinGlobal + locpos[Acts::eLOC_R] * radiusAxisGlobal.normalized());
}

// true global to local method - fully defined
bool
Acts::StraightLineSurface::globalToLocal(const Acts::Vector3D& glopos,
                                         const Acts::Vector3D& glomom,
                                         Acts::Vector2D&       locpos) const
{
  Acts::Vector3D loc3Dframe = (transform().inverse()) * glopos;
  // construct localPosition with sign*candidate.perp() and z.()
  locpos = Acts::Vector2D(loc3Dframe.perp(), loc3Dframe.z());
  Acts::Vector3D decVec(glopos - center());
  // assign the right sign
  double sign = ((lineDirection().cross(glomom)).dot(decVec) < 0.) ? -1. : 1.;
  locpos[Acts::eLOC_R] *= sign;
  return true;
}

// isOnSurface check
bool
Acts::StraightLineSurface::isOnSurface(const Acts::Vector3D& glopo,
                                       const BoundaryCheck&  bchk) const
{
  if (!bchk) return true;
  // check whether this is a boundless surface
  if (!(m_bounds.get()) && !Surface::m_associatedDetElement) return true;
  // get the standard bounds
  Acts::Vector3D loc3Dframe = (transform().inverse()) * glopo;
  Acts::Vector2D locCand(loc3Dframe.perp(), loc3Dframe.z());
  return (locCand[Acts::eLOC_R] < bounds().r() + bchk.toleranceLoc1
          && bounds().insideLoc2(locCand, bchk.toleranceLoc2));
}

// return the measurement frame
const Acts::RotationMatrix3D
Acts::StraightLineSurface::measurementFrame(const Acts::Vector3D&,
                                            const Acts::Vector3D& glomom) const
{
  Acts::RotationMatrix3D mFrame;
  // construct the measurement frame
  const Acts::Vector3D& measY = lineDirection();
  Acts::Vector3D        measX(measY.cross(glomom).unit());
  Acts::Vector3D        measDepth(measX.cross(measY));
  // assign the columnes
  mFrame.col(0) = measX;
  mFrame.col(1) = measY;
  mFrame.col(2) = measDepth;
  // return the rotation matrix
  return mFrame;
}

Acts::Intersection
Acts::StraightLineSurface::intersectionEstimate(const Acts::Vector3D& pos,
                                                const Acts::Vector3D& dir,
                                                bool                  forceDir,
                                                const BoundaryCheck& bchk) const
{
  // following nominclature found in header file and doxygen documentation
  // line one is the straight track
  const Vector3D& ma = pos;
  const Vector3D& ea = dir;
  // line two is the line surface
  const Vector3D& mb = center();
  const Vector3D& eb = lineDirection();
  // now go ahead and solve for the closest approach
  Vector3D mab(mb - ma);
  double   eaTeb = ea.dot(eb);
  double   denom = 1 - eaTeb * eaTeb;
  if (fabs(denom) > 10e-7) {
    double lambda0 = (mab.dot(ea) - mab.dot(eb) * eaTeb) / denom;
    // evaluate in terms of direction
    bool isValid = forceDir ? (lambda0 > 0.) : true;
    // evaluate validaty in terms of bounds
    Vector3D result = (ma + lambda0 * ea);
    isValid         = bchk ? (isValid && isOnSurface(result, bchk)) : isValid;
    // return the result
    return Acts::Intersection(result, lambda0, isValid);
  }
  return Acts::Intersection(pos, 0., false);
}
