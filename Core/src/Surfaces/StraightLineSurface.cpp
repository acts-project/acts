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

#include "ACTS/Surfaces/StraightLineSurface.hpp"
#include "ACTS/Utilities/Identifier.hpp"
#include <iomanip>
#include <iostream>

Acts::StraightLineSurface::StraightLineSurface(
    std::shared_ptr<Transform3D> htrans,
    double                             radius,
    double                             halez)
  : Surface(htrans)
  , m_bounds(std::make_shared<LineBounds>(radius, halez))
{
}

Acts::StraightLineSurface::StraightLineSurface(
    std::shared_ptr<Acts::Transform3D>    htrans,
    std::shared_ptr<const LineBounds> lbounds)
  : Surface(htrans), m_bounds(lbounds)
{
}

// constructor from detector elements
Acts::StraightLineSurface::StraightLineSurface(
    std::shared_ptr<const LineBounds> lbounds,
    const Acts::DetectorElementBase& detelement,
    const Identifier&                id)
  : Surface(detelement, id), m_bounds(lbounds)
{
}

// copy constructor
Acts::StraightLineSurface::StraightLineSurface(
    const Acts::StraightLineSurface& slsf)
  : Surface(slsf), m_bounds(slsf.m_bounds)
{
}

// copy constructor with shift
Acts::StraightLineSurface::StraightLineSurface(const StraightLineSurface& csf,
                                               const Acts::Transform3D& transf)
  : Surface(csf, transf), m_bounds(csf.m_bounds)
{
}

// destructor (will call destructor from base class which deletes objects)
Acts::StraightLineSurface::~StraightLineSurface()
{
}

// assignment operator
Acts::StraightLineSurface&
Acts::StraightLineSurface::operator=(const Acts::StraightLineSurface& slsf)
{
  if (this != &slsf) {
    Acts::Surface::operator=(slsf);
    m_bounds               = slsf.m_bounds;
  }
  return *this;
}

void
Acts::StraightLineSurface::localToGlobal(const Acts::Vector2D& lpos,
                                         const Acts::Vector3D& mom,
                                         Acts::Vector3D&       gpos) const
{
  // get the vector perpenticular to the momentum and the straw axis
  Acts::Vector3D radiusAxisGlobal(lineDirection().cross(mom));
  Acts::Vector3D locZinGlobal
      = transform() * Acts::Vector3D(0., 0., lpos[Acts::eLOC_Z]);
  // transform zPosition into global coordinates and add Acts::eLOC_R *
  // radiusAxis
  gpos = Acts::Vector3D(
      locZinGlobal + lpos[Acts::eLOC_R] * radiusAxisGlobal.normalized());
}

// true global to local method - fully defined
bool
Acts::StraightLineSurface::globalToLocal(const Acts::Vector3D& gpos,
                                         const Acts::Vector3D& mom,
                                         Acts::Vector2D&       lpos) const
{
  Acts::Vector3D loc3Dframe = (transform().inverse()) * gpos;
  // construct localPosition with sign*candidate.perp() and z.()
  lpos = Acts::Vector2D(loc3Dframe.perp(), loc3Dframe.z());
  Acts::Vector3D decVec(gpos - center());
  // assign the right sign
  double sign = ((lineDirection().cross(mom)).dot(decVec) < 0.) ? -1. : 1.;
  lpos[Acts::eLOC_R] *= sign;
  return true;
}

// isOnSurface check
bool
Acts::StraightLineSurface::isOnSurface(const Acts::Vector3D& gpos,
                                       const BoundaryCheck&  bchk) const
{
  if (!bchk) return true;
  // check whether this is a boundless surface
  if (!(m_bounds.get()) && !Surface::m_associatedDetElement) return true;
  // get the standard bounds
  Acts::Vector3D loc3Dframe = (transform().inverse()) * gpos;
  Acts::Vector2D locCand(loc3Dframe.perp(), loc3Dframe.z());
  return (locCand[Acts::eLOC_R] < bounds().r() + bchk.toleranceLoc0
          && bounds().insideLoc1(locCand, bchk.toleranceLoc1));
}

// return the measurement frame
const Acts::RotationMatrix3D
Acts::StraightLineSurface::measurementFrame(const Acts::Vector3D&,
                                            const Acts::Vector3D& mom) const
{
  Acts::RotationMatrix3D mFrame;
  // construct the measurement frame
  const Acts::Vector3D& measY = lineDirection();
  Acts::Vector3D        measX(measY.cross(mom).unit());
  Acts::Vector3D        measDepth(measX.cross(measY));
  // assign the columnes
  mFrame.col(0) = measX;
  mFrame.col(1) = measY;
  mFrame.col(2) = measDepth;
  // return the rotation matrix
  return std::move(mFrame);
}

Acts::Intersection
Acts::StraightLineSurface::intersectionEstimate(const Vector3D&      gpos,
                                                const Vector3D&      dir,
                                                bool                 forceDir,
                                                const BoundaryCheck& bchk) const
{
  // following nominclature found in header file and doxygen documentation
  // line one is the straight track
  const Vector3D& ma = gpos;
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
    return std::move(Intersection(result, lambda0, isValid));
  }
  return std::move(Intersection(gpos, 0., false));
}
