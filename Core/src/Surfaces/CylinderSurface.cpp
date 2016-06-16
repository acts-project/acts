// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// CylinderSurface.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Surfaces/CylinderSurface.hpp"
#include "ACTS/Surfaces/RealQuadraticEquation.hpp"
#include <assert.h>
#include <iomanip>
#include <iostream>

Acts::CylinderSurface::CylinderSurface(const CylinderSurface& csf)
  : Acts::Surface(csf), m_bounds(csf.m_bounds), m_rotSymmetryAxis(nullptr)
{
}

Acts::CylinderSurface::CylinderSurface(const CylinderSurface&   csf,
                                       const Acts::Transform3D& transf)
  : Acts::Surface(csf, transf)
  , m_bounds(csf.m_bounds)
{
}

Acts::CylinderSurface::CylinderSurface(
    std::shared_ptr<Acts::Transform3D> htrans,
    double                             radius,
    double                             hlength)
  : Acts::Surface(htrans)
  , m_bounds(std::make_shared<Acts::CylinderBounds>(radius, hlength))
{
}

Acts::CylinderSurface::CylinderSurface(
    std::shared_ptr<Acts::Transform3D> htrans,
    double                             radius,
    double                             hphi,
    double                             hlength)
  : Acts::Surface(htrans)
  , m_bounds(std::make_shared<Acts::CylinderBounds>(radius, hphi, hlength))
{
}

Acts::CylinderSurface::CylinderSurface(
    std::shared_ptr<Acts::Transform3D>          htrans,
    std::shared_ptr<const Acts::CylinderBounds> cbounds)
  : Acts::Surface(htrans)
  , m_bounds(cbounds)
{
  assert(cbounds);
}

Acts::CylinderSurface::~CylinderSurface()
{  
}

Acts::CylinderSurface&
Acts::CylinderSurface::operator=(const CylinderSurface& csf)
{
  if (this != &csf) {
    Acts::Surface::operator=(csf);
    m_bounds               = csf.m_bounds;
  }
  return *this;
}

// return the binning position for ordering in the BinnedArray
Acts::Vector3D
Acts::CylinderSurface::binningPosition(Acts::BinningValue bValue) const
{
  // special binning type for R-type methods
  if (bValue == Acts::binR || bValue == Acts::binRPhi) {
    double R   = bounds().r();
    double phi = m_bounds ? m_bounds->averagePhi() : 0.;
    return Acts::Vector3D(
        center().x() + R * cos(phi), center().y() + R * sin(phi), center().z());
  }
  // give the center as default for all of these binning types
  // binX, binY, binZ, binR, binPhi, binRPhi, binH, binEta
  return Acts::Surface::binningPosition(bValue);
}

// return the measurement frame: it's the tangential plane
const Acts::RotationMatrix3D
Acts::CylinderSurface::measurementFrame(const Acts::Vector3D& gpos,
                                        const Acts::Vector3D&) const
{
  Acts::RotationMatrix3D mFrame;
  // construct the measurement frame
  Acts::Vector3D measY(
      transform().rotation().col(2));  // measured Y is the z axis
  Acts::Vector3D measDepth
      = Acts::Vector3D(gpos.x(), gpos.y(), 0.)
            .unit();  // measured z is the position transverse normalized
  Acts::Vector3D measX(
      measY.cross(measDepth).unit());  // measured X is what comoes out of it
  // the columnes
  mFrame.col(0) = measX;
  mFrame.col(1) = measY;
  mFrame.col(2) = measDepth;
  // return the rotation matrix
  return std::move(mFrame);
}

const Acts::Vector3D
Acts::CylinderSurface::rotSymmetryAxis() const
{
  return std::move(Vector3D(transform().rotation().col(2)));
}

void
Acts::CylinderSurface::localToGlobal(const Acts::Vector2D& lpos,
                                     const Acts::Vector3D&,
                                     Acts::Vector3D& gpos) const
{
  // create the position in the local 3d frame
  double r   = bounds().r();
  double phi = lpos[Acts::eLOC_RPHI] / r;
  gpos     = Acts::Vector3D(r * cos(phi), r * sin(phi), lpos[Acts::eLOC_Z]);
  // transform it to the globalframe: CylinderSurfaces are allowed to have 0
  // pointer transform
  if (Acts::Surface::m_transform) gpos = transform() * gpos;
}

bool
Acts::CylinderSurface::globalToLocal(const Acts::Vector3D& gpos,
                                     const Acts::Vector3D&,
                                     Acts::Vector2D& lpos) const
{
  // get the transform & transform global position into cylinder frame
  // @TODO clean up intolerance parameters
  // transform it to the globalframe: CylinderSurfaces are allowed to have 0
  // pointer transform
  double radius             = 0.;
  double inttol             = bounds().r() * 0.0001;
  if (inttol < 0.01) inttol = 0.01;
  // do the transformation or not
  if (Acts::Surface::m_transform) {
    const Acts::Transform3D& surfaceTrans = transform();
    Acts::Transform3D        inverseTrans(surfaceTrans.inverse());
    Acts::Vector3D           loc3Dframe(inverseTrans * gpos);
    lpos = Acts::Vector2D(bounds().r() * loc3Dframe.phi(), loc3Dframe.z());
    radius = loc3Dframe.perp();
  } else {
    lpos = Acts::Vector2D(bounds().r() * gpos.phi(), gpos.z());
    radius = gpos.perp();
  }
  // return true or false
  return ((fabs(radius - bounds().r()) > inttol) ? false : true);
}

bool
Acts::CylinderSurface::isOnSurface(const Acts::Vector3D& gpos,
                                   const BoundaryCheck&  bchk) const
{
  Acts::Vector3D loc3Dframe
      = Acts::Surface::m_transform ? (transform().inverse()) * gpos : gpos;
  return (bchk ? bounds().inside3D(loc3Dframe, bchk) : true);
}

Acts::Intersection
Acts::CylinderSurface::intersectionEstimate(const Acts::Vector3D& gpos,
                                            const Acts::Vector3D& dir,
                                            bool                  forceDir,
                                            const BoundaryCheck&  bchk) const
{
  bool needsTransform = (m_transform || m_associatedDetElement) ? true : false;
  // create the hep points
  Acts::Vector3D point1    = gpos;
  Acts::Vector3D direction = dir;
  if (needsTransform) {
    Acts::Transform3D invTrans = transform().inverse();
    point1                     = invTrans * gpos;
    direction                  = invTrans.linear() * dir;
  }
  Acts::Vector3D point2 = point1 + dir;
  // the bounds radius
  double R  = bounds().r();
  double t1 = 0.;
  double t2 = 0.;
  if (direction.x()) {
    // get line and circle constants
    double k = (direction.y()) / (direction.x());
    double d = (point2.x() * point1.y() - point1.x() * point2.y())
        / (point2.x() - point1.x());
    // and solve the qaudratic equation
    Acts::RealQuadraticEquation pquad(1 + k * k, 2 * k * d, d * d - R * R);
    if (pquad.solutions != Acts::none) {
      // the solutions in the 3D frame of the cylinder
      t1 = (pquad.first - point1.x()) / direction.x();
      t2 = (pquad.second - point1.x()) / direction.x();
    } else  // bail out if no solution exists
      return std::move(Intersection(gpos, 0., false));
  } else {
    // x value ise th one of point1
    // x^2 + y^2 = R^2
    // y = sqrt(R^2-x^2)
    double x     = point1.x();
    double r2mx2 = R * R - x * x;
    // bail out if no solution
    if (r2mx2 < 0.) return std::move(Intersection(gpos, 0., false));
    double y = sqrt(r2mx2);
    // assign parameters and solutions
    t1 = y - point1.y();
    t2 = -y - point1.y();
  }
  Acts::Vector3D sol1raw(point1 + t1 * direction);
  Acts::Vector3D sol2raw(point1 + t2 * direction);
  // now reorder and return
  Acts::Vector3D solution(0, 0, 0);
  double         path = 0.;

  // first check the validity of the direction
  bool isValid = true;

  // both solutions are of same sign, take the smaller, but flag as false if not
  // forward
  if (t1 * t2 > 0 || !forceDir) {
    // asign validity
    isValid = forceDir ? (t1 > 0.) : true;
    // assign the right solution
    if (t1 * t1 < t2 * t2) {
      solution = sol1raw;
      path     = t1;
    } else {
      solution = sol2raw;
      path     = t2;
    }
  } else {
    if (t1 > 0.) {
      solution = sol1raw;
      path     = t1;
    } else {
      solution = sol2raw;
      path     = t2;
    }
  }
  // the solution is still in the local 3D frame, direct check
  isValid = bchk ? (isValid && m_bounds->inside3D(solution, bchk)) : isValid;

  // now return
  return needsTransform ? std::move(Intersection(transform() * solution, path, isValid))
                        : std::move(Intersection(solution, path, isValid));
}
