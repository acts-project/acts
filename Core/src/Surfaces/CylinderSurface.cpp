// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// CylinderSurface.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Surfaces/CylinderSurface.hpp"

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>

#include "Acts/Utilities/ThrowAssert.hpp"
#include "Acts/Utilities/VariantData.hpp"
#include "Acts/Utilities/detail/RealQuadraticEquation.hpp"

Acts::CylinderSurface::CylinderSurface(const CylinderSurface& other)
  : GeometryObject(), Surface(other), m_bounds(other.m_bounds)
{
}

Acts::CylinderSurface::CylinderSurface(const CylinderSurface& other,
                                       const Transform3D&     transf)
  : GeometryObject(), Surface(other, transf), m_bounds(other.m_bounds)
{
}

Acts::CylinderSurface::CylinderSurface(
    std::shared_ptr<const Transform3D> htrans,
    double                             radius,
    double                             hlength)
  : GeometryObject()
  , Surface(htrans)
  , m_bounds(std::make_shared<const CylinderBounds>(radius, hlength))
{
}

Acts::CylinderSurface::CylinderSurface(
    std::shared_ptr<const Transform3D> htrans,
    double                             radius,
    double                             hphi,
    double                             hlength)
  : GeometryObject()
  , Surface(htrans)
  , m_bounds(std::make_shared<const CylinderBounds>(radius, hphi, hlength))
{
}

Acts::CylinderSurface::CylinderSurface(
    std::shared_ptr<const CylinderBounds> cbounds,
    const Acts::DetectorElementBase&      detelement,
    const Identifier&                     identifier)
  : Surface(detelement, identifier), m_bounds(cbounds)
{
  /// surfaces representing a detector element must have bounds
  assert(cbounds);
}

Acts::CylinderSurface::CylinderSurface(
    std::shared_ptr<const Transform3D>    htrans,
    std::shared_ptr<const CylinderBounds> cbounds)
  : Surface(htrans), m_bounds(cbounds)
{
  throw_assert(cbounds, "CylinderBounds must not be nullptr");
}

Acts::CylinderSurface::CylinderSurface(const variant_data& data_)
{
  throw_assert(data_.which() == 4, "Variant data must be map");
  variant_map data = boost::get<variant_map>(data_);
  throw_assert(data.count("type"), "Variant data must have type.");
  std::string type = data.get<std::string>("type");
  throw_assert(type == "CylinderSurface",
               "Variant data type must be CylinderSurface");

  variant_map payload    = data.get<variant_map>("payload");
  variant_map var_bounds = payload.get<variant_map>("bounds");

  m_bounds = std::make_shared<const CylinderBounds>(var_bounds);

  if (payload.count("transform")) {
    // we have a transform
    auto trf = std::make_shared<const Transform3D>(
        from_variant<Transform3D>(payload.get<variant_map>("transform")));
    m_transform = trf;
  }
}

Acts::CylinderSurface::~CylinderSurface()
{
}

Acts::CylinderSurface&
Acts::CylinderSurface::operator=(const CylinderSurface& other)
{
  if (this != &other) {
    Surface::operator=(other);
    m_bounds         = other.m_bounds;
  }
  return *this;
}

// return the binning position for ordering in the BinnedArray
const Acts::Vector3D
Acts::CylinderSurface::binningPosition(BinningValue bValue) const
{
  // special binning type for R-type methods
  if (bValue == Acts::binR || bValue == Acts::binRPhi) {
    double R   = bounds().r();
    double phi = m_bounds ? m_bounds->averagePhi() : 0.;
    return Vector3D(
        center().x() + R * cos(phi), center().y() + R * sin(phi), center().z());
  }
  // give the center as default for all of these binning types
  // binX, binY, binZ, binR, binPhi, binRPhi, binH, binEta
  return center();
}

// return the measurement frame: it's the tangential plane
const Acts::RotationMatrix3D
Acts::CylinderSurface::referenceFrame(const Vector3D& gpos,
                                      const Vector3D&) const
{
  RotationMatrix3D mFrame;
  // construct the measurement frame
  // measured Y is the z axis
  Vector3D measY = rotSymmetryAxis();
  // measured z is the position normalized transverse (in local)
  Vector3D measDepth = normal(gpos);
  // measured X is what comoes out of it
  Vector3D measX(measY.cross(measDepth).unit());
  // assign the columnes
  mFrame.col(0) = measX;
  mFrame.col(1) = measY;
  mFrame.col(2) = measDepth;
  // return the rotation matrix
  return mFrame;
}

Acts::Surface::SurfaceType
Acts::CylinderSurface::type() const
{
  return Surface::Cylinder;
}

void
Acts::CylinderSurface::localToGlobal(const Vector2D& lpos,
                                     const Vector3D&,
                                     Vector3D& gpos) const
{
  // create the position in the local 3d frame
  double r   = bounds().r();
  double phi = lpos[Acts::eLOC_RPHI] / r;
  gpos       = Vector3D(r * cos(phi), r * sin(phi), lpos[Acts::eLOC_Z]);
  // transform it to the globalframe: CylinderSurfaces are allowed to have 0
  // if pointer transform exists -> port into frame
  if (Surface::m_transform) gpos = transform() * gpos;
}

bool
Acts::CylinderSurface::globalToLocal(const Vector3D& gpos,
                                     const Vector3D&,
                                     Vector2D& lpos) const
{
  // get the transform & transform global position into cylinder frame
  // @todo clean up intolerance parameters
  // transform it to the globalframe: CylinderSurfaces are allowed to have 0
  // pointer transform
  double radius             = 0.;
  double inttol             = bounds().r() * 0.0001;
  if (inttol < 0.01) inttol = 0.01;
  // do the transformation or not
  if (Surface::m_transform) {
    const Transform3D& surfaceTrans = transform();
    Transform3D        inverseTrans(surfaceTrans.inverse());
    Vector3D           loc3Dframe(inverseTrans * gpos);
    lpos   = Vector2D(bounds().r() * loc3Dframe.phi(), loc3Dframe.z());
    radius = loc3Dframe.perp();
  } else {
    lpos   = Vector2D(bounds().r() * gpos.phi(), gpos.z());
    radius = gpos.perp();
  }
  // return true or false
  return ((std::abs(radius - bounds().r()) > inttol) ? false : true);
}

bool
Acts::CylinderSurface::isOnSurface(const Acts::Vector3D& gpos,
                                   const BoundaryCheck&  bcheck) const
{
  Acts::Vector3D loc3Dframe
      = Acts::Surface::m_transform ? (transform().inverse()) * gpos : gpos;
  return (bcheck ? bounds().inside3D(loc3Dframe, bcheck) : true);
}

Acts::Intersection
Acts::CylinderSurface::intersectionEstimate(const Acts::Vector3D& gpos,
                                            const Acts::Vector3D& gdir,
                                            bool                  forceDir,
                                            const BoundaryCheck&  bcheck) const
{
  bool needsTransform
      = (Acts::Surface::m_transform || m_associatedDetElement) ? true : false;
  // create the hep points
  Vector3D point1    = gpos;
  Vector3D direction = gdir;
  if (needsTransform) {
    Transform3D invTrans = transform().inverse();
    point1               = invTrans * gpos;
    direction            = invTrans.linear() * gdir;
  }
  Acts::Vector3D point2 = point1 + direction;
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
    detail::RealQuadraticEquation pquad(1 + k * k, 2 * k * d, d * d - R * R);
    if (pquad.solutions == 2) {
      // the solutions in the 3D frame of the cylinder
      t1 = (pquad.first - point1.x()) / direction.x();
      t2 = (pquad.second - point1.x()) / direction.x();
    } else {  // bail out if no solution exists
      return Intersection(gpos, 0., false);
    }
  } else {
    // bail out if no solution exists
    if (!direction.y()) return Intersection(gpos, 0., false);
    // x value ise th one of point1
    // x^2 + y^2 = R^2
    // y = sqrt(R^2-x^2)
    double x     = point1.x();
    double r2mx2 = R * R - x * x;
    // bail out if no solution
    if (r2mx2 < 0.) return Intersection(gpos, 0., false);
    double y = sqrt(r2mx2);
    // assign parameters and solutions
    t1 = (y - point1.y()) / direction.y();
    t2 = (-y - point1.y()) / direction.y();
  }
  Vector3D sol1raw(point1 + t1 * direction);
  Vector3D sol2raw(point1 + t2 * direction);
  // now reorder and return
  Vector3D solution(0, 0, 0);
  double   path = 0.;

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
  isValid = bcheck ? (isValid && bounds().inside3D(solution, bcheck)) : isValid;

  // now return
  return needsTransform ? Intersection((transform() * solution), path, isValid)
                        : Intersection(solution, path, isValid);
}

std::string
Acts::CylinderSurface::name() const
{
  return "Acts::CylinderSurface";
}

Acts::CylinderSurface*
Acts::CylinderSurface::clone(const Acts::Transform3D* shift) const
{
  if (shift) return new CylinderSurface(*this, *shift);
  return new CylinderSurface(*this);
}

const Acts::Vector3D
Acts::CylinderSurface::normal(const Acts::Vector2D& lpos) const
{
  double   phi = lpos[Acts::eLOC_RPHI] / m_bounds->r();
  Vector3D localNormal(cos(phi), sin(phi), 0.);
  return Vector3D(transform().matrix().block<3, 3>(0, 0) * localNormal);
}

const Acts::Vector3D
Acts::CylinderSurface::normal(const Acts::Vector3D& gpos) const
{
  // get it into the cylinder frame if needed
  Vector3D pos3D          = gpos;
  bool     needsTransform = (m_transform || m_associatedDetElement);
  if (needsTransform) {
    pos3D = transform().inverse() * gpos;
  }
  // set the z coordinate to 0
  pos3D.z() = 0.;
  // normalize and rotate back into global if needed
  return needsTransform ? transform().linear() * pos3D.unit() : pos3D.unit();
}

double
Acts::CylinderSurface::pathCorrection(const Acts::Vector3D& gpos,
                                      const Acts::Vector3D& mom) const
{
  Vector3D normalT  = normal(gpos);
  double   cosAlpha = normalT.dot(mom.unit());
  return std::fabs(1. / cosAlpha);
}

const Acts::CylinderBounds&
Acts::CylinderSurface::bounds() const
{
  return (*m_bounds.get());
}

Acts::variant_data
Acts::CylinderSurface::toVariantData() const
{
  using namespace std::string_literals;

  variant_map payload;
  payload["bounds"] = m_bounds->toVariantData();

  if (m_transform) {
    payload["transform"] = to_variant(*m_transform);
  }

  variant_map data;
  data["type"]    = "CylinderSurface"s;
  data["payload"] = payload;
  return data;
}
