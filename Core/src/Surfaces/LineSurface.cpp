// This file is part of the ACTS project.
//
// Copyright (C) 2016-2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// LineSurface.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "Acts/Surfaces/LineSurface.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>

#include "Acts/Utilities/ThrowAssert.hpp"
#include "Acts/Utilities/VariantData.hpp"

Acts::LineSurface::LineSurface(std::shared_ptr<const Transform3D> htrans,
                               double                             radius,
                               double                             halez)
  : GeometryObject()
  , Surface(htrans)
  , m_bounds(std::make_shared<const LineBounds>(radius, halez))
{
}

Acts::LineSurface::LineSurface(std::shared_ptr<const Transform3D> htrans,
                               std::shared_ptr<const LineBounds>  lbounds)
  : GeometryObject(), Surface(htrans), m_bounds(lbounds)
{
}

Acts::LineSurface::LineSurface(std::shared_ptr<const LineBounds> lbounds,
                               const DetectorElementBase&        detelement,
                               const Identifier&                 id)
  : GeometryObject(), Surface(detelement, id), m_bounds(lbounds)
{
  throw_assert(lbounds, "LineBounds must not be nullptr");
}

Acts::LineSurface::LineSurface(const LineSurface& other)
  : GeometryObject(), Surface(other), m_bounds(other.m_bounds)
{
}

Acts::LineSurface::LineSurface(const LineSurface& other,
                               const Transform3D& transf)
  : GeometryObject(), Surface(other, transf), m_bounds(other.m_bounds)
{
}

Acts::LineSurface::LineSurface(const variant_data& data_) : GeometryObject()
{
  throw_assert(data_.which() == 4, "Variant data must be map");
  variant_map data = boost::get<variant_map>(data_);
  throw_assert(data.count("type"), "Variant data must have type.");
  // std::string type = boost::get<std::string>(data["type"]);
  std::string type = data.get<std::string>("type");
  throw_assert(type == "LineSurface", "Variant data type must be LineSurface");

  variant_map payload    = data.get<variant_map>("payload");
  variant_map bounds     = payload.get<variant_map>("bounds");
  std::string boundsType = bounds.get<std::string>("type");

  throw_assert(boundsType == "LineBounds",
               "Can only construct LineSurface from LineBounds");

  m_bounds = std::make_shared<const LineBounds>(bounds);

  if (payload.count("transform")) {
    // we have a transform
    auto trf = std::make_shared<const Transform3D>(
        from_variant<Transform3D>(payload.get<variant_map>("transform")));
    m_transform = trf;
  }
}

Acts::LineSurface::~LineSurface()
{
}

Acts::LineSurface&
Acts::LineSurface::operator=(const LineSurface& other)
{
  if (this != &other) {
    Surface::operator=(other);
    m_bounds         = other.m_bounds;
  }
  return *this;
}

void
Acts::LineSurface::localToGlobal(const Vector2D& lpos,
                                 const Vector3D& mom,
                                 Vector3D&       gpos) const
{
  // get the vector perpendicular to the momentum and the straw axis
  Vector3D radiusAxisGlobal(lineDirection().cross(mom));
  Vector3D locZinGlobal(0., 0., lpos[Acts::eLOC_Z]);
  // apply a transform if needed
  if (m_transform || m_associatedDetElement)
    locZinGlobal = transform() * locZinGlobal;
  // transform zPosition into global coordinates and
  // add Acts::eLOC_R * radiusAxis
  gpos = Vector3D(locZinGlobal
                  + lpos[Acts::eLOC_R] * radiusAxisGlobal.normalized());
}

bool
Acts::LineSurface::globalToLocal(const Vector3D& gpos,
                                 const Vector3D& mom,
                                 Vector2D&       lpos) const
{
  // apply the transform when needed
  Acts::Vector3D loc3Dframe = (m_transform || m_associatedDetElement)
      ? (transform().inverse()) * gpos
      : gpos;
  // construct localPosition with sign*candidate.perp() and z.()
  lpos = Acts::Vector2D(loc3Dframe.perp(), loc3Dframe.z());
  Acts::Vector3D decVec(gpos - center());
  // assign the right sign
  double sign = ((lineDirection().cross(mom)).dot(decVec) < 0.) ? -1. : 1.;
  lpos[Acts::eLOC_R] *= sign;
  return true;
}

bool
Acts::LineSurface::isOnSurface(const Vector3D&      gpos,
                               const BoundaryCheck& bcheck) const
{
  if (!bcheck) return true;
  // check whether this is a boundless surface
  if (!m_bounds && !Surface::m_associatedDetElement) return true;
  // get the standard bounds
  Vector3D loc3Dframe = (transform().inverse()) * gpos;
  Vector2D locCand(loc3Dframe.perp(), loc3Dframe.z());
  return bounds().inside(locCand, bcheck);
}

std::string
Acts::LineSurface::name() const
{
  return "Acts::LineSurface";
}

const Acts::RotationMatrix3D
Acts::LineSurface::referenceFrame(const Vector3D&, const Vector3D& mom) const
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
  return mFrame;
}

Acts::Intersection
Acts::LineSurface::intersectionEstimate(const Vector3D&      gpos,
                                        const Vector3D&      gdir,
                                        bool                 forceDir,
                                        const BoundaryCheck& bcheck) const
{
  // following nominclature found in header file and doxygen documentation
  // line one is the straight track
  const Vector3D& ma = gpos;
  const Vector3D& ea = gdir;
  // line two is the line surface
  const Vector3D& mb = center();
  const Vector3D  eb = lineDirection();
  // now go ahead and solve for the closest approach
  Vector3D mab(mb - ma);
  double   eaTeb = ea.dot(eb);
  double   denom = 1 - eaTeb * eaTeb;
  if (std::abs(denom) > 10e-7) {
    double lambda0 = (mab.dot(ea) - mab.dot(eb) * eaTeb) / denom;
    // evaluate in terms of direction
    bool isValid = forceDir ? (lambda0 > 0.) : true;
    // evaluate validaty in terms of bounds
    Vector3D result = (ma + lambda0 * ea);
    isValid = bcheck ? (isValid && isOnSurface(result, bcheck)) : isValid;
    // return the result
    return Intersection(result, lambda0, isValid);
  }
  return Intersection(gpos, std::numeric_limits<double>::max(), false);
}

double
Acts::LineSurface::pathCorrection(const Acts::Vector3D&,
                                  const Acts::Vector3D&) const
{
  return 1.;
}

const Acts::Vector3D
    Acts::LineSurface::binningPosition(Acts::BinningValue /*bValue*/) const
{
  return center();
}

const Acts::Vector3D
Acts::LineSurface::normal(const Acts::Vector2D& /*lpos*/) const
{
  // the normal is conceptionally closest to the line direction
  return lineDirection();
}

const Acts::SurfaceBounds&
Acts::LineSurface::bounds() const
{
  if (m_bounds) return (*m_bounds.get());
  return s_noBounds;
}

Acts::variant_data
Acts::LineSurface::toVariantData() const
{
  using namespace std::string_literals;

  variant_map payload;

  variant_data bounds = m_bounds->toVariantData();
  payload["bounds"]   = bounds;

  if (m_transform) {
    payload["transform"] = to_variant(*m_transform);
  }

  variant_map data;
  data["type"]    = "LineSurface"s;
  data["payload"] = payload;
  return data;
}
