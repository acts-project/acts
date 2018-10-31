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
#include "Acts/Surfaces/PolyhedronRepresentation.hpp"

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <utility>

#include "Acts/Utilities/ThrowAssert.hpp"
#include "Acts/Utilities/VariantData.hpp"

using Acts::VectorHelpers::phi;
using Acts::VectorHelpers::perp;

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
  , Surface(std::move(htrans))
  , m_bounds(std::make_shared<const CylinderBounds>(radius, hlength))
{
}

Acts::CylinderSurface::CylinderSurface(
    std::shared_ptr<const Transform3D> htrans,
    double                             radius,
    double                             hphi,
    double                             hlength)
  : GeometryObject()
  , Surface(std::move(htrans))
  , m_bounds(std::make_shared<const CylinderBounds>(radius, hphi, hlength))
{
}

Acts::CylinderSurface::CylinderSurface(
    std::shared_ptr<const CylinderBounds> cbounds,
    const Acts::DetectorElementBase&      detelement)
  : Surface(detelement), m_bounds(std::move(cbounds))
{
  /// surfaces representing a detector element must have bounds
  assert(cbounds);
}

Acts::CylinderSurface::CylinderSurface(
    std::shared_ptr<const Transform3D>           htrans,
    const std::shared_ptr<const CylinderBounds>& cbounds)
  : Surface(std::move(htrans)), m_bounds(cbounds)
{
  throw_assert(cbounds, "CylinderBounds must not be nullptr");
}

Acts::CylinderSurface::CylinderSurface(const variant_data& vardata)
{
  throw_assert(vardata.which() == 4, "Variant data must be map");
  variant_map data = boost::get<variant_map>(vardata);
  throw_assert(data.count("type"), "Variant data must have type.");
  std::string type = data.get<std::string>("type");
  throw_assert(type == "CylinderSurface",
               "Variant data type must be CylinderSurface");

  variant_map payload    = data.get<variant_map>("payload");
  variant_map var_bounds = payload.get<variant_map>("bounds");

  m_bounds = std::make_shared<const CylinderBounds>(var_bounds);

  if (payload.count("transform") != 0u) {
    // we have a transform
    auto trf = std::make_shared<const Transform3D>(
        from_variant<Transform3D>(payload.get<variant_map>("transform")));
    m_transform = trf;
  }
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
                                      const Vector3D& /*gmom*/) const
{
  RotationMatrix3D mFrame;
  // construct the measurement frame
  // measured Y is the z axis
  Vector3D measY = rotSymmetryAxis();
  // measured z is the position normalized transverse (in local)
  Vector3D measDepth = normal(gpos);
  // measured X is what comoes out of it
  Vector3D measX(measY.cross(measDepth).normalized());
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
                                     const Vector3D& /*gmom*/,
                                     Vector3D& gpos) const
{
  // create the position in the local 3d frame
  double r   = bounds().r();
  double phi = lpos[Acts::eLOC_RPHI] / r;
  gpos       = Vector3D(r * cos(phi), r * sin(phi), lpos[Acts::eLOC_Z]);
  // transform it to the globalframe: CylinderSurfaces are allowed to have 0
  // if pointer transform exists -> port into frame
  if (Surface::m_transform) {
    gpos = transform() * gpos;
  }
}

bool
Acts::CylinderSurface::globalToLocal(const Vector3D& gpos,
                                     const Vector3D& /*gmom*/,
                                     Vector2D& lpos) const
{
  // get the transform & transform global position into cylinder frame
  // @todo clean up intolerance parameters
  // transform it to the globalframe: CylinderSurfaces are allowed to have 0
  // pointer transform
  double radius = 0.;
  double inttol = bounds().r() * 0.0001;
  if (inttol < 0.01) {
    inttol = 0.01;
  }
  // do the transformation or not
  if (Surface::m_transform) {
    const Transform3D& surfaceTrans = transform();
    Transform3D        inverseTrans(surfaceTrans.inverse());
    Vector3D           loc3Dframe(inverseTrans * gpos);
    lpos   = Vector2D(bounds().r() * phi(loc3Dframe), loc3Dframe.z());
    radius = perp(loc3Dframe);
  } else {
    lpos   = Vector2D(bounds().r() * phi(gpos), gpos.z());
    radius = perp(gpos);
  }
  // return true or false
  return ((std::abs(radius - bounds().r()) > inttol) ? false : true);
}

std::string
Acts::CylinderSurface::name() const
{
  return "Acts::CylinderSurface";
}

std::shared_ptr<Acts::CylinderSurface>
Acts::CylinderSurface::clone(const Transform3D* shift) const
{
  return std::shared_ptr<CylinderSurface>(this->clone_impl(shift));
}

Acts::CylinderSurface*
Acts::CylinderSurface::clone_impl(const Transform3D* shift) const
{
  if (shift != nullptr) {
    return new CylinderSurface(*this, *shift);
  }
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
  Vector3D pos3D      = gpos;
  bool needsTransform = (m_transform || (m_associatedDetElement != nullptr));
  if (needsTransform) {
    pos3D = transform().inverse() * gpos;
  }
  // set the z coordinate to 0
  pos3D.z() = 0.;
  // normalize and rotate back into global if needed
  return needsTransform ? transform().linear() * pos3D.normalized()
                        : pos3D.normalized();
}

double
Acts::CylinderSurface::pathCorrection(const Acts::Vector3D& gpos,
                                      const Acts::Vector3D& mom) const
{
  Vector3D normalT  = normal(gpos);
  double   cosAlpha = normalT.dot(mom.normalized());
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

Acts::PolyhedronRepresentation
Acts::CylinderSurface::polyhedronRepresentation(size_t l0div,
                                                size_t /*l1div*/) const
{
  std::vector<Vector3D>            vertices;
  std::vector<std::vector<size_t>> faces;

  if (l0div <= 1) {
    throw std::domain_error(
        "Polyhedron repr of cylinder with 1 div is undefined");
  }

  double phistep = 2 * M_PI / l0div;
  double hlZ     = bounds().halflengthZ();
  double r       = bounds().r();

  Vector3D left(r, 0, -hlZ);
  Vector3D right(r, 0, hlZ);

  for (size_t i = 0; i < l0div; i++) {
    Transform3D rot(AngleAxis3D(i * phistep, Vector3D::UnitZ()));
    vertices.push_back(transform() * rot * left);
    vertices.push_back(transform() * rot * right);
  }

  for (size_t v = 0; v < vertices.size() - 2; v = v + 2) {

    faces.push_back({v, v + 1, v + 3, v + 2});
  }
  if (l0div > 2) {
    faces.push_back({vertices.size() - 2, vertices.size() - 1, 1, 0});
  }

  return PolyhedronRepresentation(vertices, faces);
}
