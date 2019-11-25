// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// ConeSurface.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Surfaces/ConeSurface.hpp"

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <utility>

#include "Acts/Utilities/ThrowAssert.hpp"
#include "Acts/Utilities/detail/RealQuadraticEquation.hpp"

using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;

Acts::ConeSurface::ConeSurface(const ConeSurface& other)
    : GeometryObject(), Surface(other), m_bounds(other.m_bounds) {}

Acts::ConeSurface::ConeSurface(const GeometryContext& gctx,
                               const ConeSurface& other,
                               const Transform3D& transf)
    : GeometryObject(),
      Surface(gctx, other, transf),
      m_bounds(other.m_bounds) {}

Acts::ConeSurface::ConeSurface(std::shared_ptr<const Transform3D> htrans,
                               double alpha, bool symmetric)
    : GeometryObject(),
      Surface(std::move(htrans)),
      m_bounds(std::make_shared<const ConeBounds>(alpha, symmetric)) {}

Acts::ConeSurface::ConeSurface(std::shared_ptr<const Transform3D> htrans,
                               double alpha, double zmin, double zmax,
                               double halfPhi)
    : GeometryObject(),
      Surface(std::move(htrans)),
      m_bounds(std::make_shared<const ConeBounds>(alpha, zmin, zmax, halfPhi)) {
}

Acts::ConeSurface::ConeSurface(std::shared_ptr<const Transform3D> htrans,
                               const std::shared_ptr<const ConeBounds>& cbounds)
    : GeometryObject(), Surface(std::move(htrans)), m_bounds(cbounds) {
  throw_assert(cbounds, "ConeBounds must not be nullptr");
}

const Acts::Vector3D Acts::ConeSurface::binningPosition(
    const GeometryContext& gctx, Acts::BinningValue bValue) const {
  const Vector3D& sfCenter = center(gctx);

  // special binning type for R-type methods
  if (bValue == Acts::binR || bValue == Acts::binRPhi) {
    return Vector3D(sfCenter.x() + bounds().r(sfCenter.z()), sfCenter.y(),
                    sfCenter.z());
  }
  // give the center as default for all of these binning types
  // binX, binY, binZ, binR, binPhi, binRPhi, binH, binEta
  return sfCenter;
}

Acts::Surface::SurfaceType Acts::ConeSurface::type() const {
  return Surface::Cone;
}

Acts::ConeSurface& Acts::ConeSurface::operator=(const ConeSurface& other) {
  if (this != &other) {
    Surface::operator=(other);
    m_bounds = other.m_bounds;
  }
  return *this;
}

const Acts::Vector3D Acts::ConeSurface::rotSymmetryAxis(
    const GeometryContext& gctx) const {
  return std::move(transform(gctx).matrix().block<3, 1>(0, 2));
}

const Acts::RotationMatrix3D Acts::ConeSurface::referenceFrame(
    const GeometryContext& gctx, const Vector3D& position,
    const Vector3D& /*unused*/) const {
  RotationMatrix3D mFrame;
  // construct the measurement frame
  // measured Y is the local z axis
  Vector3D measY = rotSymmetryAxis(gctx);
  // measured z is the position transverse normalized
  Vector3D measDepth = Vector3D(position.x(), position.y(), 0.).normalized();
  // measured X is what comoes out of it
  Acts::Vector3D measX(measY.cross(measDepth).normalized());
  // the columnes
  mFrame.col(0) = measX;
  mFrame.col(1) = measY;
  mFrame.col(2) = measDepth;
  // return the rotation matrix
  //!< @todo fold in alpha
  // return it
  return mFrame;
}

void Acts::ConeSurface::localToGlobal(const GeometryContext& gctx,
                                      const Vector2D& lposition,
                                      const Vector3D& /*unused*/,
                                      Vector3D& position) const {
  // create the position in the local 3d frame
  double r = lposition[Acts::eLOC_Z] * bounds().tanAlpha();
  double phi = lposition[Acts::eLOC_RPHI] / r;
  Vector3D loc3Dframe(r * cos(phi), r * sin(phi), lposition[Acts::eLOC_Z]);
  // transport it to the globalframe
  if (m_transform) {
    position = transform(gctx) * loc3Dframe;
  }
}

bool Acts::ConeSurface::globalToLocal(const GeometryContext& gctx,
                                      const Vector3D& position,
                                      const Vector3D& /*unused*/,
                                      Vector2D& lposition) const {
  Vector3D loc3Dframe =
      m_transform ? (transform(gctx).inverse() * position) : position;
  double r = loc3Dframe.z() * bounds().tanAlpha();
  lposition =
      Vector2D(r * atan2(loc3Dframe.y(), loc3Dframe.x()), loc3Dframe.z());
  // now decide on the quility of the transformation
  double inttol = r * 0.0001;
  inttol = (inttol < 0.01) ? 0.01 : 0.01;  // ?
  return ((std::abs(perp(loc3Dframe) - r) > inttol) ? false : true);
}

double Acts::ConeSurface::pathCorrection(const GeometryContext& gctx,
                                         const Vector3D& position,
                                         const Vector3D& direction) const {
  // (cos phi cos alpha, sin phi cos alpha, sgn z sin alpha)
  Vector3D posLocal =
      m_transform ? transform(gctx).inverse() * position : position;
  double phi = VectorHelpers::phi(posLocal);
  double sgn = posLocal.z() > 0. ? -1. : +1.;
  Vector3D normalC(cos(phi) * bounds().cosAlpha(),
                   sin(phi) * bounds().cosAlpha(), sgn * bounds().sinAlpha());
  if (m_transform) {
    normalC = transform(gctx) * normalC;
  }
  // Back to the global frame
  double cAlpha = normalC.dot(direction);
  return std::abs(1. / cAlpha);
}

std::string Acts::ConeSurface::name() const {
  return "Acts::ConeSurface";
}

std::shared_ptr<Acts::ConeSurface> Acts::ConeSurface::clone(
    const GeometryContext& gctx, const Transform3D& shift) const {
  return std::shared_ptr<ConeSurface>(this->clone_impl(gctx, shift));
}

Acts::ConeSurface* Acts::ConeSurface::clone_impl(
    const GeometryContext& gctx, const Transform3D& shift) const {
  return new ConeSurface(gctx, *this, shift);
}

const Acts::Vector3D Acts::ConeSurface::normal(
    const GeometryContext& gctx, const Acts::Vector2D& lposition) const {
  // (cos phi cos alpha, sin phi cos alpha, sgn z sin alpha)
  double phi =
             lposition[Acts::eLOC_RPHI] / (bounds().r(lposition[Acts::eLOC_Z])),
         sgn = lposition[Acts::eLOC_Z] > 0 ? -1. : +1.;
  Vector3D localNormal(cos(phi) * bounds().cosAlpha(),
                       sin(phi) * bounds().cosAlpha(),
                       sgn * bounds().sinAlpha());
  return m_transform ? Vector3D(transform(gctx).linear() * localNormal)
                     : localNormal;
}

const Acts::Vector3D Acts::ConeSurface::normal(
    const GeometryContext& gctx, const Acts::Vector3D& position) const {
  // get it into the cylinder frame if needed
  // @todo respect opening angle
  Vector3D pos3D = position;
  if (m_transform || (m_associatedDetElement != nullptr)) {
    pos3D = transform(gctx).inverse() * position;
    pos3D.z() = 0;
  }
  return pos3D.normalized();
}

const Acts::ConeBounds& Acts::ConeSurface::bounds() const {
  // is safe because no constructor w/o bounds exists
  return (*m_bounds.get());
}
