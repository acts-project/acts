// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/TrapezoidVolumeBounds.hpp"
#include "Acts/Geometry/GeometryStatics.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Utilities/BoundingBox.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>

Acts::TrapezoidVolumeBounds::TrapezoidVolumeBounds(double minhalex,
                                                   double maxhalex,
                                                   double haley, double halez)
    : VolumeBounds() {
  m_values[eHalfLengthXnegY] = minhalex;
  m_values[eHalfLengthXposY] = maxhalex;
  m_values[eHalfLengthY] = haley;
  m_values[eHalfLengthZ] = halez;
  m_values[eAlpha] = atan2(2 * haley, (maxhalex - minhalex));
  m_values[eBeta] = M_PI - get(eAlpha);
  checkConsistency();
  buildSurfaceBounds();
}

Acts::TrapezoidVolumeBounds::TrapezoidVolumeBounds(double minhalex,
                                                   double haley, double halez,
                                                   double alpha, double beta)
    : VolumeBounds() {
  m_values[eHalfLengthXnegY] = minhalex;
  m_values[eHalfLengthY] = haley;
  m_values[eHalfLengthZ] = halez;
  m_values[eAlpha] = alpha;
  m_values[eBeta] = beta;
  // now calculate the remaining max half X
  double gamma = (alpha > beta) ? (alpha - 0.5 * M_PI) : (beta - 0.5 * M_PI);
  m_values[eHalfLengthXposY] = minhalex + (2. * haley) * tan(gamma);

  checkConsistency();
  buildSurfaceBounds();
}

Acts::OrientedSurfaces Acts::TrapezoidVolumeBounds::orientedSurfaces(
    const Transform3D* transformPtr) const {
  Transform3D transform =
      (transformPtr == nullptr) ? Transform3D::Identity() : (*transformPtr);

  OrientedSurfaces oSurfaces;
  oSurfaces.reserve(6);

  // Face surfaces xy
  RotationMatrix3D trapezoidRotation(transform.rotation());
  Vector3D trapezoidX(trapezoidRotation.col(0));
  Vector3D trapezoidY(trapezoidRotation.col(1));
  Vector3D trapezoidZ(trapezoidRotation.col(2));
  Vector3D trapezoidCenter(transform.translation());

  //   (1) - At negative local z
  auto nzTransform = std::make_shared<Transform3D>(
      transform * Translation3D(0., 0., -get(eHalfLengthZ)));
  auto sf =
      Surface::makeShared<PlaneSurface>(nzTransform, m_faceXYTrapezoidBounds);
  oSurfaces.push_back(OrientedSurface(std::move(sf), forward));
  //   (2) - At positive local z
  auto pzTransform = std::make_shared<Transform3D>(
      transform * Translation3D(0., 0., get(eHalfLengthZ)));
  sf = Surface::makeShared<PlaneSurface>(pzTransform, m_faceXYTrapezoidBounds);
  oSurfaces.push_back(OrientedSurface(std::move(sf), backward));

  double poshOffset = get(eHalfLengthY) / std::tan(get(eAlpha));
  double neghOffset = get(eHalfLengthY) / std::tan(get(eBeta));
  double topShift = poshOffset + neghOffset;

  // Face surfaces yz
  // (3) - At point B, attached to beta opening angle
  Vector3D fbPosition(-get(eHalfLengthXnegY) + neghOffset, 0., 0.);
  auto fbTransform = std::make_shared<Transform3D>(
      transform * Translation3D(fbPosition) *
      AngleAxis3D(-0.5 * M_PI + get(eBeta), Vector3D(0., 0., 1.)) * s_planeYZ);
  sf =
      Surface::makeShared<PlaneSurface>(fbTransform, m_faceBetaRectangleBounds);
  oSurfaces.push_back(OrientedSurface(std::move(sf), forward));

  // (4) - At point A, attached to alpha opening angle
  Vector3D faPosition(get(eHalfLengthXnegY) + poshOffset, 0., 0.);
  auto faTransform = std::make_shared<Transform3D>(
      transform * Translation3D(faPosition) *
      AngleAxis3D(-0.5 * M_PI + get(eAlpha), Vector3D(0., 0., 1.)) * s_planeYZ);
  sf = Surface::makeShared<PlaneSurface>(faTransform,
                                         m_faceAlphaRectangleBounds);
  oSurfaces.push_back(OrientedSurface(std::move(sf), backward));

  // Face surfaces zx
  //   (5) - At negative local y
  auto nxTransform = std::make_shared<Transform3D>(
      transform * Translation3D(0., -get(eHalfLengthY), 0.) * s_planeZX);
  sf = Surface::makeShared<PlaneSurface>(nxTransform,
                                         m_faceZXRectangleBoundsBottom);
  oSurfaces.push_back(OrientedSurface(std::move(sf), forward));
  //   (6) - At positive local y
  auto pxTransform = std::make_shared<Transform3D>(
      transform * Translation3D(topShift, get(eHalfLengthY), 0.) * s_planeZX);
  sf = Surface::makeShared<PlaneSurface>(pxTransform,
                                         m_faceZXRectangleBoundsTop);
  oSurfaces.push_back(OrientedSurface(std::move(sf), backward));

  return oSurfaces;
}

void Acts::TrapezoidVolumeBounds::buildSurfaceBounds() {
  m_faceXYTrapezoidBounds = std::make_shared<const TrapezoidBounds>(
      get(eHalfLengthXnegY), get(eHalfLengthXposY), get(eHalfLengthY));

  m_faceAlphaRectangleBounds = std::make_shared<const RectangleBounds>(
      get(eHalfLengthY) / cos(get(eAlpha) - 0.5 * M_PI), get(eHalfLengthZ));

  m_faceBetaRectangleBounds = std::make_shared<const RectangleBounds>(
      get(eHalfLengthY) / cos(get(eBeta) - 0.5 * M_PI), get(eHalfLengthZ));

  m_faceZXRectangleBoundsBottom = std::make_shared<const RectangleBounds>(
      get(eHalfLengthZ), get(eHalfLengthXnegY));

  m_faceZXRectangleBoundsTop = std::make_shared<const RectangleBounds>(
      get(eHalfLengthZ), get(eHalfLengthXposY));
}

bool Acts::TrapezoidVolumeBounds::inside(const Vector3D& pos,
                                         double tol) const {
  if (std::abs(pos.z()) > get(eHalfLengthZ) + tol) {
    return false;
  }
  if (std::abs(pos.y()) > get(eHalfLengthY) + tol) {
    return false;
  }
  Vector2D locp(pos.x(), pos.y());
  bool inside(m_faceXYTrapezoidBounds->inside(
      locp, BoundaryCheck(true, true, tol, tol)));
  return inside;
}

std::ostream& Acts::TrapezoidVolumeBounds::toStream(std::ostream& sl) const {
  return dumpT<std::ostream>(sl);
}

Acts::Volume::BoundingBox Acts::TrapezoidVolumeBounds::boundingBox(
    const Acts::Transform3D* trf, const Vector3D& envelope,
    const Volume* entity) const {
  double minx = get(eHalfLengthXnegY);
  double maxx = get(eHalfLengthXposY);
  double haley = get(eHalfLengthY);
  double halez = get(eHalfLengthZ);

  std::array<Vector3D, 8> vertices = {{{-minx, -haley, -halez},
                                       {+minx, -haley, -halez},
                                       {-maxx, +haley, -halez},
                                       {+maxx, +haley, -halez},
                                       {-minx, -haley, +halez},
                                       {+minx, -haley, +halez},
                                       {-maxx, +haley, +halez},
                                       {+maxx, +haley, +halez}}};

  Transform3D transform = Transform3D::Identity();
  if (trf != nullptr) {
    transform = *trf;
  }

  Vector3D vmin = transform * vertices[0];
  Vector3D vmax = transform * vertices[0];

  for (size_t i = 1; i < 8; i++) {
    const Vector3D vtx = transform * vertices[i];
    vmin = vmin.cwiseMin(vtx);
    vmax = vmax.cwiseMax(vtx);
  }

  return {entity, vmin - envelope, vmax + envelope};
}
