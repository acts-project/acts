// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// TrapezoidVolumeBounds.cpp, Acts project
///////////////////////////////////////////////////////////////////

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
  m_values[eAlpha] =
      atan((m_values[eHalfLengthXposY] - m_values[eHalfLengthXnegY]) / 2 /
           m_values[eHalfLengthY]) +
      0.5 * M_PI;
  m_values[eBeta] = m_values[eAlpha];
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

std::vector<std::shared_ptr<const Acts::Surface>>
Acts::TrapezoidVolumeBounds::decomposeToSurfaces(
    const Transform3D* transformPtr) const {
  std::vector<std::shared_ptr<const Surface>> rSurfaces;

  Transform3D transform =
      (transformPtr == nullptr) ? Transform3D::Identity() : (*transformPtr);
  const Transform3D* tTransform = nullptr;
  // face surfaces xy
  RotationMatrix3D trapezoidRotation(transform.rotation());
  Vector3D trapezoidX(trapezoidRotation.col(0));
  Vector3D trapezoidY(trapezoidRotation.col(1));
  Vector3D trapezoidZ(trapezoidRotation.col(2));
  Vector3D trapezoidCenter(transform.translation());

  //   (1) - at negative local z
  tTransform =
      new Transform3D(transform * AngleAxis3D(M_PI, Vector3D(0., 1., 0.)) *
                      Translation3D(Vector3D(0., 0., get(eHalfLengthZ))));

  rSurfaces.push_back(Surface::makeShared<PlaneSurface>(
      std::shared_ptr<const Transform3D>(tTransform), m_faceXYTrapezoidBounds));
  //   (2) - at positive local z
  tTransform = new Transform3D(
      transform * Translation3D(Vector3D(0., 0., get(eHalfLengthZ))));
  rSurfaces.push_back(Surface::makeShared<PlaneSurface>(
      std::shared_ptr<const Transform3D>(tTransform), m_faceXYTrapezoidBounds));

  // face surfaces yz
  // transmute cyclical
  //   (3) - at point A, attached to alpha opening angle
  Vector3D A(get(eHalfLengthXnegY), get(eHalfLengthY), trapezoidCenter.z());
  RotationMatrix3D alphaZRotation =
      (s_idRotation *
       AngleAxis3D(get(eAlpha) - 0.5 * M_PI, Vector3D(0., 0., 1.)))
          .toRotationMatrix();
  RotationMatrix3D faceAlphaRotation;
  faceAlphaRotation.col(0) = alphaZRotation.col(1);
  faceAlphaRotation.col(1) = -alphaZRotation.col(2);
  faceAlphaRotation.col(2) = -alphaZRotation.col(0);

  // Vector3D
  // faceAlphaPosition(A+faceAlphaRotation.colX()*m_faceAlphaRectangleBounds->halflengthX());
  Vector3D faceAlphaPosition0(
      -0.5 * (get(eHalfLengthXnegY) + get(eHalfLengthXposY)), 0., 0.);
  Vector3D faceAlphaPosition = transform * faceAlphaPosition0;
  tTransform = new Transform3D(Translation3D(faceAlphaPosition) *
                               (trapezoidRotation * faceAlphaRotation));
  rSurfaces.push_back(Surface::makeShared<PlaneSurface>(
      std::shared_ptr<const Transform3D>(tTransform),
      m_faceAlphaRectangleBounds));

  //   (4) - at point B, attached to beta opening angle
  Vector3D B(get(eHalfLengthXnegY), -get(eHalfLengthY), trapezoidCenter.z());
  RotationMatrix3D betaZRotation =
      (s_idRotation *
       AngleAxis3D(-(get(eBeta) - 0.5 * M_PI), Vector3D(0., 0., 1.)))
          .toRotationMatrix();
  RotationMatrix3D faceBetaRotation;
  faceBetaRotation.col(0) = betaZRotation.col(1);
  faceBetaRotation.col(1) = betaZRotation.col(2);
  faceBetaRotation.col(2) = betaZRotation.col(0);
  // Vector3D
  // faceBetaPosition(B+faceBetaRotation.colX()*m_faceBetaRectangleBounds->halflengthX());
  Vector3D faceBetaPosition0(
      0.5 * (get(eHalfLengthXnegY) + get(eHalfLengthXposY)), 0., 0.);
  Vector3D faceBetaPosition = transform * faceBetaPosition0;
  tTransform = new Transform3D(Translation3D(faceBetaPosition) *
                               (trapezoidRotation * faceBetaRotation));
  rSurfaces.push_back(Surface::makeShared<PlaneSurface>(
      std::shared_ptr<const Transform3D>(tTransform),
      m_faceBetaRectangleBounds));

  // face surfaces zx
  //   (5) - at negative local y
  tTransform =
      new Transform3D(transform * AngleAxis3D(M_PI, Vector3D(1., 0., 0.)) *
                      Translation3D(Vector3D(0., get(eHalfLengthY), 0.)) *
                      AngleAxis3D(-0.5 * M_PI, Vector3D(0., 1., 0.)) *
                      AngleAxis3D(-0.5 * M_PI, Vector3D(1., 0., 0.)));
  rSurfaces.push_back(Surface::makeShared<PlaneSurface>(
      std::shared_ptr<const Transform3D>(tTransform),
      std::shared_ptr<const PlanarBounds>(m_faceZXRectangleBoundsBottom)));
  //   (6) - at positive local y
  tTransform = new Transform3D(
      transform * Translation3D(Vector3D(0., get(eHalfLengthY), 0.)) *
      AngleAxis3D(-0.5 * M_PI, Vector3D(0., 1., 0.)) *
      AngleAxis3D(-0.5 * M_PI, Vector3D(1., 0., 0.)));
  rSurfaces.push_back(Surface::makeShared<PlaneSurface>(
      std::shared_ptr<const Transform3D>(tTransform),
      std::shared_ptr<const PlanarBounds>(m_faceZXRectangleBoundsTop)));

  return rSurfaces;
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
