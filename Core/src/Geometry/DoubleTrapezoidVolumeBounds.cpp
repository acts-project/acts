// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// DoubleTrapezoidVolumeBounds.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Geometry/DoubleTrapezoidVolumeBounds.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include "Acts/Surfaces/DiamondBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"

Acts::DoubleTrapezoidVolumeBounds::DoubleTrapezoidVolumeBounds()
    : VolumeBounds(), m_boundValues(bv_length, 0.) {}

Acts::DoubleTrapezoidVolumeBounds::DoubleTrapezoidVolumeBounds(
    double minhalex, double medhalex, double maxhalex, double haley1,
    double haley2, double halez)
    : VolumeBounds(), m_boundValues(bv_length, 0.) {
  m_boundValues.at(bv_minHalfX) = minhalex;
  m_boundValues.at(bv_medHalfX) = medhalex;
  m_boundValues.at(bv_maxHalfX) = maxhalex;
  m_boundValues.at(bv_halfY1) = haley1;
  m_boundValues.at(bv_halfY2) = haley2;
  m_boundValues.at(bv_halfZ) = halez;
  m_boundValues.at(bv_alpha1) =
      atan2(m_boundValues.at(bv_medHalfX) - m_boundValues.at(bv_minHalfX),
            2. * m_boundValues.at(bv_halfY1));
  m_boundValues.at(bv_alpha2) =
      atan2(m_boundValues.at(bv_medHalfX) - m_boundValues.at(bv_maxHalfX),
            2. * m_boundValues.at(bv_halfY2));
}

Acts::DoubleTrapezoidVolumeBounds::DoubleTrapezoidVolumeBounds(
    const Acts::DoubleTrapezoidVolumeBounds& trabo)
    : VolumeBounds(), m_boundValues(trabo.m_boundValues) {}

Acts::DoubleTrapezoidVolumeBounds::~DoubleTrapezoidVolumeBounds() = default;

Acts::DoubleTrapezoidVolumeBounds& Acts::DoubleTrapezoidVolumeBounds::operator=(
    const Acts::DoubleTrapezoidVolumeBounds& trabo) {
  if (this != &trabo) {
    m_boundValues = trabo.m_boundValues;
  }
  return *this;
}

std::vector<std::shared_ptr<const Acts::Surface>>
Acts::DoubleTrapezoidVolumeBounds::decomposeToSurfaces(
    const Transform3D* transformPtr) const {
  std::vector<std::shared_ptr<const Surface>> rSurfaces;

  // the transform
  Transform3D transform =
      (transformPtr == nullptr) ? Transform3D::Identity() : (*transformPtr);

  // face surfaces xy
  RotationMatrix3D diamondRotation(transform.rotation());
  Vector3D diamondX(diamondRotation.col(0));
  Vector3D diamondY(diamondRotation.col(1));
  Vector3D diamondZ(diamondRotation.col(2));
  Vector3D diamondCenter(transform.translation());

  const Transform3D* tTransform = nullptr;

  //   (1) - at negative local z
  tTransform =
      new Transform3D(transform * AngleAxis3D(M_PI, Vector3D(0., 1., 0.)) *
                      Translation3D(Vector3D(0., 0., halflengthZ())));
  rSurfaces.push_back(Surface::makeShared<PlaneSurface>(
      std::shared_ptr<const Transform3D>(tTransform),
      std::shared_ptr<const PlanarBounds>(faceXYDiamondBounds())));
  //   (2) - at positive local z
  tTransform = new Transform3D(transform *
                               Translation3D(Vector3D(0., 0., halflengthZ())));
  rSurfaces.push_back(Surface::makeShared<PlaneSurface>(
      std::shared_ptr<const Transform3D>(tTransform),
      std::shared_ptr<const PlanarBounds>(faceXYDiamondBounds())));
  // face surfaces yz
  // transmute cyclical
  //   (3) - at point A, attached to alpha opening angle
  Vector3D A(diamondCenter - minHalflengthX() * diamondX -
             2 * halflengthY1() * diamondY);
  AngleAxis3D alpha1ZRotation(alpha1(), Vector3D(0., 0., 1.));
  RotationMatrix3D alpha1Rotation(
      diamondRotation * alpha1ZRotation *
      AngleAxis3D(-0.5 * M_PI, Vector3D(0., 1., 0.)) *
      AngleAxis3D(0.5 * M_PI, Vector3D(0., 0., 1.)));
  RectangleBounds* faceAlpha1Bounds = faceAlpha1RectangleBounds();
  Vector3D faceAlpha1Position(A + alpha1Rotation.col(0) *
                                      faceAlpha1Bounds->halflengthX());
  tTransform =
      new Transform3D(Translation3D(faceAlpha1Position) * alpha1Rotation);
  rSurfaces.push_back(Surface::makeShared<PlaneSurface>(
      std::shared_ptr<const Transform3D>(tTransform),
      std::shared_ptr<const PlanarBounds>(faceAlpha1Bounds)));
  //   (4) - at point B, attached to beta opening angle
  Vector3D B(diamondCenter + minHalflengthX() * diamondX -
             2 * halflengthY1() * diamondY);
  AngleAxis3D beta1ZRotation(-alpha1(), Vector3D(0., 0., 1.));
  RotationMatrix3D beta1Rotation(diamondRotation * beta1ZRotation *
                                 AngleAxis3D(0.5 * M_PI, Vector3D(0., 1., 0.)) *
                                 AngleAxis3D(0.5 * M_PI, Vector3D(0., 0., 1.)));
  RectangleBounds* faceBeta1Bounds = faceBeta1RectangleBounds();
  Vector3D faceBeta1Position(B + beta1Rotation.col(0) *
                                     faceBeta1Bounds->halflengthX());
  tTransform =
      new Transform3D(Translation3D(faceBeta1Position) * beta1Rotation);
  rSurfaces.push_back(Surface::makeShared<PlaneSurface>(
      std::shared_ptr<const Transform3D>(tTransform),
      std::shared_ptr<const PlanarBounds>(faceBeta1Bounds)));
  // face surfaces yz
  // transmute cyclical
  //   (5) - at point A', attached to alpha opening angle
  Vector3D AA(diamondCenter - maxHalflengthX() * diamondX +
              2 * halflengthY2() * diamondY);
  AngleAxis3D alpha2ZRotation(-alpha2(), Vector3D(0., 0., 1.));
  RotationMatrix3D alpha2Rotation(
      diamondRotation * alpha2ZRotation *
      AngleAxis3D(-0.5 * M_PI, Vector3D(0., 1., 0.)) *
      AngleAxis3D(-0.5 * M_PI, Vector3D(0., 0., 1.)));
  RectangleBounds* faceAlpha2Bounds = faceAlpha2RectangleBounds();
  Vector3D faceAlpha2Position(AA + alpha2Rotation.col(0) *
                                       faceAlpha2Bounds->halflengthX());
  tTransform =
      new Transform3D(Translation3D(faceAlpha2Position) * alpha2Rotation);
  rSurfaces.push_back(Surface::makeShared<PlaneSurface>(
      std::shared_ptr<const Transform3D>(tTransform),
      std::shared_ptr<const PlanarBounds>(faceAlpha2Bounds)));
  //   (6) - at point B', attached to beta opening angle
  Vector3D BB(diamondCenter + maxHalflengthX() * diamondX +
              2 * halflengthY2() * diamondY);
  AngleAxis3D beta2ZRotation(alpha2(), Vector3D(0., 0., 1.));
  RotationMatrix3D beta2Rotation(
      diamondRotation * beta2ZRotation *
      AngleAxis3D(0.5 * M_PI, Vector3D(0., 1., 0.)) *
      AngleAxis3D(-0.5 * M_PI, Vector3D(0., 0., 1.)));
  RectangleBounds* faceBeta2Bounds = faceBeta2RectangleBounds();
  Vector3D faceBeta2Position(BB + beta2Rotation.col(0) *
                                      faceBeta2Bounds->halflengthX());
  tTransform =
      new Transform3D(Translation3D(faceBeta2Position) * beta2Rotation);
  rSurfaces.push_back(Surface::makeShared<PlaneSurface>(
      std::shared_ptr<const Transform3D>(tTransform),
      std::shared_ptr<const PlanarBounds>(faceBeta2Bounds)));
  // face surfaces zx
  //   (7) - at negative local y
  tTransform =
      new Transform3D(transform * AngleAxis3D(M_PI, Vector3D(1., 0., 0.)) *
                      Translation3D(Vector3D(0., 2 * halflengthY1(), 0.)) *
                      AngleAxis3D(-0.5 * M_PI, Vector3D(0., 1., 0.)) *
                      AngleAxis3D(-0.5 * M_PI, Vector3D(1., 0., 0.)));
  rSurfaces.push_back(Surface::makeShared<PlaneSurface>(
      std::shared_ptr<const Transform3D>(tTransform),
      std::shared_ptr<const PlanarBounds>(faceZXRectangleBoundsBottom())));
  //   (8) - at positive local y
  tTransform = new Transform3D(
      transform * Translation3D(Vector3D(0., 2 * halflengthY2(), 0.)) *
      AngleAxis3D(-0.5 * M_PI, Vector3D(0., 1., 0.)) *
      AngleAxis3D(-0.5 * M_PI, Vector3D(1., 0., 0.)));
  rSurfaces.push_back(Surface::makeShared<PlaneSurface>(
      std::shared_ptr<const Transform3D>(tTransform),
      std::shared_ptr<const PlanarBounds>(faceZXRectangleBoundsTop())));

  return rSurfaces;
}

// faces in xy
Acts::DiamondBounds* Acts::DoubleTrapezoidVolumeBounds::faceXYDiamondBounds()
    const {
  return new DiamondBounds(
      m_boundValues.at(bv_minHalfX), m_boundValues.at(bv_medHalfX),
      m_boundValues.at(bv_maxHalfX), 2 * m_boundValues.at(bv_halfY1),
      2 * m_boundValues.at(bv_halfY2));
}

Acts::RectangleBounds*
Acts::DoubleTrapezoidVolumeBounds::faceAlpha1RectangleBounds() const {
  return new RectangleBounds(
      m_boundValues.at(bv_halfY1) / cos(m_boundValues.at(bv_alpha1)),
      m_boundValues.at(bv_halfZ));
}

Acts::RectangleBounds*
Acts::DoubleTrapezoidVolumeBounds::faceAlpha2RectangleBounds() const {
  return new RectangleBounds(
      m_boundValues.at(bv_halfY2) / cos(m_boundValues.at(bv_alpha2)),
      m_boundValues.at(bv_halfZ));
}

Acts::RectangleBounds*
Acts::DoubleTrapezoidVolumeBounds::faceBeta1RectangleBounds() const {
  return new RectangleBounds(
      m_boundValues.at(bv_halfY1) / cos(m_boundValues.at(bv_alpha1)),
      m_boundValues.at(bv_halfZ));
}

Acts::RectangleBounds*
Acts::DoubleTrapezoidVolumeBounds::faceBeta2RectangleBounds() const {
  return new RectangleBounds(
      m_boundValues.at(bv_halfY2) / cos(m_boundValues.at(bv_alpha2)),
      m_boundValues.at(bv_halfZ));
}

Acts::RectangleBounds*
Acts::DoubleTrapezoidVolumeBounds::faceZXRectangleBoundsBottom() const {
  return new RectangleBounds(m_boundValues.at(bv_halfZ),
                             m_boundValues.at(bv_minHalfX));
}

Acts::RectangleBounds*
Acts::DoubleTrapezoidVolumeBounds::faceZXRectangleBoundsTop() const {
  return new RectangleBounds(m_boundValues.at(bv_halfZ),
                             m_boundValues.at(bv_maxHalfX));
}

bool Acts::DoubleTrapezoidVolumeBounds::inside(const Vector3D& pos,
                                               double tol) const {
  if (std::abs(pos.z()) > m_boundValues.at(bv_halfZ) + tol) {
    return false;
  }
  if (pos.y() < -2 * m_boundValues.at(bv_halfY1) - tol) {
    return false;
  }
  if (pos.y() > 2 * m_boundValues.at(bv_halfY2) - tol) {
    return false;
  }
  DiamondBounds* faceXYBounds = faceXYDiamondBounds();
  Vector2D locp(pos.x(), pos.y());
  bool inside(faceXYBounds->inside(locp, BoundaryCheck(true, true, tol, tol)));
  delete faceXYBounds;
  return inside;
}

// ostream operator overload
std::ostream& Acts::DoubleTrapezoidVolumeBounds::toStream(
    std::ostream& sl) const {
  return dumpT<std::ostream>(sl);
}

Acts::Volume::BoundingBox Acts::DoubleTrapezoidVolumeBounds::boundingBox(
    const Transform3D* trf, const Vector3D& envelope,
    const Volume* entity) const {
  float minx = minHalflengthX();
  float medx = medHalflengthX();
  float maxx = maxHalflengthX();
  float haley1 = 2 * halflengthY1();
  float haley2 = 2 * halflengthY2();
  float halez = halflengthZ();

  std::array<Vector3D, 12> vertices = {{
      {-minx, -haley1, -halez},
      {+minx, -haley1, -halez},
      {+medx, 0, -halez},
      {-medx, 0, -halez},
      {-maxx, +haley2, -halez},
      {+maxx, +haley2, -halez},
      {-minx, -haley1, +halez},
      {+minx, -haley1, +halez},
      {+medx, 0, +halez},
      {-medx, 0, +halez},
      {-maxx, +haley2, +halez},
      {+maxx, +haley2, +halez},
  }};

  Transform3D transform = Transform3D::Identity();
  if (trf != nullptr) {
    transform = *trf;
  }

  Vector3D vmin = transform * vertices[0];
  Vector3D vmax = transform * vertices[0];

  for (size_t i = 1; i < 12; i++) {
    const Vector3D vtx = transform * vertices[i];
    vmin = vmin.cwiseMin(vtx);
    vmax = vmax.cwiseMax(vtx);
  }

  return {entity, vmin - envelope, vmax + envelope};
}
