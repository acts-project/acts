// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// DoubleTrapezoidVolumeBounds.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Volumes/DoubleTrapezoidVolumeBounds.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include "Acts/Surfaces/DiamondBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"

Acts::DoubleTrapezoidVolumeBounds::DoubleTrapezoidVolumeBounds()
  : VolumeBounds(), m_valueStore(bv_length, 0.)
{
}

Acts::DoubleTrapezoidVolumeBounds::DoubleTrapezoidVolumeBounds(double minhalex,
                                                               double medhalex,
                                                               double maxhalex,
                                                               double haley1,
                                                               double haley2,
                                                               double halez)
  : VolumeBounds(), m_valueStore(bv_length, 0.)
{
  m_valueStore.at(bv_minHalfX) = minhalex;
  m_valueStore.at(bv_medHalfX) = medhalex;
  m_valueStore.at(bv_maxHalfX) = maxhalex;
  m_valueStore.at(bv_halfY1)   = haley1;
  m_valueStore.at(bv_halfY2)   = haley2;
  m_valueStore.at(bv_halfZ)    = halez;
  m_valueStore.at(bv_alpha1)
      = atan2(m_valueStore.at(bv_medHalfX) - m_valueStore.at(bv_minHalfX),
              2. * m_valueStore.at(bv_halfY1));
  m_valueStore.at(bv_alpha2)
      = atan2(m_valueStore.at(bv_medHalfX) - m_valueStore.at(bv_maxHalfX),
              2. * m_valueStore.at(bv_halfY2));
}

Acts::DoubleTrapezoidVolumeBounds::DoubleTrapezoidVolumeBounds(
    const Acts::DoubleTrapezoidVolumeBounds& trabo)
  : VolumeBounds(), m_valueStore(trabo.m_valueStore)
{
}

Acts::DoubleTrapezoidVolumeBounds::~DoubleTrapezoidVolumeBounds() = default;

Acts::DoubleTrapezoidVolumeBounds&
Acts::DoubleTrapezoidVolumeBounds::
operator=(const Acts::DoubleTrapezoidVolumeBounds& trabo)
{
  if (this != &trabo) {
    m_valueStore = trabo.m_valueStore;
  }
  return *this;
}

std::vector<std::shared_ptr<const Acts::Surface>>
Acts::DoubleTrapezoidVolumeBounds::decomposeToSurfaces(
    std::shared_ptr<const Transform3D> transformPtr) const
{
  std::vector<std::shared_ptr<const Surface>> rSurfaces;

  // the transform
  Transform3D transform = (transformPtr == nullptr) ? Transform3D::Identity()
                                                    : (*transformPtr.get());

  // face surfaces xy
  RotationMatrix3D diamondRotation(transform.rotation());
  Vector3D         diamondX(diamondRotation.col(0));
  Vector3D         diamondY(diamondRotation.col(1));
  Vector3D         diamondZ(diamondRotation.col(2));
  Vector3D         diamondCenter(transform.translation());

  const Transform3D* tTransform = nullptr;

  //   (1) - at negative local z
  tTransform
      = new Transform3D(transform * AngleAxis3D(M_PI, Vector3D(0., 1., 0.))
                        * Translation3D(Vector3D(0., 0., halflengthZ())));
  rSurfaces.push_back(Surface::makeShared<PlaneSurface>(
      std::shared_ptr<const Transform3D>(tTransform),
      std::shared_ptr<const PlanarBounds>(faceXYDiamondBounds())));
  //   (2) - at positive local z
  tTransform = new Transform3D(
      transform * Translation3D(Vector3D(0., 0., halflengthZ())));
  rSurfaces.push_back(Surface::makeShared<PlaneSurface>(
      std::shared_ptr<const Transform3D>(tTransform),
      std::shared_ptr<const PlanarBounds>(faceXYDiamondBounds())));
  // face surfaces yz
  // transmute cyclical
  //   (3) - at point A, attached to alpha opening angle
  Vector3D A(diamondCenter - minHalflengthX() * diamondX
             - 2 * halflengthY1() * diamondY);
  AngleAxis3D      alpha1ZRotation(alpha1(), Vector3D(0., 0., 1.));
  RotationMatrix3D alpha1Rotation(
      diamondRotation * alpha1ZRotation
      * AngleAxis3D(-0.5 * M_PI, Vector3D(0., 1., 0.))
      * AngleAxis3D(0.5 * M_PI, Vector3D(0., 0., 1.)));
  RectangleBounds* faceAlpha1Bounds = faceAlpha1RectangleBounds();
  Vector3D         faceAlpha1Position(
      A + alpha1Rotation.col(0) * faceAlpha1Bounds->halflengthX());
  tTransform
      = new Transform3D(alpha1Rotation * Translation3D(faceAlpha1Position));
  rSurfaces.push_back(Surface::makeShared<PlaneSurface>(
      std::shared_ptr<const Transform3D>(tTransform),
      std::shared_ptr<const PlanarBounds>(faceAlpha1Bounds)));
  //   (4) - at point B, attached to beta opening angle
  Vector3D B(diamondCenter + minHalflengthX() * diamondX
             - 2 * halflengthY1() * diamondY);
  AngleAxis3D      beta1ZRotation(-alpha1(), Vector3D(0., 0., 1.));
  RotationMatrix3D beta1Rotation(
      diamondRotation * beta1ZRotation
      * AngleAxis3D(0.5 * M_PI, Vector3D(0., 1., 0.))
      * AngleAxis3D(0.5 * M_PI, Vector3D(0., 0., 1.)));
  RectangleBounds* faceBeta1Bounds = faceBeta1RectangleBounds();
  Vector3D         faceBeta1Position(
      B + beta1Rotation.col(0) * faceBeta1Bounds->halflengthX());
  tTransform
      = new Transform3D(beta1Rotation * Translation3D(faceBeta1Position));
  rSurfaces.push_back(Surface::makeShared<PlaneSurface>(
      std::shared_ptr<const Transform3D>(tTransform),
      std::shared_ptr<const PlanarBounds>(faceBeta1Bounds)));
  // face surfaces yz
  // transmute cyclical
  //   (5) - at point A', attached to alpha opening angle
  Vector3D AA(diamondCenter - maxHalflengthX() * diamondX
              + 2 * halflengthY2() * diamondY);
  AngleAxis3D      alpha2ZRotation(-alpha2(), Vector3D(0., 0., 1.));
  RotationMatrix3D alpha2Rotation(
      diamondRotation * alpha2ZRotation
      * AngleAxis3D(-0.5 * M_PI, Vector3D(0., 1., 0.))
      * AngleAxis3D(-0.5 * M_PI, Vector3D(0., 0., 1.)));
  RectangleBounds* faceAlpha2Bounds = faceAlpha2RectangleBounds();
  Vector3D         faceAlpha2Position(
      AA + alpha2Rotation.col(0) * faceAlpha2Bounds->halflengthX());
  tTransform
      = new Transform3D(alpha2Rotation * Translation3D(faceAlpha2Position));
  rSurfaces.push_back(Surface::makeShared<PlaneSurface>(
      std::shared_ptr<const Transform3D>(tTransform),
      std::shared_ptr<const PlanarBounds>(faceAlpha2Bounds)));
  //   (6) - at point B', attached to beta opening angle
  Vector3D BB(diamondCenter + maxHalflengthX() * diamondX
              + 2 * halflengthY2() * diamondY);
  AngleAxis3D      beta2ZRotation(alpha2(), Vector3D(0., 0., 1.));
  RotationMatrix3D beta2Rotation(
      diamondRotation * beta2ZRotation
      * AngleAxis3D(0.5 * M_PI, Vector3D(0., 1., 0.))
      * AngleAxis3D(-0.5 * M_PI, Vector3D(0., 0., 1.)));
  RectangleBounds* faceBeta2Bounds = faceBeta2RectangleBounds();
  Vector3D         faceBeta2Position(
      BB + beta2Rotation.col(0) * faceBeta2Bounds->halflengthX());
  tTransform
      = new Transform3D(beta2Rotation * Translation3D(faceBeta2Position));
  rSurfaces.push_back(Surface::makeShared<PlaneSurface>(
      std::shared_ptr<const Transform3D>(tTransform),
      std::shared_ptr<const PlanarBounds>(faceBeta2Bounds)));
  // face surfaces zx
  //   (7) - at negative local y
  tTransform
      = new Transform3D(transform * AngleAxis3D(M_PI, Vector3D(1., 0., 0.))
                        * Translation3D(Vector3D(0., 2 * halflengthY1(), 0.))
                        * AngleAxis3D(-0.5 * M_PI, Vector3D(0., 1., 0.))
                        * AngleAxis3D(-0.5 * M_PI, Vector3D(1., 0., 0.)));
  rSurfaces.push_back(Surface::makeShared<PlaneSurface>(
      std::shared_ptr<const Transform3D>(tTransform),
      std::shared_ptr<const PlanarBounds>(faceZXRectangleBoundsBottom())));
  //   (8) - at positive local y
  tTransform = new Transform3D(
      transform * Translation3D(Vector3D(0., halflengthY2(), 0.))
      * AngleAxis3D(-0.5 * M_PI, Vector3D(0., 1., 0.))
      * AngleAxis3D(-0.5 * M_PI, Vector3D(1., 0., 0.)));
  rSurfaces.push_back(Surface::makeShared<PlaneSurface>(
      std::shared_ptr<const Transform3D>(tTransform),
      std::shared_ptr<const PlanarBounds>(faceZXRectangleBoundsTop())));

  return rSurfaces;
}

// faces in xy
Acts::DiamondBounds*
Acts::DoubleTrapezoidVolumeBounds::faceXYDiamondBounds() const
{
  return new DiamondBounds(m_valueStore.at(bv_minHalfX),
                           m_valueStore.at(bv_medHalfX),
                           m_valueStore.at(bv_maxHalfX),
                           m_valueStore.at(bv_halfY1),
                           m_valueStore.at(bv_halfY2));
}

Acts::RectangleBounds*
Acts::DoubleTrapezoidVolumeBounds::faceAlpha1RectangleBounds() const
{
  return new RectangleBounds(m_valueStore.at(bv_halfY1)
                                 / cos(m_valueStore.at(bv_alpha1)),
                             m_valueStore.at(bv_halfZ));
}

Acts::RectangleBounds*
Acts::DoubleTrapezoidVolumeBounds::faceAlpha2RectangleBounds() const
{
  return new RectangleBounds(m_valueStore.at(bv_halfY2)
                                 / cos(m_valueStore.at(bv_alpha2)),
                             m_valueStore.at(bv_halfZ));
}

Acts::RectangleBounds*
Acts::DoubleTrapezoidVolumeBounds::faceBeta1RectangleBounds() const
{
  return new RectangleBounds(m_valueStore.at(bv_halfY1)
                                 / cos(m_valueStore.at(bv_alpha1)),
                             m_valueStore.at(bv_halfZ));
}

Acts::RectangleBounds*
Acts::DoubleTrapezoidVolumeBounds::faceBeta2RectangleBounds() const
{
  return new RectangleBounds(m_valueStore.at(bv_halfY2)
                                 / cos(m_valueStore.at(bv_alpha2)),
                             m_valueStore.at(bv_halfZ));
}

Acts::RectangleBounds*
Acts::DoubleTrapezoidVolumeBounds::faceZXRectangleBoundsBottom() const
{
  return new RectangleBounds(m_valueStore.at(bv_halfZ),
                             m_valueStore.at(bv_minHalfX));
}

Acts::RectangleBounds*
Acts::DoubleTrapezoidVolumeBounds::faceZXRectangleBoundsTop() const
{
  return new RectangleBounds(m_valueStore.at(bv_halfZ),
                             m_valueStore.at(bv_maxHalfX));
}

bool
Acts::DoubleTrapezoidVolumeBounds::inside(const Vector3D& pos, double tol) const
{
  if (std::abs(pos.z()) > m_valueStore.at(bv_halfZ) + tol) {
    return false;
  }
  if (pos.y() < -2 * m_valueStore.at(bv_halfY1) - tol) {
    return false;
  }
  if (pos.y() > 2 * m_valueStore.at(bv_halfY2) - tol) {
    return false;
  }
  DiamondBounds* faceXYBounds = faceXYDiamondBounds();
  Vector2D       locp(pos.x(), pos.y());
  bool inside(faceXYBounds->inside(locp, BoundaryCheck(true, true, tol, tol)));
  delete faceXYBounds;
  return inside;
}

// ostream operator overload
std::ostream&
Acts::DoubleTrapezoidVolumeBounds::dump(std::ostream& sl) const
{
  return dumpT<std::ostream>(sl);
}
