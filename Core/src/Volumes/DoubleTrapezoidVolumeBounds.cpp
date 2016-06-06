// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// DoubleTrapezoidVolumeBounds.cpp, ACTS project
///////////////////////////////////////////////////////////////////

// Geometry module
#include "ACTS/Volumes/DoubleTrapezoidVolumeBounds.hpp"
#include "ACTS/Surfaces/DiamondBounds.hpp"
#include "ACTS/Surfaces/PlaneSurface.hpp"
#include "ACTS/Surfaces/RectangleBounds.hpp"
// STD/STL
#include <iomanip>
#include <iostream>
#include <math.h>

Acts::DoubleTrapezoidVolumeBounds::DoubleTrapezoidVolumeBounds()
  : VolumeBounds(), m_boundValues(bv_length, 0.)
{
}

Acts::DoubleTrapezoidVolumeBounds::DoubleTrapezoidVolumeBounds(double minhalex,
                                                               double medhalex,
                                                               double maxhalex,
                                                               double haley1,
                                                               double haley2,
                                                               double halez)
  : VolumeBounds(), m_boundValues(bv_length, 0.)
{
  m_boundValues.at(bv_minHalfX) = minhalex;
  m_boundValues.at(bv_medHalfX) = medhalex;
  m_boundValues.at(bv_maxHalfX) = maxhalex;
  m_boundValues.at(bv_halfY1)   = haley1;
  m_boundValues.at(bv_halfY2)   = haley2;
  m_boundValues.at(bv_halfZ)    = halez;
  m_boundValues.at(bv_alpha1)
      = atan2(m_boundValues.at(bv_medHalfX) - m_boundValues.at(bv_minHalfX),
              2. * m_boundValues.at(bv_halfY1));
  m_boundValues.at(bv_alpha2)
      = atan2(m_boundValues.at(bv_medHalfX) - m_boundValues.at(bv_maxHalfX),
              2. * m_boundValues.at(bv_halfY2));
}

Acts::DoubleTrapezoidVolumeBounds::DoubleTrapezoidVolumeBounds(
    const Acts::DoubleTrapezoidVolumeBounds& trabo)
  : VolumeBounds(), m_boundValues(trabo.m_boundValues)
{
}

Acts::DoubleTrapezoidVolumeBounds::~DoubleTrapezoidVolumeBounds()
{
}

Acts::DoubleTrapezoidVolumeBounds&
Acts::DoubleTrapezoidVolumeBounds::
operator=(const Acts::DoubleTrapezoidVolumeBounds& trabo)
{
  if (this != &trabo) m_boundValues = trabo.m_boundValues;
  return *this;
}

const std::vector<const Acts::Surface*>*
Acts::DoubleTrapezoidVolumeBounds::decomposeToSurfaces(
    std::shared_ptr<Acts::Transform3D> transformPtr) const
{
  std::vector<const Acts::Surface*>* retsf
      = new std::vector<const Acts::Surface*>;

  // the transform
  Acts::Transform3D transform = (transformPtr == nullptr)
      ? Acts::Transform3D::Identity()
      : (*transformPtr.get());

  // face surfaces xy
  Acts::RotationMatrix3D diamondRotation(transform.rotation());
  Acts::Vector3D         diamondX(diamondRotation.col(0));
  Acts::Vector3D         diamondY(diamondRotation.col(1));
  Acts::Vector3D         diamondZ(diamondRotation.col(2));
  Acts::Vector3D         diamondCenter(transform.translation());

  Acts::Transform3D* tTransform = nullptr;

  //   (1) - at negative local z
  tTransform = new Acts::Transform3D(
      transform * Acts::AngleAxis3D(M_PI, Acts::Vector3D(0., 1., 0.))
      * Acts::Translation3D(Acts::Vector3D(0., 0., halflengthZ())));
  retsf->push_back(new Acts::PlaneSurface(
      std::shared_ptr<Acts::Transform3D>(tTransform), faceXYDiamondBounds()));
  //   (2) - at positive local z
  tTransform = new Acts::Transform3D(
      transform * Acts::Translation3D(Acts::Vector3D(0., 0., halflengthZ())));
  retsf->push_back(new Acts::PlaneSurface(
      std::shared_ptr<Acts::Transform3D>(tTransform), faceXYDiamondBounds()));
  // face surfaces yz
  // transmute cyclical
  //   (3) - at point A, attached to alpha opening angle
  Acts::Vector3D A(diamondCenter - minHalflengthX() * diamondX
                   - 2 * halflengthY1() * diamondY);
  Acts::AngleAxis3D      alpha1ZRotation(alpha1(), Acts::Vector3D(0., 0., 1.));
  Acts::RotationMatrix3D alpha1Rotation(
      diamondRotation * alpha1ZRotation
      * Acts::AngleAxis3D(-0.5 * M_PI, Acts::Vector3D(0., 1., 0.))
      * Acts::AngleAxis3D(0.5 * M_PI, Acts::Vector3D(0., 0., 1.)));
  RectangleBounds* faceAlpha1Bounds = faceAlpha1RectangleBounds();
  Acts::Vector3D   faceAlpha1Position(
      A + alpha1Rotation.col(0) * faceAlpha1Bounds->halflengthX());
  tTransform = new Acts::Transform3D(alpha1Rotation
                                     * Acts::Translation3D(faceAlpha1Position));
  retsf->push_back(new Acts::PlaneSurface(
      std::shared_ptr<Acts::Transform3D>(tTransform), faceAlpha1Bounds));
  //   (4) - at point B, attached to beta opening angle
  Acts::Vector3D B(diamondCenter + minHalflengthX() * diamondX
                   - 2 * halflengthY1() * diamondY);
  Acts::AngleAxis3D      beta1ZRotation(-alpha1(), Acts::Vector3D(0., 0., 1.));
  Acts::RotationMatrix3D beta1Rotation(
      diamondRotation * beta1ZRotation
      * Acts::AngleAxis3D(0.5 * M_PI, Acts::Vector3D(0., 1., 0.))
      * Acts::AngleAxis3D(0.5 * M_PI, Acts::Vector3D(0., 0., 1.)));
  RectangleBounds* faceBeta1Bounds = faceBeta1RectangleBounds();
  Acts::Vector3D   faceBeta1Position(
      B + beta1Rotation.col(0) * faceBeta1Bounds->halflengthX());
  tTransform = new Acts::Transform3D(beta1Rotation
                                     * Acts::Translation3D(faceBeta1Position));
  retsf->push_back(new Acts::PlaneSurface(
      std::shared_ptr<Acts::Transform3D>(tTransform), faceBeta1Bounds));
  // face surfaces yz
  // transmute cyclical
  //   (5) - at point A', attached to alpha opening angle
  Acts::Vector3D AA(diamondCenter - maxHalflengthX() * diamondX
                    + 2 * halflengthY2() * diamondY);
  Acts::AngleAxis3D      alpha2ZRotation(-alpha2(), Acts::Vector3D(0., 0., 1.));
  Acts::RotationMatrix3D alpha2Rotation(
      diamondRotation * alpha2ZRotation
      * Acts::AngleAxis3D(-0.5 * M_PI, Acts::Vector3D(0., 1., 0.))
      * Acts::AngleAxis3D(-0.5 * M_PI, Acts::Vector3D(0., 0., 1.)));
  RectangleBounds* faceAlpha2Bounds = faceAlpha2RectangleBounds();
  Acts::Vector3D   faceAlpha2Position(
      AA + alpha2Rotation.col(0) * faceAlpha2Bounds->halflengthX());
  tTransform = new Acts::Transform3D(alpha2Rotation
                                     * Acts::Translation3D(faceAlpha2Position));
  retsf->push_back(new Acts::PlaneSurface(
      std::shared_ptr<Acts::Transform3D>(tTransform), faceAlpha2Bounds));
  //   (6) - at point B', attached to beta opening angle
  Acts::Vector3D BB(diamondCenter + maxHalflengthX() * diamondX
                    + 2 * halflengthY2() * diamondY);
  Acts::AngleAxis3D      beta2ZRotation(alpha2(), Acts::Vector3D(0., 0., 1.));
  Acts::RotationMatrix3D beta2Rotation(
      diamondRotation * beta2ZRotation
      * Acts::AngleAxis3D(0.5 * M_PI, Acts::Vector3D(0., 1., 0.))
      * Acts::AngleAxis3D(-0.5 * M_PI, Acts::Vector3D(0., 0., 1.)));
  RectangleBounds* faceBeta2Bounds = faceBeta2RectangleBounds();
  Acts::Vector3D   faceBeta2Position(
      BB + beta2Rotation.col(0) * faceBeta2Bounds->halflengthX());
  tTransform = new Acts::Transform3D(beta2Rotation
                                     * Acts::Translation3D(faceBeta2Position));
  retsf->push_back(new Acts::PlaneSurface(
      std::shared_ptr<Acts::Transform3D>(tTransform), faceBeta2Bounds));
  // face surfaces zx
  //   (7) - at negative local y
  tTransform = new Acts::Transform3D(
      transform * Acts::AngleAxis3D(M_PI, Acts::Vector3D(1., 0., 0.))
      * Acts::Translation3D(Acts::Vector3D(0., 2 * halflengthY1(), 0.))
      * Acts::AngleAxis3D(-0.5 * M_PI, Acts::Vector3D(0., 1., 0.))
      * Acts::AngleAxis3D(-0.5 * M_PI, Acts::Vector3D(1., 0., 0.)));
  retsf->push_back(
      new Acts::PlaneSurface(std::shared_ptr<Acts::Transform3D>(tTransform),
                             faceZXRectangleBoundsBottom()));
  //   (8) - at positive local y
  tTransform = new Acts::Transform3D(
      transform * Acts::Translation3D(Acts::Vector3D(0., halflengthY2(), 0.))
      * Acts::AngleAxis3D(-0.5 * M_PI, Acts::Vector3D(0., 1., 0.))
      * Acts::AngleAxis3D(-0.5 * M_PI, Acts::Vector3D(1., 0., 0.)));
  retsf->push_back(
      new Acts::PlaneSurface(std::shared_ptr<Acts::Transform3D>(tTransform),
                             faceZXRectangleBoundsTop()));

  return retsf;
}

// faces in xy
Acts::DiamondBounds*
Acts::DoubleTrapezoidVolumeBounds::faceXYDiamondBounds() const
{
  return new Acts::DiamondBounds(m_boundValues.at(bv_minHalfX),
                                 m_boundValues.at(bv_medHalfX),
                                 m_boundValues.at(bv_maxHalfX),
                                 m_boundValues.at(bv_halfY1),
                                 m_boundValues.at(bv_halfY2));
}

Acts::RectangleBounds*
Acts::DoubleTrapezoidVolumeBounds::faceAlpha1RectangleBounds() const
{
  return new Acts::RectangleBounds(m_boundValues.at(bv_halfY1)
                                       / cos(m_boundValues.at(bv_alpha1)),
                                   m_boundValues.at(bv_halfZ));
}

Acts::RectangleBounds*
Acts::DoubleTrapezoidVolumeBounds::faceAlpha2RectangleBounds() const
{
  return new Acts::RectangleBounds(m_boundValues.at(bv_halfY2)
                                       / cos(m_boundValues.at(bv_alpha2)),
                                   m_boundValues.at(bv_halfZ));
}

Acts::RectangleBounds*
Acts::DoubleTrapezoidVolumeBounds::faceBeta1RectangleBounds() const
{
  return new Acts::RectangleBounds(m_boundValues.at(bv_halfY1)
                                       / cos(m_boundValues.at(bv_alpha1)),
                                   m_boundValues.at(bv_halfZ));
}

Acts::RectangleBounds*
Acts::DoubleTrapezoidVolumeBounds::faceBeta2RectangleBounds() const
{
  return new Acts::RectangleBounds(m_boundValues.at(bv_halfY2)
                                       / cos(m_boundValues.at(bv_alpha2)),
                                   m_boundValues.at(bv_halfZ));
}

Acts::RectangleBounds*
Acts::DoubleTrapezoidVolumeBounds::faceZXRectangleBoundsBottom() const
{
  return new Acts::RectangleBounds(m_boundValues.at(bv_halfZ),
                                   m_boundValues.at(bv_minHalfX));
}

Acts::RectangleBounds*
Acts::DoubleTrapezoidVolumeBounds::faceZXRectangleBoundsTop() const
{
  return new Acts::RectangleBounds(m_boundValues.at(bv_halfZ),
                                   m_boundValues.at(bv_maxHalfX));
}

bool
Acts::DoubleTrapezoidVolumeBounds::inside(const Acts::Vector3D& pos,
                                          double                tol) const
{
  if (fabs(pos.z()) > m_boundValues.at(bv_halfZ) + tol) return false;
  if (pos.y() < -2 * m_boundValues.at(bv_halfY1) - tol) return false;
  if (pos.y() > 2 * m_boundValues.at(bv_halfY2) - tol) return false;
  Acts::DiamondBounds* faceXYBounds = faceXYDiamondBounds();
  Acts::Vector2D       locp(pos.x(), pos.y());
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
