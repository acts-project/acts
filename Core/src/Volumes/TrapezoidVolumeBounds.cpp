// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// TrapezoidVolumeBounds.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Volumes/TrapezoidVolumeBounds.hpp"
#include "ACTS/Surfaces/PlaneSurface.hpp"
#include "ACTS/Surfaces/RectangleBounds.hpp"
#include "ACTS/Surfaces/TrapezoidBounds.hpp"
#include "ACTS/Utilities/GeometryStatics.hpp"
#include <iomanip>
#include <iostream>
#include <math.h>

Acts::TrapezoidVolumeBounds::TrapezoidVolumeBounds()
  : VolumeBounds()
  , m_valueStore(bv_length, 0.)
{
}

Acts::TrapezoidVolumeBounds::TrapezoidVolumeBounds(double minhalex,
                                                   double maxhalex,
                                                   double haley,
                                                   double halez)
  : VolumeBounds()
  , m_valueStore(bv_length, 0.)
{
  m_valueStore.at(bv_minHalfX) = minhalex;
  m_valueStore.at(bv_maxHalfX) = maxhalex;
  m_valueStore.at(bv_halfY)    = haley;
  m_valueStore.at(bv_halfZ)    = halez;
  m_valueStore.at(bv_alpha)
      = atan((m_valueStore.at(bv_maxHalfX) - m_valueStore.at(bv_minHalfX)) / 2
             / m_valueStore.at(bv_halfY))
      + 0.5 * M_PI;
  m_valueStore.at(bv_beta) = m_valueStore.at(bv_alpha);
}

Acts::TrapezoidVolumeBounds::TrapezoidVolumeBounds(double minhalex,
                                                   double haley,
                                                   double halez,
                                                   double alpha,
                                                   double beta)
  : VolumeBounds()
  , m_valueStore(bv_length, 0.)
{
  m_valueStore.at(bv_minHalfX) = minhalex;
  m_valueStore.at(bv_halfY)    = haley;
  m_valueStore.at(bv_halfZ)    = halez;
  m_valueStore.at(bv_alpha)    = alpha;
  m_valueStore.at(bv_beta)     = beta;
  // now calculate the remaining max half X
  double gamma = (alpha > beta) ? (alpha - 0.5 * M_PI) : (beta - 0.5 * M_PI);
  m_valueStore.at(bv_maxHalfX) = minhalex + (2. * haley) * tan(gamma);
}

Acts::TrapezoidVolumeBounds::TrapezoidVolumeBounds(
    const TrapezoidVolumeBounds& trabo)
  : VolumeBounds(), m_valueStore(trabo.m_valueStore)
{
}

Acts::TrapezoidVolumeBounds::~TrapezoidVolumeBounds()
{
}

TrapezoidVolumeBounds&
TrapezoidVolumeBounds::operator=(const TrapezoidVolumeBounds& trabo)
{
  if (this != &trabo) m_valueStore = trabo.m_valueStore;
  return *this;
}

const std::vector<const Acts::Surface*>*
Acts::TrapezoidVolumeBounds::decomposeToSurfaces(
    std::shared_ptr<Transform3D> transformPtr) const
{
  std::vector<const Surface*> rSurfaces;

  Transform3D transform = (transformPtr == nullptr)
      ? Transform3D::Identity()
      : (*transformPtr.get());
  Transform3D* tTransform = nullptr;
  // face surfaces xy
  RotationMatrix3D trapezoidRotation(transform.rotation());
  Vector3D         trapezoidX(trapezoidRotation.col(0));
  Vector3D         trapezoidY(trapezoidRotation.col(1));
  Vector3D         trapezoidZ(trapezoidRotation.col(2));
  Vector3D         trapezoidCenter(transform.translation());

  //   (1) - at negative local z
  std::shared_ptr<const PlanarBounds> xytBounds(faceXYTrapezoidBounds());
  tTransform = new Transform3D(
      transform * AngleAxis3D(M_PI, Vector3D(0., 1., 0.))
      * Translation3D(Vector3D(0., 0., halflengthZ())));

  rSurfaces.push_back(new PlaneSurface(
      std::shared_ptr<Transform3D>(tTransform), xytBounds));
  //   (2) - at positive local z
  tTransform = new Transform3D(
      transform * Translation3D(Vector3D(0., 0., halflengthZ())));
  rSurfaces.push_back(new PlaneSurface(
      std::shared_ptr<Transform3D>(tTransform), xytBounds));
  // face surfaces yz
  // transmute cyclical
  //   (3) - at point A, attached to alpha opening angle
  Vector3D A(minHalflengthX(), halflengthY(), trapezoidCenter.z());
  RotationMatrix3D alphaZRotation
      = (s_idRotation
         * AngleAxis3D(alpha() - 0.5 * M_PI, Vector3D(0., 0., 1.)))
            .toRotationMatrix();
  // CLHEP::HepRotation  alphaRotation(alphaZRotation*trapezoidRotation);
  RotationMatrix3D faceAlphaRotation;
  faceAlphaRotation.col(0)         = alphaZRotation.col(1);
  faceAlphaRotation.col(1)         = -alphaZRotation.col(2);
  faceAlphaRotation.col(2)         = -alphaZRotation.col(0);
  RectangleBounds* faceAlphaBounds = faceAlphaRectangleBounds();
  // Vector3D
  // faceAlphaPosition(A+faceAlphaRotation.colX()*faceAlphaBounds->halflengthX());
  Vector3D faceAlphaPosition0(
      -0.5 * (minHalflengthX() + maxHalflengthX()), 0., 0.);
  Vector3D faceAlphaPosition = transform * faceAlphaPosition0;
  tTransform = new Transform3D((trapezoidRotation * faceAlphaRotation)
                                     * Translation3D(faceAlphaPosition));
  rSurfaces.push_back(new PlaneSurface(
      std::shared_ptr<Transform3D>(tTransform), faceAlphaBounds));
  //   (4) - at point B, attached to beta opening angle
  Vector3D B(minHalflengthX(), -halflengthY(), trapezoidCenter.z());
  RotationMatrix3D betaZRotation
      = (s_idRotation * AngleAxis3D(-(beta() - 0.5 * M_PI),
                                          Vector3D(0., 0., 1.)))
            .toRotationMatrix();
  // CLHEP::HepRotation  betaRotation(betaZRotation*trapezoidRotation);
  RotationMatrix3D faceBetaRotation;
  faceBetaRotation.col(0)         = betaZRotation.col(1);
  faceBetaRotation.col(1)         = betaZRotation.col(2);
  faceBetaRotation.col(2)         = betaZRotation.col(0);
  RectangleBounds* faceBetaBounds = faceBetaRectangleBounds();
  // Vector3D
  // faceBetaPosition(B+faceBetaRotation.colX()*faceBetaBounds->halflengthX());
  Vector3D faceBetaPosition0(
      0.5 * (minHalflengthX() + maxHalflengthX()), 0., 0.);
  Vector3D faceBetaPosition = transform * faceBetaPosition0;
  tTransform = new Transform3D(trapezoidRotation * faceBetaRotation
                                     * Translation3D(faceBetaPosition));
  rSurfaces.push_back(new PlaneSurface(
      std::shared_ptr<Transform3D>(tTransform), faceBetaBounds));
  // face surfaces zx
  //   (5) - at negative local x
  tTransform = new Transform3D(
      transform * AngleAxis3D(M_PI, Vector3D(1., 0., 0.))
      * Translation3D(Vector3D(0., halflengthY(), 0.))
      * AngleAxis3D(-0.5 * M_PI, Vector3D(0., 1., 0.))
      * AngleAxis3D(-0.5 * M_PI, Vector3D(1., 0., 0.)));
  rSurfaces.push_back(
      new PlaneSurface(std::shared_ptr<Transform3D>(tTransform),
                             faceZXRectangleBoundsBottom()));
  //   (6) - at positive local x
  tTransform = new Transform3D(
      transform * Translation3D(Vector3D(0., halflengthY(), 0.))
      * AngleAxis3D(-0.5 * M_PI, Vector3D(0., 1., 0.))
      * AngleAxis3D(-0.5 * M_PI, Vector3D(1., 0., 0.)));
  rSurfaces.push_back(
      new PlaneSurface(std::shared_ptr<Transform3D>(tTransform),
                             faceZXRectangleBoundsTop()));

  return rSurfaces;
}

TrapezoidBounds*
TrapezoidVolumeBounds::faceXYTrapezoidBounds() const
{
  return new TrapezoidBounds(m_valueStore.at(bv_minHalfX),
                                   m_valueStore.at(bv_maxHalfX),
                                   m_valueStore.at(bv_halfY));
}

RectangleBounds*
TrapezoidVolumeBounds::faceAlphaRectangleBounds() const
{
  return new RectangleBounds(
      m_valueStore.at(bv_halfY) / cos(m_valueStore.at(bv_alpha) - 0.5 * M_PI),
      m_valueStore.at(bv_halfZ));
}

RectangleBounds*
TrapezoidVolumeBounds::faceBetaRectangleBounds() const
{
  return new RectangleBounds(
      m_valueStore.at(bv_halfY) / cos(m_valueStore.at(bv_beta) - 0.5 * M_PI),
      m_valueStore.at(bv_halfZ));
}

RectangleBounds*
TrapezoidVolumeBounds::faceZXRectangleBoundsBottom() const
{
  return new RectangleBounds(m_valueStore.at(bv_halfZ),
                                   m_valueStore.at(bv_minHalfX));
}

RectangleBounds*
TrapezoidVolumeBounds::faceZXRectangleBoundsTop() const
{
  // double delta = (m_valueStore.at(bv_alpha) < m_valueStore.at(bv_beta)) ?
  // m_valueStore.at(bv_alpha) - M_PI/2. : m_valueStore.at(bv_beta) - M_PI/2.;
  // return new RectangleBounds(m_valueStore.at(bv_halfZ),
  // 0.5*(m_valueStore.at(bv_minHalfX)+m_valueStore.at(bv_minHalfX)+2.*m_valueStore.at(bv_halfY)/cos(delta)));
  return new RectangleBounds(m_valueStore.at(bv_halfZ),
                             m_valueStore.at(bv_maxHalfX));
}

bool
TrapezoidVolumeBounds::inside(const Vector3D& pos, double tol) const
{
  if (fabs(pos.z()) > m_valueStore.at(bv_halfZ) + tol) return false;
  if (fabs(pos.y()) > m_valueStore.at(bv_halfY) + tol) return false;
  TrapezoidBounds* faceXYBounds = faceXYTrapezoidBounds();
  Vector2D         locp(pos.x(), pos.y());
  bool inside(faceXYBounds->inside(locp, BoundaryCheck(true, true, tol, tol)));
  delete faceXYBounds;
  return inside;
}

std::ostream&
TrapezoidVolumeBounds::dump(std::ostream& sl) const
{
  return dumpT<std::ostream>(sl);
}
