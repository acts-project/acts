// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// EllipseBounds.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Surfaces/EllipseBounds.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>

Acts::EllipseBounds::EllipseBounds(double minradX,
                                   double minradY,
                                   double maxradX,
                                   double maxradY,
                                   double avephi,
                                   double hphisec)
  : PlanarBounds(EllipseBounds::bv_length), m_boundingBox(0., 0.)
{
  m_valueStore.at(EllipseBounds::bv_rMinX)         = minradX;
  m_valueStore.at(EllipseBounds::bv_rMinY)         = minradY;
  m_valueStore.at(EllipseBounds::bv_rMaxX)         = maxradX;
  m_valueStore.at(EllipseBounds::bv_rMaxY)         = maxradY;
  m_valueStore.at(EllipseBounds::bv_averagePhi)    = avephi;
  m_valueStore.at(EllipseBounds::bv_halfPhiSector) = hphisec;
  double mx     = minradX > maxradX ? minradX : maxradX;
  double my     = minradY > maxradY ? minradY : maxradY;
  m_boundingBox = RectangleBounds(mx, my);
}

Acts::EllipseBounds::~EllipseBounds()
{
}

Acts::EllipseBounds&
Acts::EllipseBounds::operator=(const EllipseBounds& ebo)
{
  if (this != &ebo) {
    PlanarBounds::operator=(ebo);
    m_boundingBox         = ebo.m_boundingBox;
  }
  return *this;
}

// For ellipse bound this is only approximation which is valid
// only if m_valueStore.at(EllipseBounds::bv_rMinX) ~=
// m_valueStore.at(EllipseBounds::bv_rMinY)
// and m_valueStore.at(EllipseBounds::bv_rMaxX) ~=
// m_valueStore.at(EllipseBounds::bv_rMaxY)
//
double
Acts::EllipseBounds::distanceToBoundary(const Vector2D& lpos) const
{
  const double pi2 = 2. * M_PI;

  double r = sqrt(lpos[0] * lpos[0] + lpos[1] * lpos[1]);
  if (r == 0.) {
    if (m_valueStore.at(EllipseBounds::bv_rMinX)
        <= m_valueStore.at(EllipseBounds::bv_rMinY))
      return m_valueStore.at(EllipseBounds::bv_rMinX);
    return m_valueStore.at(EllipseBounds::bv_rMinY);
  }

  const double inv_r = 1. / r;
  double       sn    = lpos[1] * inv_r;
  double       cs    = lpos[0] * inv_r;
  double       sf    = 0.;
  double       dF    = 0.;

  if (m_valueStore.at(EllipseBounds::bv_halfPhiSector) < M_PI) {
    dF = atan2(cs, sn) - m_valueStore.at(EllipseBounds::bv_averagePhi);
    dF += (dF > M_PI) ? -pi2 : (dF < -M_PI) ? pi2 : 0;
    double df = std::abs(dF) - m_valueStore.at(EllipseBounds::bv_halfPhiSector);
    sf        = r * sin(df);
    if (df > 0.) r *= cos(df);
  } else {
    sf = -1.e+10;
  }

  if (sf <= 0.) {
    double a   = cs / m_valueStore.at(EllipseBounds::bv_rMaxX);
    double b   = sn / m_valueStore.at(EllipseBounds::bv_rMaxY);
    double sr0 = r - 1. / sqrt(a * a + b * b);
    if (sr0 >= 0.) return sr0;
    a          = cs / m_valueStore.at(EllipseBounds::bv_rMinX);
    b          = sn / m_valueStore.at(EllipseBounds::bv_rMinY);
    double sr1 = 1. / sqrt(a * a + b * b) - r;
    if (sr1 >= 0.) return sr1;
    if (sf < sr0) sf = sr0;
    if (sf < sr1) sf = sr1;
    return sf;
  }

  double fb;
  fb = (dF > 0.)
      ? m_valueStore.at(EllipseBounds::bv_averagePhi)
          + m_valueStore.at(EllipseBounds::bv_halfPhiSector)
      : m_valueStore.at(EllipseBounds::bv_averagePhi)
          - m_valueStore.at(EllipseBounds::bv_halfPhiSector);
  sn         = sin(fb);
  cs         = cos(fb);
  double a   = cs / m_valueStore.at(EllipseBounds::bv_rMaxX);
  double b   = sn / m_valueStore.at(EllipseBounds::bv_rMaxY);
  double sr0 = r - 1. / sqrt(a * a + b * b);
  if (sr0 >= 0.) return sqrt(sr0 * sr0 + sf * sf);
  a          = cs / m_valueStore.at(EllipseBounds::bv_rMinX);
  b          = sn / m_valueStore.at(EllipseBounds::bv_rMinY);
  double sr1 = 1. / sqrt(a * a + b * b) - r;
  if (sr1 >= 0.) return sqrt(sr1 * sr1 + sf * sf);
  return sf;
}

// ostream operator overload
std::ostream&
Acts::EllipseBounds::dump(std::ostream& sl) const
{
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::EllipseBounds:  (innerRadiusX, innerRadiusY, outerRadiusX, "
        "outerRadiusY, hPhiSector) = ";
  sl << "(" << rMinX() << ", " << rMinY() << ", " << rMaxX() << ", " << rMaxY()
     << ", " << averagePhi() << ", " << halfPhiSector() << ")";
  sl << std::setprecision(-1);
  return sl;
}
