// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// TrapezoidBounds.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Surfaces/TrapezoidBounds.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>

Acts::TrapezoidBounds::TrapezoidBounds(double minhalex,
                                       double maxhalex,
                                       double haley)
  : PlanarBounds(TrapezoidBounds::bv_length)
  , m_alpha(0.)
  , m_beta(0.)
  , m_boundingBox(0.,0.)
{
  m_valueStore.at(TrapezoidBounds::bv_minHalfX) = std::abs(minhalex);
  m_valueStore.at(TrapezoidBounds::bv_maxHalfX) = std::abs(maxhalex);
  m_valueStore.at(TrapezoidBounds::bv_halfY)    = std::abs(haley);
  // find the maximum at for the bounding box  
  double mx = m_valueStore.at(TrapezoidBounds::bv_minHalfX) > 
              m_valueStore.at(TrapezoidBounds::bv_maxHalfX) ?
              m_valueStore.at(TrapezoidBounds::bv_minHalfX) :
              m_valueStore.at(TrapezoidBounds::bv_maxHalfX);
  m_boundingBox = RectangleBounds(mx, m_valueStore.at(TrapezoidBounds::bv_halfY));
}

Acts::TrapezoidBounds::TrapezoidBounds(double minhalex,
                                       double haley,
                                       double alpha,
                                       double beta)
  : PlanarBounds(TrapezoidBounds::bv_length)
  , m_alpha(alpha)
  , m_beta(beta)
  , m_boundingBox(0.,0.)
{
  double gamma = (alpha > beta) ? (alpha - 0.5 * M_PI) : (beta - 0.5 * M_PI);
  // now fill them
  m_valueStore.at(TrapezoidBounds::bv_minHalfX) = std::abs(minhalex);
  m_valueStore.at(TrapezoidBounds::bv_maxHalfX) = minhalex
      + (2. * m_valueStore.at(TrapezoidBounds::bv_halfY)) * tan(gamma);
  m_valueStore.at(TrapezoidBounds::bv_halfY) = std::abs(haley);
  // find the maximum for the bounding box
  double mx = m_valueStore.at(TrapezoidBounds::bv_minHalfX) > 
              m_valueStore.at(TrapezoidBounds::bv_maxHalfX) ?
              m_valueStore.at(TrapezoidBounds::bv_minHalfX) :
              m_valueStore.at(TrapezoidBounds::bv_maxHalfX);
  m_boundingBox = RectangleBounds(mx, m_valueStore.at(TrapezoidBounds::bv_halfY));
}

Acts::TrapezoidBounds::~TrapezoidBounds()
{
}

Acts::TrapezoidBounds&
Acts::TrapezoidBounds::operator=(const TrapezoidBounds& trabo)
{
  if (this != &trabo) {
    PlanarBounds::operator=(trabo);
    m_alpha               = trabo.m_alpha;
    m_beta                = trabo.m_beta;
    m_boundingBox         = trabo.m_boundingBox;
  }
  return *this;
}

bool
Acts::TrapezoidBounds::inside(const Acts::Vector2D& lpos,
                              double                tol0,
                              double                tol1) const
{
  if (m_alpha == 0.) return insideFull(lpos, tol0, tol1);
  return (insideFull(lpos, tol0, tol1) && !insideExclude(lpos, tol0, tol1));
}

bool
Acts::TrapezoidBounds::insideFull(const Acts::Vector2D& lpos,
                                  double                tol0,
                                  double                tol1) const
{
  // the cases:
  // the cases:
  double std::absX = std::abs(lpos[Acts::eLOC_X]);
  double std::absY = std::abs(lpos[Acts::eLOC_Y]);
  // (1) a fast FALSE
  if (std::absY > (m_valueStore.at(TrapezoidBounds::bv_halfY) + tol1)) return false;
  // (2) a fast FALSE
  if (std::absX > (m_valueStore.at(TrapezoidBounds::bv_maxHalfX) + tol0))
    return false;
  // (3) a fast TRUE
  if (std::absX < (m_valueStore.at(TrapezoidBounds::bv_minHalfX) - tol0))
    return true;
  // (4) particular case - a rectangle
  if (m_valueStore.at(TrapezoidBounds::bv_maxHalfX)
      == m_valueStore.at(TrapezoidBounds::bv_minHalfX))
    return true;
  // (5) /** use basic calculation of a straight line */
  double k = 2.0 * m_valueStore.at(TrapezoidBounds::bv_halfY)
      / (m_valueStore.at(TrapezoidBounds::bv_maxHalfX)
         - m_valueStore.at(TrapezoidBounds::bv_minHalfX))
      * ((lpos[Acts::eLOC_X] > 0.) ? 1.0 : -1.0);
  double d
      = -std::abs(k) * 0.5 * (m_valueStore.at(TrapezoidBounds::bv_maxHalfX)
                              + m_valueStore.at(TrapezoidBounds::bv_minHalfX));
  return (isAbove(lpos, tol0, tol1, k, d));
}

bool
Acts::TrapezoidBounds::insideExclude(const Acts::Vector2D& lpos,
                                     double                tol0,
                                     double                tol1) const
{
  // line a
  bool   alphaBiggerBeta(m_alpha > m_beta);
  double ka   = alphaBiggerBeta ? tan(M_PI - m_alpha) : tan(m_alpha);
  double kb   = alphaBiggerBeta ? tan(M_PI - m_beta) : tan(m_beta);
  double sign = alphaBiggerBeta ? -1. : 1.;
  double da   = -m_valueStore.at(TrapezoidBounds::bv_halfY)
      + sign * ka * m_valueStore.at(TrapezoidBounds::bv_minHalfX);
  double db = -m_valueStore.at(TrapezoidBounds::bv_halfY)
      + sign * kb * m_valueStore.at(TrapezoidBounds::bv_minHalfX);

  return (isAbove(lpos, tol0, tol1, ka, da)
          && isAbove(lpos, tol0, tol1, kb, db));
}

bool
Acts::TrapezoidBounds::isAbove(const Acts::Vector2D& lpos,
                               double                tol0,
                               double                tol1,
                               double                k,
                               double                d) const
{
  // the most tolerant approach for tol0 and tol1
  double sign = k > 0. ? -1. : +1.;
  return (lpos[Acts::eLOC_Y] + tol1
          > (k * (lpos[Acts::eLOC_X] + sign * tol0) + d));
}

double
Acts::TrapezoidBounds::distanceToBoundary(const Acts::Vector2D& lpos) const
{
  const int Np = 4;

  double xl = -m_valueStore.at(TrapezoidBounds::bv_maxHalfX);
  double xr = m_valueStore.at(TrapezoidBounds::bv_maxHalfX);
  if (m_alpha != 0.) {
    xl = -m_valueStore.at(TrapezoidBounds::bv_minHalfX)
        - 2. * tan(m_alpha) * m_valueStore.at(TrapezoidBounds::bv_halfY);
  } else if (m_beta != 0.) {
    xr = m_valueStore.at(TrapezoidBounds::bv_minHalfX)
        + 2. * tan(m_beta) * m_valueStore.at(TrapezoidBounds::bv_halfY);
  }
  double X[4] = {-m_valueStore.at(TrapezoidBounds::bv_minHalfX),
                 xl,
                 xr,
                 m_valueStore.at(TrapezoidBounds::bv_minHalfX)};
  double Y[4] = {-m_valueStore.at(TrapezoidBounds::bv_halfY),
                 m_valueStore.at(TrapezoidBounds::bv_halfY),
                 m_valueStore.at(TrapezoidBounds::bv_halfY),
                 -m_valueStore.at(TrapezoidBounds::bv_halfY)};

  double dm = 1.e+20;
  double Ao = 0.;
  bool   in = true;

  for (int i = 0; i != Np; ++i) {
    int j          = i + 1;
    if (j == Np) j = 0;

    double x  = X[i] - lpos[0];
    double y  = Y[i] - lpos[1];
    double dx = X[j] - X[i];
    double dy = Y[j] - Y[i];
    double A  = x * dy - y * dx;
    double S  = -(x * dx + y * dy);

    if (S <= 0.) {
      double d       = x * x + y * y;
      if (d < dm) dm = d;
    } else {
      double a = dx * dx + dy * dy;
      if (S <= a) {
        double d       = (A * A) / a;
        if (d < dm) dm = d;
      }
    }
    if (i && in && Ao * A < 0.) in = false;
    Ao                             = A;
  }
  if (in)
    return -sqrt(dm);
  else
    return sqrt(dm);
}

std::ostream&
Acts::TrapezoidBounds::dump(std::ostream& sl) const
{
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::TrapezoidBounds:  (minHlenghtX, maxHlengthX, hlengthY) = "
     << "(" << m_valueStore.at(TrapezoidBounds::bv_minHalfX) << ", "
     << m_valueStore.at(TrapezoidBounds::bv_maxHalfX) << ", "
     << m_valueStore.at(TrapezoidBounds::bv_halfY) << ")";
  sl << std::setprecision(-1);
  return sl;
}
