// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// DiamondBounds.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Surfaces/DiamondBounds.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>

Acts::DiamondBounds::DiamondBounds(double minhalex,
                                   double medhalex,
                                   double maxhalex,
                                   double haley1,
                                   double haley2)
  : PlanarBounds(DiamondBounds::bv_length), m_alpha1(0.), m_alpha2(0.)
{
  m_valueStore.at(DiamondBounds::bv_minHalfX) = minhalex;
  m_valueStore.at(DiamondBounds::bv_medHalfX) = medhalex;
  m_valueStore.at(DiamondBounds::bv_maxHalfX) = maxhalex;
  m_valueStore.at(DiamondBounds::bv_halfY1)   = haley1;
  m_valueStore.at(DiamondBounds::bv_halfY2)   = haley2;
  if (minhalex > maxhalex)
    std::swap(m_valueStore.at(DiamondBounds::bv_minHalfX),
              m_valueStore.at(DiamondBounds::bv_maxHalfX));
  initCache();
}

// copy constructor
Acts::DiamondBounds::DiamondBounds(const DiamondBounds& diabo)
  : PlanarBounds(DiamondBounds::bv_length)
  , m_valueStore(diabo.m_valueStore)
  , m_alpha1(diabo.m_alpha1)
  , m_alpha2(diabo.m_alpha2)
{
}

// destructor
Acts::DiamondBounds::~DiamondBounds()
{
}

Acts::DiamondBounds&
Acts::DiamondBounds::operator=(const DiamondBounds& diabo)
{
  if (this != &diabo) {
    m_valueStore = diabo.m_valueStore;
    m_alpha1     = diabo.m_alpha1;
    m_alpha2     = diabo.m_alpha2;
  }
  return *this;
}

bool
Acts::DiamondBounds::operator==(const Acts::SurfaceBounds& sbo) const
{
  // fast exit
  if (&sbo == this) return true;
  // check the type first not to compare apples with oranges
  const Acts::DiamondBounds* diabo
      = dynamic_cast<const Acts::DiamondBounds*>(&sbo);
  if (!diabo) return false;
  return (m_valueStore == diabo->m_valueStore);
}

// checking if inside bounds (Full symmetrical Diamond)
bool
Acts::DiamondBounds::inside(const Acts::Vector2D& locpo,
                            double                tol1,
                            double                tol2) const
{
  return insideFull(locpo, tol1, tol2);
}

// checking if inside bounds (Full symmetrical Diamond)
bool
Acts::DiamondBounds::insideFull(const Acts::Vector2D& locpo,
                                double                tol1,
                                double                tol2) const
{
  // the cases:
  // (0)
  if (!m_valueStore.at(DiamondBounds::bv_halfY1)
      && !m_valueStore.at(DiamondBounds::bv_minHalfX))
    return false;
  // (1)
  if (locpo[Acts::eLOC_Y]
      < -2. * m_valueStore.at(DiamondBounds::bv_halfY1) - tol2)
    return false;
  if (locpo[Acts::eLOC_Y]
      > 2. * m_valueStore.at(DiamondBounds::bv_halfY2) + tol2)
    return false;
  // (2)
  if (std::abs(locpo[Acts::eLOC_X])
      > (m_valueStore.at(DiamondBounds::bv_medHalfX) + tol1))
    return false;
  // (3)
  if (std::abs(locpo[Acts::eLOC_X])
      < (fmin(m_valueStore.at(DiamondBounds::bv_minHalfX),
              m_valueStore.at(DiamondBounds::bv_maxHalfX))
         - tol1))
    return true;
  // (4)
  /** use basic calculation of a straight line */
  if (locpo[Acts::eLOC_Y] < 0) {
    double k = m_valueStore.at(DiamondBounds::bv_halfY1) > 0.
        ? (m_valueStore.at(DiamondBounds::bv_medHalfX)
           - m_valueStore.at(DiamondBounds::bv_minHalfX))
            / 2 / m_valueStore.at(DiamondBounds::bv_halfY1)
        : 0.;
    return (std::abs(locpo[Acts::eLOC_X])
            <= m_valueStore.at(DiamondBounds::bv_medHalfX)
                - k * std::abs(locpo[Acts::eLOC_Y]));
  } else {
    double k = m_valueStore.at(DiamondBounds::bv_halfY2) > 0.
        ? (m_valueStore.at(DiamondBounds::bv_medHalfX)
           - m_valueStore.at(DiamondBounds::bv_maxHalfX))
            / 2 / m_valueStore.at(DiamondBounds::bv_halfY2)
        : 0.;
    return (std::abs(locpo[Acts::eLOC_X])
            <= m_valueStore.at(DiamondBounds::bv_medHalfX)
                - k * std::abs(locpo[Acts::eLOC_Y]));
  }
}

double
Acts::DiamondBounds::alpha1() const
{
  return m_alpha1;
}

double
Acts::DiamondBounds::alpha2() const
{
  return m_alpha2;
}

double
Acts::DiamondBounds::distanceToBoundary(const Acts::Vector2D& pos) const
{
  const int Np = 6;

  double y1 = 2. * m_valueStore.at(DiamondBounds::bv_halfY1);
  double y2 = 2. * m_valueStore.at(DiamondBounds::bv_halfY2);

  double X[6] = {-m_valueStore.at(DiamondBounds::bv_minHalfX),
                 -m_valueStore.at(DiamondBounds::bv_medHalfX),
                 -m_valueStore.at(DiamondBounds::bv_maxHalfX),
                 m_valueStore.at(DiamondBounds::bv_maxHalfX),
                 m_valueStore.at(DiamondBounds::bv_medHalfX),
                 m_valueStore.at(DiamondBounds::bv_minHalfX)};
  double Y[6] = {-y1, 0., y2, y2, 0., -y1};

  double dm = 1.e+20;
  double Ao = 0.;
  bool   in = true;

  for (int i = 0; i != Np; ++i) {
    int j          = i + 1;
    if (j == Np) j = 0;

    double x  = X[i] - pos[0];
    double y  = Y[i] - pos[1];
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
  if (in) return -sqrt(dm);
  return sqrt(dm);
}

// ostream operator overload
std::ostream&
Acts::DiamondBounds::dump(std::ostream& sl) const
{
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::DiamondBounds:  (minHlenghtX, medHlengthX, maxHlengthX, "
        "hlengthY1, hlengthY2 ) = ";
  sl << "(" << m_valueStore.at(DiamondBounds::bv_minHalfX) << ", "
     << m_valueStore.at(DiamondBounds::bv_medHalfX) << ", "
     << m_valueStore.at(DiamondBounds::bv_maxHalfX) << ", "
     << m_valueStore.at(DiamondBounds::bv_halfY1) << ", "
     << m_valueStore.at(DiamondBounds::bv_halfY2) << ")";
  sl << std::setprecision(-1);
  return sl;
}
