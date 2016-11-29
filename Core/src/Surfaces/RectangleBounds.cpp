// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// RectangleBounds.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Surfaces/RectangleBounds.hpp"
#include <iomanip>
#include <iostream>
#include <cmath>

Acts::RectangleBounds::RectangleBounds(double halex, double haley)
  : Acts::PlanarBounds(RectangleBounds::bv_length)
{
  m_valueStore.at(RectangleBounds::bv_halfX) = halex;
  m_valueStore.at(RectangleBounds::bv_halfY) = haley;
}

Acts::RectangleBounds::~RectangleBounds()
{
}

Acts::RectangleBounds&
Acts::RectangleBounds::operator=(const RectangleBounds& recbo)
{
  if (this != &recbo) PlanarBounds::operator=(recbo);
  return *this;
}

double
Acts::RectangleBounds::distanceToBoundary(const Acts::Vector2D& lpos) const
{
  double dx = std::abs(lpos[0]) - m_valueStore.at(RectangleBounds::bv_halfX);
  double dy = std::abs(lpos[1]) - m_valueStore.at(RectangleBounds::bv_halfY);

  if (dx <= 0. || dy <= 0.) {
    if (dx > dy)
      return dx;
    else
      return dy;
  }
  return std::sqrt(dx * dx + dy * dy);
}

// ostream operator overload
std::ostream&
Acts::RectangleBounds::dump(std::ostream& sl) const
{
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::RectangleBounds:  (halflenghtX, halflengthY) = "
     << "(" << m_valueStore.at(RectangleBounds::bv_halfX) << ", "
     << m_valueStore.at(RectangleBounds::bv_halfY) << ")";
  sl << std::setprecision(-1);
  return sl;
}
