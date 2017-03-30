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
  : PlanarBounds(TrapezoidBounds::bv_length), m_boundingBox(0.0, 0.0)
{
  m_valueStore[TrapezoidBounds::bv_minHalfX] = std::abs(minhalex);
  m_valueStore[TrapezoidBounds::bv_maxHalfX] = std::abs(maxhalex);
  m_valueStore[TrapezoidBounds::bv_halfY]    = std::abs(haley);
  m_boundingBox = RectangleBounds(std::max(minhalex, maxhalex), haley);
}

Acts::TrapezoidBounds::~TrapezoidBounds()
{
}

Acts::TrapezoidBounds&
Acts::TrapezoidBounds::operator=(const TrapezoidBounds& trabo)
{
  if (this != &trabo) {
    PlanarBounds::operator=(trabo);
    m_boundingBox         = trabo.m_boundingBox;
  }
  return *this;
}

double
Acts::TrapezoidBounds::distanceToBoundary(const Acts::Vector2D& lpos) const
{
  return BoundaryCheck(true).distanceToBoundary(lpos, vertices());
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
