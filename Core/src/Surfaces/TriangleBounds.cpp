// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// TriangleBounds.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Surfaces/TriangleBounds.hpp"
#include <iomanip>
#include <iostream>

Acts::TriangleBounds::TriangleBounds(const std::vector<Vector2D>& vertices)
  : PlanarBounds(), m_boundingBox(0., 0.)
{
  m_valueStore.reserve(6);
  double mx, my = 0.;
  for (auto& v : vertices) {
    m_valueStore.push_back(v.x());
    m_valueStore.push_back(v.y());
    mx = (fabs(v.x()) > mx) ? fabs(v.x()) : mx;
    my = (fabs(v.y()) > my) ? fabs(v.y()) : my;
  }
  m_boundingBox = RectangleBounds(mx, my);
}

Acts::TriangleBounds::~TriangleBounds()
{
}

Acts::TriangleBounds&
Acts::TriangleBounds::operator=(const TriangleBounds& tribo)
{
  if (this != &tribo) {
    PlanarBounds::operator=(tribo);
    m_boundingBox         = tribo.m_boundingBox;
  }
  return *this;
}

double
Acts::TriangleBounds::distanceToBoundary(const Acts::Vector2D& lpos) const
{
  return BoundaryCheck(true).distanceToBoundary(lpos, vertices());
}

// ostream operator overload
std::ostream&
Acts::TriangleBounds::dump(std::ostream& sl) const
{
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::TriangleBounds:  generating vertices (X, Y)";
  sl << "(" << m_valueStore.at(TriangleBounds::bv_x1) << " , "
     << m_valueStore.at(TriangleBounds::bv_y1) << ") " << '\n';
  sl << "(" << m_valueStore.at(TriangleBounds::bv_x2) << " , "
     << m_valueStore.at(TriangleBounds::bv_y2) << ") " << '\n';
  sl << "(" << m_valueStore.at(TriangleBounds::bv_x3) << " , "
     << m_valueStore.at(TriangleBounds::bv_y3) << ") ";
  sl << std::setprecision(-1);
  return sl;
}
