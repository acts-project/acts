// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// LineBounds.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Surfaces/LineBounds.hpp"
#include <iomanip>
#include <iostream>

Acts::LineBounds::LineBounds(double radius, double halez)
  : SurfaceBounds()
{
  if (radius != 0. || halez != 0.){
    m_valueStore.push_back(radius);
    m_valueStore.push_back(halez);
  }
}

Acts::LineBounds::~LineBounds()
{
}

Acts::LineBounds&
Acts::LineBounds::operator=(const Acts::LineBounds& libo)
{
  if (this != &libo)
    SurfaceBounds::operator=(libo);
  return *this;
}

double
Acts::LineBounds::minDistance(const Acts::Vector2D& lpos) const
{
  // per definition the min Distance of a correct local position on the line is r
  return lpos[Acts::eLOC_R];
}

// ostream operator overload
std::ostream&
Acts::LineBounds::dump(std::ostream& sl) const
{
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::LineBounds: (radius, halflengthInZ) = ";
  sl << "(" << r() << ", " << halflengthZ() << ")";
  sl << std::setprecision(-1);
  return sl;
}
