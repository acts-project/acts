// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/DiamondBounds.hpp"
#include "Acts/Utilities/ThrowAssert.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>

Acts::DiamondBounds::DiamondBounds(double x1, double x2, double x3, double y1,
                                   double y2)
    : m_x1(std::abs(x1)),
      m_x2(std::abs(x2)),
      m_x3(std::abs(x3)),
      m_y1(std::abs(y1)),
      m_y2(std::abs(y2)),
      m_boundingBox(std::max(std::max(x1, x2), x3), std::max(y1, y2)) {
  throw_assert((x1 <= x2), "Hexagon must be a convex polygon");
  throw_assert((x3 <= x2), "Hexagon must be a convex polygon");
}

Acts::DiamondBounds* Acts::DiamondBounds::clone() const {
  return new DiamondBounds(*this);
}

Acts::SurfaceBounds::BoundsType Acts::DiamondBounds::type() const {
  return SurfaceBounds::Diamond;
}

std::vector<TDD_real_t> Acts::DiamondBounds::valueStore() const {
  std::vector<TDD_real_t> values(DiamondBounds::bv_length);
  values[DiamondBounds::bv_x1] = x1();
  values[DiamondBounds::bv_x2] = x2();
  values[DiamondBounds::bv_x3] = x3();
  values[DiamondBounds::bv_y1] = y1();
  values[DiamondBounds::bv_y2] = y2();
  return values;
}

bool Acts::DiamondBounds::inside(const Acts::Vector2D& lposition,
                                 const Acts::BoundaryCheck& bcheck) const {
  return bcheck.isInside(lposition, vertices());
}

double Acts::DiamondBounds::distanceToBoundary(
    const Acts::Vector2D& lposition) const {
  return BoundaryCheck(true).distance(lposition, vertices());
}

std::vector<Acts::Vector2D> Acts::DiamondBounds::vertices(
    unsigned int /*lseg*/) const {
  // Vertices starting at lower left (min rel. phi)
  // counter-clockwise
  return {{-x1(), -y1()}, {x1(), -y1()}, {x2(), 0.},
          {x3(), y2()},   {-x3(), y2()}, {-x2(), 0.}};
}

const Acts::RectangleBounds& Acts::DiamondBounds::boundingBox() const {
  return m_boundingBox;
}

// ostream operator overload
std::ostream& Acts::DiamondBounds::toStream(std::ostream& sl) const {
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::DiamondBounds:  (x1, x2, x3, "
        "y1, y2 ) = ";
  sl << "(" << x1() << ", " << x2() << ", " << x3() << ", " << y1() << ", "
     << y2() << ")";
  sl << std::setprecision(-1);
  return sl;
}
