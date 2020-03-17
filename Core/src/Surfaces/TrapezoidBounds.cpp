// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/TrapezoidBounds.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>

Acts::TrapezoidBounds::TrapezoidBounds(double minhalex, double maxhalex,
                                       double haley)
    : m_minHalfX(std::abs(minhalex)),
      m_maxHalfX(std::abs(maxhalex)),
      m_halfY(std::abs(haley)),
      m_boundingBox(std::max(minhalex, maxhalex), haley) {}

Acts::TrapezoidBounds::~TrapezoidBounds() = default;

Acts::TrapezoidBounds* Acts::TrapezoidBounds::clone() const {
  return new TrapezoidBounds(*this);
}

Acts::SurfaceBounds::BoundsType Acts::TrapezoidBounds::type() const {
  return SurfaceBounds::Trapezoid;
}

std::vector<TDD_real_t> Acts::TrapezoidBounds::valueStore() const {
  std::vector<TDD_real_t> values(TrapezoidBounds::bv_length);
  values[TrapezoidBounds::bv_minHalfX] = minHalflengthX();
  values[TrapezoidBounds::bv_maxHalfX] = maxHalflengthX();
  values[TrapezoidBounds::bv_halfY] = halflengthY();
  return values;
}

bool Acts::TrapezoidBounds::inside(const Acts::Vector2D& lposition,
                                   const Acts::BoundaryCheck& bcheck) const {
  return bcheck.isInside(lposition, vertices());
}

double Acts::TrapezoidBounds::distanceToBoundary(
    const Acts::Vector2D& lposition) const {
  return BoundaryCheck(true).distance(lposition, vertices());
}

std::vector<Acts::Vector2D> Acts::TrapezoidBounds::vertices(
    unsigned int /*lseg*/) const {
  // counter-clockwise from bottom-right corner
  return {{minHalflengthX(), -halflengthY()},
          {maxHalflengthX(), halflengthY()},
          {-maxHalflengthX(), halflengthY()},
          {-minHalflengthX(), -halflengthY()}};
}

const Acts::RectangleBounds& Acts::TrapezoidBounds::boundingBox() const {
  return m_boundingBox;
}

std::ostream& Acts::TrapezoidBounds::toStream(std::ostream& sl) const {
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::TrapezoidBounds:  (minHlengthX, maxHlengthX, hlengthY) = "
     << "(" << minHalflengthX() << ", " << maxHalflengthX() << ", "
     << halflengthY() << ")";
  sl << std::setprecision(-1);
  return sl;
}
