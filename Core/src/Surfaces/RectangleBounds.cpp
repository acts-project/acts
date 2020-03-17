// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/ThrowAssert.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>

Acts::RectangleBounds::RectangleBounds(double halex, double haley)
    : m_min(-halex, -haley), m_max(halex, haley) {}

Acts::RectangleBounds::RectangleBounds(const Vector2D& vmin,
                                       const Vector2D& vmax)
    : m_min(vmin), m_max(vmax) {}

Acts::RectangleBounds* Acts::RectangleBounds::clone() const {
  return new RectangleBounds(*this);
}

std::vector<TDD_real_t> Acts::RectangleBounds::valueStore() const {
  std::vector<TDD_real_t> values(RectangleBounds::bv_length);
  values = {m_min.x(), m_min.y(), m_max.x(), m_max.y()};
  return values;
}

bool Acts::RectangleBounds::inside(const Acts::Vector2D& lposition,
                                   const Acts::BoundaryCheck& bcheck) const {
  return bcheck.isInside(lposition, m_min, m_max);
}

double Acts::RectangleBounds::distanceToBoundary(
    const Acts::Vector2D& lposition) const {
  return BoundaryCheck(true).distance(lposition, m_min, m_max);
}

std::vector<Acts::Vector2D> Acts::RectangleBounds::vertices(
    unsigned int /*lseg*/) const {
  // counter-clockwise starting from bottom-left corner
  return {m_min, {m_max.x(), m_min.y()}, m_max, {m_min.x(), m_max.y()}};
}

const Acts::RectangleBounds& Acts::RectangleBounds::boundingBox() const {
  return (*this);
}

// ostream operator overload
std::ostream& Acts::RectangleBounds::toStream(std::ostream& sl) const {
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::RectangleBounds:  (hlX, hlY) = "
     << "(" << halflengthX() << ", " << halflengthY() << ")";
  sl << "\n(lower left, upper right):\n";
  sl << m_min.transpose() << "\n" << m_max.transpose();
  sl << std::setprecision(-1);
  return sl;
}
