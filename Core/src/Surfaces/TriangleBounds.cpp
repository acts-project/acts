// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/TriangleBounds.hpp"

#include <iomanip>
#include <iostream>

Acts::TriangleBounds::TriangleBounds(const std::array<Vector2D, 3>& vertices)
    : m_vertices(vertices), m_boundingBox(0, 0) {
  double mx = 0;
  double my = 0;
  for (auto& v : vertices) {
    mx = std::max(mx, std::abs(v.x()));
    my = std::max(my, std::abs(v.y()));
  }
  m_boundingBox = RectangleBounds(mx, my);
}

Acts::TriangleBounds::~TriangleBounds() = default;

Acts::TriangleBounds* Acts::TriangleBounds::clone() const {
  return new TriangleBounds(*this);
}

Acts::SurfaceBounds::BoundsType Acts::TriangleBounds::type() const {
  return SurfaceBounds::Triangle;
}

std::vector<TDD_real_t> Acts::TriangleBounds::valueStore() const {
  std::vector<TDD_real_t> values(TriangleBounds::bv_length);
  values[TriangleBounds::bv_x1] = m_vertices[0].x();
  values[TriangleBounds::bv_y1] = m_vertices[0].y();
  values[TriangleBounds::bv_x2] = m_vertices[1].x();
  values[TriangleBounds::bv_y2] = m_vertices[1].y();
  values[TriangleBounds::bv_x3] = m_vertices[2].x();
  values[TriangleBounds::bv_y3] = m_vertices[2].y();
  return values;
}

bool Acts::TriangleBounds::inside(const Acts::Vector2D& lposition,
                                  const Acts::BoundaryCheck& bcheck) const {
  return bcheck.isInside(lposition, m_vertices);
}

double Acts::TriangleBounds::distanceToBoundary(
    const Acts::Vector2D& lposition) const {
  return BoundaryCheck(true).distance(lposition, m_vertices);
}

std::vector<Acts::Vector2D> Acts::TriangleBounds::vertices(
    unsigned int /*lseg*/) const {
  return {std::begin(m_vertices), std::end(m_vertices)};
}

const Acts::RectangleBounds& Acts::TriangleBounds::boundingBox() const {
  return m_boundingBox;
}

// ostream operator overload
std::ostream& Acts::TriangleBounds::toStream(std::ostream& sl) const {
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::TriangleBounds:  generating vertices (X, Y)";
  sl << "(" << m_vertices[0].x() << " , " << m_vertices[1].y() << ") " << '\n';
  sl << "(" << m_vertices[1].x() << " , " << m_vertices[1].y() << ") " << '\n';
  sl << "(" << m_vertices[2].x() << " , " << m_vertices[2].y() << ") ";
  sl << std::setprecision(-1);
  return sl;
}
