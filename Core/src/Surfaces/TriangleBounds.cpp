// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// TriangleBounds.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Surfaces/TriangleBounds.hpp"
#include "Acts/Utilities/VariantData.hpp"

#include <iomanip>
#include <iostream>

Acts::TriangleBounds::TriangleBounds(const std::array<Vector2D, 3>& vertices)
  : m_vertices(vertices), m_boundingBox(0, 0)
{
  double mx = 0;
  double my = 0;
  for (auto& v : vertices) {
    mx = std::max(mx, std::abs(v.x()));
    my = std::max(my, std::abs(v.y()));
  }
  m_boundingBox = RectangleBounds(mx, my);
}

Acts::TriangleBounds::TriangleBounds(const variant_data& vardata)
  : m_vertices(), m_boundingBox(0, 0)
{
  throw_assert(vardata.which() == 4, "Variant data must be map");
  const variant_map& data = boost::get<variant_map>(vardata);
  std::string        type = data.get<std::string>("type");
  throw_assert(type == "TriangleBounds", "Type must be TriangleBounds");

  const variant_map&    payload  = data.get<variant_map>("payload");
  const variant_vector& vertices = payload.get<variant_vector>("vertices");
  throw_assert(vertices.size() == 3,
               "Vertices for triangle must be exactly 3.");

  double mx = 0, my = 0;
  for (size_t i = 0; i < 3; i++) {
    Vector2D vtx     = from_variant<Vector2D>(vertices.at(i));
    mx               = std::max(mx, std::abs(vtx.x()));
    my               = std::max(my, std::abs(vtx.y()));
    m_vertices.at(i) = vtx;
  }

  m_boundingBox = RectangleBounds(mx, my);
}

Acts::TriangleBounds::~TriangleBounds() = default;

Acts::TriangleBounds*
Acts::TriangleBounds::clone() const
{
  return new TriangleBounds(*this);
}

Acts::SurfaceBounds::BoundsType
Acts::TriangleBounds::type() const
{
  return SurfaceBounds::Triangle;
}

std::vector<TDD_real_t>
Acts::TriangleBounds::valueStore() const
{
  std::vector<TDD_real_t> values(TriangleBounds::bv_length);
  values[TriangleBounds::bv_x1] = m_vertices[0].x();
  values[TriangleBounds::bv_y1] = m_vertices[0].y();
  values[TriangleBounds::bv_x2] = m_vertices[1].x();
  values[TriangleBounds::bv_y2] = m_vertices[1].y();
  values[TriangleBounds::bv_x3] = m_vertices[2].x();
  values[TriangleBounds::bv_y3] = m_vertices[2].y();
  return values;
}

bool
Acts::TriangleBounds::inside(const Acts::Vector2D&      lpos,
                             const Acts::BoundaryCheck& bcheck) const
{
  return bcheck.isInside(lpos, m_vertices);
}

double
Acts::TriangleBounds::distanceToBoundary(const Acts::Vector2D& lpos) const
{
  return BoundaryCheck(true).distance(lpos, m_vertices);
}

std::vector<Acts::Vector2D>
Acts::TriangleBounds::vertices() const
{
  return {std::begin(m_vertices), std::end(m_vertices)};
}

const Acts::RectangleBounds&
Acts::TriangleBounds::boundingBox() const
{
  return m_boundingBox;
}

// ostream operator overload
std::ostream&
Acts::TriangleBounds::dump(std::ostream& sl) const
{
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::TriangleBounds:  generating vertices (X, Y)";
  sl << "(" << m_vertices[0].x() << " , " << m_vertices[1].y() << ") " << '\n';
  sl << "(" << m_vertices[1].x() << " , " << m_vertices[1].y() << ") " << '\n';
  sl << "(" << m_vertices[2].x() << " , " << m_vertices[2].y() << ") ";
  sl << std::setprecision(-1);
  return sl;
}

Acts::variant_data
Acts::TriangleBounds::toVariantData() const
{
  using namespace std::string_literals;

  variant_map payload;

  variant_vector vertices;
  for (const auto& vtx : m_vertices) {
    vertices.push_back(to_variant(vtx));
  }

  variant_map data;
  data["type"]    = "TriangleBounds"s;
  data["payload"] = variant_map({{"vertices", vertices}});

  return data;
}
