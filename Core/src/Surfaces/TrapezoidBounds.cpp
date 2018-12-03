// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// TrapezoidBounds.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Surfaces/TrapezoidBounds.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include "Acts/Utilities/VariantData.hpp"

Acts::TrapezoidBounds::TrapezoidBounds(double minhalex,
                                       double maxhalex,
                                       double haley)
  : m_minHalfX(std::abs(minhalex))
  , m_maxHalfX(std::abs(maxhalex))
  , m_halfY(std::abs(haley))
  , m_boundingBox(std::max(minhalex, maxhalex), haley)
{
}

Acts::TrapezoidBounds::TrapezoidBounds(const variant_data& vardata)
  : m_boundingBox(0, 0)
{
  throw_assert(vardata.which() == 4, "Variant data must be map");
  const variant_map& data = boost::get<variant_map>(vardata);
  std::string        type = data.get<std::string>("type");
  throw_assert(type == "TrapezoidBounds", "Type must be TrapezoidBounds");

  const variant_map& payload = data.get<variant_map>("payload");

  m_minHalfX = payload.get<double>("minHalfX");
  m_maxHalfX = payload.get<double>("maxHalfX");
  m_halfY    = payload.get<double>("halfY");

  m_boundingBox = RectangleBounds(m_maxHalfX, m_halfY);
}

Acts::TrapezoidBounds::~TrapezoidBounds() = default;

Acts::TrapezoidBounds*
Acts::TrapezoidBounds::clone() const
{
  return new TrapezoidBounds(*this);
}

Acts::SurfaceBounds::BoundsType
Acts::TrapezoidBounds::type() const
{
  return SurfaceBounds::Trapezoid;
}

std::vector<TDD_real_t>
Acts::TrapezoidBounds::valueStore() const
{
  std::vector<TDD_real_t> values(TrapezoidBounds::bv_length);
  values[TrapezoidBounds::bv_minHalfX] = minHalflengthX();
  values[TrapezoidBounds::bv_maxHalfX] = maxHalflengthX();
  values[TrapezoidBounds::bv_halfY]    = halflengthY();
  return values;
}

bool
Acts::TrapezoidBounds::inside(const Acts::Vector2D&      lpos,
                              const Acts::BoundaryCheck& bcheck) const
{
  return bcheck.isInside(lpos, vertices());
}

double
Acts::TrapezoidBounds::distanceToBoundary(const Acts::Vector2D& lpos) const
{
  return BoundaryCheck(true).distance(lpos, vertices());
}

std::vector<Acts::Vector2D>
Acts::TrapezoidBounds::vertices() const
{
  // counter-clockwise from bottom-right corner
  return {{minHalflengthX(), -halflengthY()},
          {maxHalflengthX(), halflengthY()},
          {-maxHalflengthX(), halflengthY()},
          {-minHalflengthX(), -halflengthY()}};
}

const Acts::RectangleBounds&
Acts::TrapezoidBounds::boundingBox() const
{
  return m_boundingBox;
}

std::ostream&
Acts::TrapezoidBounds::dump(std::ostream& sl) const
{
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::TrapezoidBounds:  (minHlengthX, maxHlengthX, hlengthY) = "
     << "(" << minHalflengthX() << ", " << maxHalflengthX() << ", "
     << halflengthY() << ")";
  sl << std::setprecision(-1);
  return sl;
}

Acts::variant_data
Acts::TrapezoidBounds::toVariantData() const
{
  using namespace std::string_literals;

  variant_map payload;
  payload["minHalfX"] = m_minHalfX;
  payload["maxHalfX"] = m_maxHalfX;
  payload["halfY"]    = m_halfY;

  variant_map data;
  data["type"]    = "TrapezoidBounds"s;
  data["payload"] = payload;

  return data;
}
