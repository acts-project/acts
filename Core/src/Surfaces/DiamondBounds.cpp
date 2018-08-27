// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// DiamondBounds.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Surfaces/DiamondBounds.hpp"
#include "Acts/Utilities/ThrowAssert.hpp"
#include "Acts/Utilities/VariantData.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>

Acts::DiamondBounds::DiamondBounds(double minhalex,
                                   double medhalex,
                                   double maxhalex,
                                   double haley1,
                                   double haley2)
  : m_minHalfX(std::abs(minhalex))
  , m_medHalfX(std::abs(medhalex))
  , m_maxHalfX(std::abs(maxhalex))
  , m_minY(std::abs(haley1))
  , m_maxY(std::abs(haley2))
  , m_boundingBox(std::max(std::max(minhalex, medhalex), maxhalex),
                  std::max(haley1, haley2))
{
  throw_assert((minhalex <= medhalex), "Hexagon must be a convex polygon");
  throw_assert((maxhalex <= medhalex), "Hexagon must be a convex polygon");
}

Acts::DiamondBounds::DiamondBounds(const variant_data& vardata)
  : m_boundingBox(0, 0)
{
  throw_assert(vardata.which() == 4, "Variant data must be map");
  const variant_map& data = boost::get<variant_map>(vardata);
  std::string        type = data.get<std::string>("type");
  throw_assert(type == "DiamondBounds", "Type must be DiamondBounds");

  const variant_map& payload = data.get<variant_map>("payload");

  m_minHalfX = payload.get<double>("minHalfX");
  m_medHalfX = payload.get<double>("medHalfX");
  m_maxHalfX = payload.get<double>("maxHalfX");
  m_minY     = payload.get<double>("minY");
  m_maxY     = payload.get<double>("maxY");

  m_boundingBox
      = RectangleBounds(std::max(std::max(m_minHalfX, m_medHalfX), m_maxHalfX),
                        std::max(m_minY, m_maxY));
}

Acts::DiamondBounds::~DiamondBounds() = default;

Acts::DiamondBounds*
Acts::DiamondBounds::clone() const
{
  return new DiamondBounds(*this);
}

Acts::SurfaceBounds::BoundsType
Acts::DiamondBounds::type() const
{
  return SurfaceBounds::Diamond;
}

std::vector<TDD_real_t>
Acts::DiamondBounds::valueStore() const
{
  std::vector<TDD_real_t> values(DiamondBounds::bv_length);
  values[DiamondBounds::bv_minHalfX] = minHalflengthX();
  values[DiamondBounds::bv_medHalfX] = medHalflengthX();
  values[DiamondBounds::bv_maxHalfX] = maxHalflengthX();
  values[DiamondBounds::bv_halfY1]   = halflengthY1();
  values[DiamondBounds::bv_halfY2]   = halflengthY2();
  return values;
}

bool
Acts::DiamondBounds::inside(const Acts::Vector2D&      lpos,
                            const Acts::BoundaryCheck& bcheck) const
{
  return bcheck.isInside(lpos, vertices());
}

double
Acts::DiamondBounds::distanceToBoundary(const Acts::Vector2D& lpos) const
{
  return BoundaryCheck(true).distance(lpos, vertices());
}

std::vector<Acts::Vector2D>
Acts::DiamondBounds::vertices() const
{
  // vertices starting from lower right in clock-wise order
  return {{minHalflengthX(), -halflengthY1()},
          {medHalflengthX(), 0},
          {maxHalflengthX(), halflengthY2()},
          {-maxHalflengthX(), halflengthY2()},
          {-medHalflengthX(), 0},
          {-minHalflengthX(), -halflengthY1()}};
}

const Acts::RectangleBounds&
Acts::DiamondBounds::boundingBox() const
{
  return m_boundingBox;
}

// ostream operator overload
std::ostream&
Acts::DiamondBounds::dump(std::ostream& sl) const
{
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::DiamondBounds:  (minHlengthX, medHlengthX, maxHlengthX, "
        "hlengthY1, hlengthY2 ) = ";
  sl << "(" << minHalflengthX() << ", " << medHalflengthX() << ", "
     << maxHalflengthX() << ", " << halflengthY1() << ", " << halflengthY2()
     << ")";
  sl << std::setprecision(-1);
  return sl;
}

Acts::variant_data
Acts::DiamondBounds::toVariantData() const
{
  using namespace std::string_literals;

  variant_map payload;
  payload["minHalfX"] = m_minHalfX;
  payload["medHalfX"] = m_medHalfX;
  payload["maxHalfX"] = m_maxHalfX;
  payload["minY"]     = m_minY;
  payload["maxY"]     = m_maxY;

  variant_map data;
  data["type"]    = "DiamondBounds"s;
  data["payload"] = payload;

  return data;
}
