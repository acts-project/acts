// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// RectangleBounds.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/ThrowAssert.hpp"
#include "Acts/Utilities/VariantData.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>

Acts::RectangleBounds::RectangleBounds(double halex, double haley)
  : m_halfX(std::abs(halex)), m_halfY(std::abs(haley))
{
}

Acts::RectangleBounds::RectangleBounds(const variant_data& vardata)
{
  throw_assert(vardata.which() == 4, "Variant data must be map");
  const variant_map& data = boost::get<variant_map>(vardata);
  std::string        type = data.get<std::string>("type");
  throw_assert(type == "RectangleBounds", "Type must be RectangleBounds");

  const variant_map& payload = data.get<variant_map>("payload");

  m_halfX = payload.get<double>("halflengthX");
  m_halfY = payload.get<double>("halflengthY");
}

Acts::RectangleBounds::~RectangleBounds() = default;

Acts::RectangleBounds*
Acts::RectangleBounds::clone() const
{
  return new RectangleBounds(*this);
}

std::vector<TDD_real_t>
Acts::RectangleBounds::valueStore() const
{
  std::vector<TDD_real_t> values(RectangleBounds::bv_length);
  values[RectangleBounds::bv_halfX] = halflengthX();
  values[RectangleBounds::bv_halfY] = halflengthY();
  return values;
}

bool
Acts::RectangleBounds::inside(const Acts::Vector2D&      lpos,
                              const Acts::BoundaryCheck& bcheck) const
{
  return bcheck.isInside(lpos,
                         Vector2D(-halflengthX(), -halflengthY()),
                         Vector2D(halflengthX(), halflengthY()));
}

double
Acts::RectangleBounds::distanceToBoundary(const Acts::Vector2D& lpos) const
{
  return BoundaryCheck(true).distance(lpos,
                                      Vector2D(-halflengthX(), -halflengthY()),
                                      Vector2D(halflengthX(), halflengthY()));
}

std::vector<Acts::Vector2D>
Acts::RectangleBounds::vertices() const
{
  // counter-clockwise starting from bottom-right corner
  return {{halflengthX(), -halflengthY()},
          {halflengthX(), halflengthY()},
          {-halflengthX(), halflengthY()},
          {-halflengthX(), -halflengthY()}};
}

const Acts::RectangleBounds&
Acts::RectangleBounds::boundingBox() const
{
  return (*this);
}

// ostream operator overload
std::ostream&
Acts::RectangleBounds::dump(std::ostream& sl) const
{
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::RectangleBounds:  (halflengthX, halflengthY) = "
     << "(" << halflengthX() << ", " << halflengthY() << ")";
  sl << std::setprecision(-1);
  return sl;
}

Acts::variant_data
Acts::RectangleBounds::toVariantData() const
{
  using namespace std::string_literals;

  variant_map payload;
  payload["halflengthX"] = m_halfX;
  payload["halflengthY"] = m_halfY;

  variant_map data({{"type", "RectangleBounds"s}, {"payload", payload}});

  return data;
}
