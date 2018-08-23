// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// LineBounds.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Surfaces/LineBounds.hpp"
#include "Acts/Utilities/VariantData.hpp"

#include <iomanip>
#include <iostream>

Acts::LineBounds::LineBounds(double radius, double halez)
  : m_radius(std::abs(radius)), m_halfZ(std::abs(halez))
{
}

Acts::LineBounds::LineBounds(const variant_data& vardata)
{

  throw_assert(vardata.which() == 4, "Variant data must be map");
  variant_map data = boost::get<variant_map>(vardata);
  throw_assert(data.count("type"), "Variant data must have type.");
  // std::string type = boost::get<std::string>(data["type"]);
  std::string type = data.get<std::string>("type");
  throw_assert(type == "LineBounds", "Variant data type must be LineBounds");

  variant_map payload = data.get<variant_map>("payload");

  m_radius = payload.get<double>("radius");
  m_halfZ  = payload.get<double>("halfZ");
}

Acts::LineBounds::~LineBounds() = default;

Acts::LineBounds*
Acts::LineBounds::clone() const
{
  return new LineBounds(*this);
}

Acts::SurfaceBounds::BoundsType
Acts::LineBounds::type() const
{
  return SurfaceBounds::Line;
}

std::vector<TDD_real_t>
Acts::LineBounds::valueStore() const
{
  std::vector<TDD_real_t> values(LineBounds::bv_length);
  values[LineBounds::bv_radius] = r();
  values[LineBounds::bv_halfZ]  = halflengthZ();
  return values;
}

bool
Acts::LineBounds::inside(const Acts::Vector2D&      lpos,
                         const Acts::BoundaryCheck& bcheck) const
{
  return bcheck.isInside(
      lpos, Vector2D(0, -halflengthZ()), Vector2D(r(), halflengthZ()));
}

double
Acts::LineBounds::distanceToBoundary(const Acts::Vector2D& lpos) const
{
  // per definition the min Distance of a correct local position is r
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

Acts::variant_data
Acts::LineBounds::toVariantData() const
{
  using namespace std::string_literals;
  variant_map payload;

  payload["radius"] = m_radius;
  payload["halfZ"]  = m_halfZ;

  variant_map data;
  data["type"]    = "LineBounds"s;
  data["payload"] = payload;
  return data;
}
