// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// RadialBounds.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Surfaces/RadialBounds.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>

#include "Acts/Utilities/VariantData.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

Acts::RadialBounds::RadialBounds(double minrad, double maxrad, double hphisec)
  : RadialBounds(minrad, maxrad, 0, hphisec)
{
}

Acts::RadialBounds::RadialBounds(double minrad,
                                 double maxrad,
                                 double avephi,
                                 double hphisec)
  : m_rMin(std::min(minrad, maxrad))
  , m_rMax(std::max(minrad, maxrad))
  , m_avgPhi(detail::radian_sym(avephi))
  , m_halfPhi(std::abs(hphisec))
{
}

Acts::RadialBounds::RadialBounds(const variant_data& vardata)
{

  throw_assert(vardata.which() == 4, "Variant data must be map");
  const variant_map& data = boost::get<variant_map>(vardata);
  std::string        type = data.get<std::string>("type");
  throw_assert(type == "RadialBounds", "Type must be RadialBounds");

  const variant_map& payload = data.get<variant_map>("payload");

  m_rMin    = payload.get<double>("rMin");
  m_rMax    = payload.get<double>("rMax");
  m_avgPhi  = payload.get<double>("avgPhi");
  m_halfPhi = payload.get<double>("halfPhi");
}

Acts::RadialBounds::~RadialBounds() = default;

Acts::RadialBounds*
Acts::RadialBounds::clone() const
{
  return new RadialBounds(*this);
}

Acts::SurfaceBounds::BoundsType
Acts::RadialBounds::type() const
{
  return SurfaceBounds::Disc;
}

std::vector<TDD_real_t>
Acts::RadialBounds::valueStore() const
{
  std::vector<TDD_real_t> values(RadialBounds::bv_length);
  values[RadialBounds::bv_rMin]          = rMin();
  values[RadialBounds::bv_rMax]          = rMax();
  values[RadialBounds::bv_averagePhi]    = averagePhi();
  values[RadialBounds::bv_halfPhiSector] = halfPhiSector();
  return values;
}

Acts::Vector2D
Acts::RadialBounds::shifted(const Acts::Vector2D& lpos) const
{
  Vector2D tmp;
  tmp[eLOC_R]   = lpos[eLOC_R];
  tmp[eLOC_PHI] = detail::radian_sym(lpos[eLOC_PHI] - averagePhi());
  return tmp;
}

bool
Acts::RadialBounds::inside(const Acts::Vector2D&      lpos,
                           const Acts::BoundaryCheck& bcheck) const
{
  return bcheck.isInside(shifted(lpos),
                         Vector2D(rMin(), -halfPhiSector()),
                         Vector2D(rMax(), halfPhiSector()));
}

double
Acts::RadialBounds::distanceToBoundary(const Acts::Vector2D& lpos) const
{
  return BoundaryCheck(true).distance(shifted(lpos),
                                      Vector2D(rMin(), -halfPhiSector()),
                                      Vector2D(rMax(), halfPhiSector()));
}

// ostream operator overload
std::ostream&
Acts::RadialBounds::dump(std::ostream& sl) const
{
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::RadialBounds:  (innerRadius, outerRadius, hPhiSector) = ";
  sl << "(" << rMin() << ", " << rMax() << ", " << averagePhi() << ", "
     << halfPhiSector() << ")";
  sl << std::setprecision(-1);
  return sl;
}

Acts::variant_data
Acts::RadialBounds::toVariantData() const
{
  using namespace std::string_literals;

  variant_map payload;
  payload["rMin"]    = m_rMin;
  payload["rMax"]    = m_rMax;
  payload["avgPhi"]  = m_avgPhi;
  payload["halfPhi"] = m_halfPhi;

  return variant_map({{"type", "RadialBounds"s}, {"payload", payload}});
}
