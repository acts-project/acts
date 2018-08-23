// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// ConeBounds.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Surfaces/ConeBounds.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

#include "Acts/Utilities/VariantData.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

Acts::ConeBounds::ConeBounds(double alpha,
                             bool   symm,
                             double halfphi,
                             double avphi)
  : ConeBounds(alpha,
               symm ? -std::numeric_limits<double>::infinity() : 0,
               std::numeric_limits<double>::infinity(),
               halfphi,
               avphi)
{
}

Acts::ConeBounds::ConeBounds(double alpha,
                             double zmin,
                             double zmax,
                             double halfphi,
                             double avphi)
  : m_alpha(alpha)
  , m_tanAlpha(std::tan(alpha))
  , m_zMin(zmin)
  , m_zMax(zmax)
  , m_avgPhi(detail::radian_sym(avphi))
  , m_halfPhi(std::abs(halfphi))
{
}

Acts::ConeBounds::ConeBounds(const variant_data& vardata)
{
  throw_assert(vardata.which() == 4, "Variant data must be map");
  const variant_map& data = boost::get<variant_map>(vardata);
  std::string        type = data.get<std::string>("type");
  throw_assert(type == "ConeBounds", "Type must be ConeBounds");

  const variant_map& payload = data.get<variant_map>("payload");

  m_alpha    = payload.get<double>("alpha");
  m_tanAlpha = std::tan(m_alpha);
  m_zMin     = payload.get<double>("zMin");
  m_zMax     = payload.get<double>("zMax");
  m_avgPhi   = payload.get<double>("avgPhi");
  m_halfPhi  = payload.get<double>("halfPhi");
}

Acts::ConeBounds::~ConeBounds() = default;

Acts::ConeBounds*
Acts::ConeBounds::clone() const
{
  return new ConeBounds(*this);
}

Acts::SurfaceBounds::BoundsType
Acts::ConeBounds::type() const
{
  return SurfaceBounds::Cone;
}

std::vector<TDD_real_t>
Acts::ConeBounds::valueStore() const
{
  std::vector<TDD_real_t> values(ConeBounds::bv_length);
  values[ConeBounds::bv_alpha]         = alpha();
  values[ConeBounds::bv_minZ]          = minZ();
  values[ConeBounds::bv_maxZ]          = maxZ();
  values[ConeBounds::bv_averagePhi]    = averagePhi();
  values[ConeBounds::bv_halfPhiSector] = halfPhiSector();
  return values;
}

/// Shift r-phi coordinate to be centered around the average phi.
Acts::Vector2D
Acts::ConeBounds::shifted(const Acts::Vector2D& lpos) const
{
  using Acts::detail::radian_sym;

  auto     x = r(lpos[eLOC_Z]);  // cone radius at the local position
  Vector2D shifted;
  shifted[eLOC_Z]    = lpos[eLOC_Z];
  shifted[eLOC_RPHI] = std::isnormal(x)
      ? (x * radian_sym((lpos[eLOC_RPHI] / x) - averagePhi()))
      : lpos[eLOC_RPHI];
  return shifted;
}

bool
Acts::ConeBounds::inside(const Acts::Vector2D&      lpos,
                         const Acts::BoundaryCheck& bcheck) const
{
  auto rphiHalf = r(lpos[eLOC_Z]) * halfPhiSector();
  return bcheck.isInside(
      shifted(lpos), Vector2D(-rphiHalf, minZ()), Vector2D(rphiHalf, maxZ()));
}

double
Acts::ConeBounds::distanceToBoundary(const Acts::Vector2D& lpos) const
{
  auto rphiHalf = r(lpos[eLOC_Z]) * halfPhiSector();
  return BoundaryCheck(true).distance(
      shifted(lpos), Vector2D(-rphiHalf, minZ()), Vector2D(rphiHalf, maxZ()));
}

std::ostream&
Acts::ConeBounds::dump(std::ostream& sl) const
{
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::ConeBounds: (tanAlpha, minZ, maxZ, averagePhi, halfPhiSector) "
        "= ";
  sl << "(" << m_tanAlpha << ", " << m_zMin << ", " << m_zMax << ", "
     << m_avgPhi << ", " << m_halfPhi << ")";
  sl << std::setprecision(-1);
  return sl;
}

Acts::variant_data
Acts::ConeBounds::toVariantData() const
{
  using namespace std::string_literals;

  variant_map payload;
  payload["alpha"]   = m_alpha;
  payload["zMin"]    = m_zMin;
  payload["zMax"]    = m_zMax;
  payload["avgPhi"]  = m_avgPhi;
  payload["halfPhi"] = m_halfPhi;

  return variant_map({{"type", "ConeBounds"s}, {"payload", payload}});
}
