// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/ConeBounds.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

Acts::ConeBounds::ConeBounds(double alpha, bool symm, double halfphi,
                             double avphi)
    : m_values({alpha, symm ? -std::numeric_limits<double>::infinity() : 0,
                std::numeric_limits<double>::infinity(), std::abs(halfphi),
                detail::radian_sym(avphi)}),
      m_tanAlpha(std::tan(alpha)) {}

Acts::ConeBounds::ConeBounds(double alpha, double minz, double maxz,
                             double halfphi, double avphi)
    : m_values(
          {alpha, minz, maxz, std::abs(halfphi), detail::radian_sym(avphi)}),
      m_tanAlpha(std::tan(alpha)) {}

Acts::ConeBounds::ConeBounds(const std::array<double, eSize>& values)
    : m_values(values), m_tanAlpha(std::tan(values[eAlpha])) {}

Acts::ConeBounds* Acts::ConeBounds::clone() const {
  return new ConeBounds(*this);
}

Acts::SurfaceBounds::BoundsType Acts::ConeBounds::type() const {
  return SurfaceBounds::eCone;
}

/// Shift r-phi coordinate to be centered around the average phi.
Acts::Vector2D Acts::ConeBounds::shifted(
    const Acts::Vector2D& lposition) const {
  using Acts::detail::radian_sym;

  auto x = r(lposition[eLOC_Z]);  // cone radius at the local position
  Vector2D shifted;
  shifted[eLOC_Z] = lposition[eLOC_Z];
  shifted[eLOC_RPHI] =
      std::isnormal(x)
          ? (x * radian_sym((lposition[eLOC_RPHI] / x) - get(eAveragePhi)))
          : lposition[eLOC_RPHI];
  return shifted;
}

bool Acts::ConeBounds::inside(const Acts::Vector2D& lposition,
                              const Acts::BoundaryCheck& bcheck) const {
  auto rphiHalf = r(lposition[eLOC_Z]) * get(eHalfPhiSector);
  return bcheck.isInside(shifted(lposition), Vector2D(-rphiHalf, get(eMinZ)),
                         Vector2D(rphiHalf, get(eMaxZ)));
}

double Acts::ConeBounds::distanceToBoundary(
    const Acts::Vector2D& lposition) const {
  auto rphiHalf = r(lposition[eLOC_Z]) * get(eHalfPhiSector);
  return BoundaryCheck(true).distance(shifted(lposition),
                                      Vector2D(-rphiHalf, get(eMinZ)),
                                      Vector2D(rphiHalf, get(eMaxZ)));
}

std::ostream& Acts::ConeBounds::toStream(std::ostream& sl) const {
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::ConeBounds: (tanAlpha, minZ, maxZ, halfPhiSector, averagePhi) "
        "= ";
  sl << "(" << m_tanAlpha << ", " << get(eMinZ) << ", " << get(eMaxZ) << ", "
     << get(eHalfPhiSector) << ", " << get(eAveragePhi) << ")";
  sl << std::setprecision(-1);
  return sl;
}
