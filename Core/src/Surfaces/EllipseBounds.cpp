// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/EllipseBounds.hpp"
#include "Acts/Surfaces/detail/VerticesHelper.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>

#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;

Acts::EllipseBounds::EllipseBounds(double minRadius0, double minRadius1,
                                   double maxRadius0, double maxRadius1,
                                   double averagePhi, double halfPhi)
    : m_rMinX(std::abs(minRadius0)),
      m_rMinY(std::abs(minRadius1)),
      m_rMaxX(std::abs(maxRadius0)),
      m_rMaxY(std::abs(maxRadius1)),
      m_avgPhi(detail::radian_sym(averagePhi)),
      m_halfPhi(std::abs(halfPhi)),
      m_boundingBox(std::max(minRadius0, maxRadius0),
                    std::max(minRadius1, maxRadius1)) {}

Acts::EllipseBounds* Acts::EllipseBounds::clone() const {
  return new EllipseBounds(*this);
}

Acts::SurfaceBounds::BoundsType Acts::EllipseBounds::type() const {
  return SurfaceBounds::Ellipse;
}

std::vector<TDD_real_t> Acts::EllipseBounds::valueStore() const {
  std::vector<TDD_real_t> values(EllipseBounds::bv_length);
  values[EllipseBounds::bv_rMinX] = m_rMinX;
  values[EllipseBounds::bv_rMinY] = m_rMinY;
  values[EllipseBounds::bv_rMaxX] = m_rMaxX;
  values[EllipseBounds::bv_rMaxY] = m_rMaxY;
  values[EllipseBounds::bv_averagePhi] = m_avgPhi;
  values[EllipseBounds::bv_halfPhiSector] = m_halfPhi;
  return values;
}

static inline double square(double x) {
  return x * x;
}

/// @warning This **only** works for tolerance-based checks
bool Acts::EllipseBounds::inside(const Acts::Vector2D& lposition,
                                 const Acts::BoundaryCheck& bcheck) const {
  double tol0 = bcheck.m_tolerance[0];
  double tol1 = bcheck.m_tolerance[1];
  double phi = detail::radian_sym(VectorHelpers::phi(lposition) - averagePhi());
  double phiHalf = halfPhiSector() + tol1;

  bool insidePhi = (-phiHalf <= phi) && (phi < phiHalf);
  bool insideInner = (rMinX() <= tol0) || (rMinY() <= tol0) ||
                     (1 < (square(lposition[Acts::eLOC_X] / (rMinX() - tol0)) +
                           square(lposition[Acts::eLOC_Y] / (rMinY() - tol0))));
  bool insideOuter = ((square(lposition[Acts::eLOC_X] / (rMaxX() + tol0)) +
                       square(lposition[Acts::eLOC_Y] / (rMaxY() + tol0))) < 1);
  return (insidePhi && insideInner && insideOuter);
}

// For ellipse bound this is only approximation which is valid
// only if m_valueStore.at(EllipseBounds::bv_rMinX) ~=
// m_valueStore.at(EllipseBounds::bv_rMinY)
// and m_valueStore.at(EllipseBounds::bv_rMaxX) ~=
// m_valueStore.at(EllipseBounds::bv_rMaxY)
//
double Acts::EllipseBounds::distanceToBoundary(
    const Vector2D& lposition) const {
  double r = perp(lposition);
  if (r == 0) {
    return std::min(rMinX(), rMinY());
  }

  double sn = lposition[eLOC_X] / r;
  double cs = lposition[eLOC_Y] / r;
  double dF = detail::radian_sym(phi(lposition) - m_avgPhi);
  double sf = 0.;

  if (m_halfPhi < M_PI) {
    double df = std::abs(dF) - m_halfPhi;
    sf = r * std::sin(df);
    if (df > 0.) {
      r *= std::cos(df);
    }
  } else {
    sf = -1.e+10;
  }

  if (sf <= 0.) {
    double a = cs / m_rMaxX;
    double b = sn / m_rMaxY;
    double sr0 = r - 1. / std::hypot(a, b);
    if (sr0 >= 0.) {
      return sr0;
    }
    a = cs / m_rMinX;
    b = sn / m_rMinY;
    double sr1 = 1. / std::hypot(a, b) - r;
    if (sr1 >= 0.) {
      return sr1;
    }
    if (sf < sr0) {
      sf = sr0;
    }
    if (sf < sr1) {
      sf = sr1;
    }
    return sf;
  }

  double fb;
  fb = (dF > 0.) ? (m_avgPhi + m_halfPhi) : (m_avgPhi - m_halfPhi);
  sn = sin(fb);
  cs = cos(fb);
  double a = cs / m_rMaxX;
  double b = sn / m_rMaxY;
  double sr0 = r - 1. / std::hypot(a, b);
  if (sr0 >= 0.) {
    return std::hypot(sr0, sf);
  }
  a = cs / m_rMinX;
  b = sn / m_rMinY;
  double sr1 = (1. / std::hypot(a, b)) - r;
  if (sr1 >= 0.) {
    return std::hypot(sr1, sf);
  }
  return sf;
}

std::vector<Acts::Vector2D> Acts::EllipseBounds::vertices(
    unsigned int lseg) const {
  return detail::VerticesHelper::ellispoidVertices(
      m_rMinX, m_rMinY, m_rMaxX, m_rMaxY, m_avgPhi, m_halfPhi, lseg);
}

const Acts::RectangleBounds& Acts::EllipseBounds::boundingBox() const {
  return m_boundingBox;
}

// ostream operator overload
std::ostream& Acts::EllipseBounds::toStream(std::ostream& sl) const {
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::EllipseBounds:  (innerRadiusX, innerRadiusY, outerRadiusX, "
        "outerRadiusY, hPhiSector) = ";
  sl << "(" << rMinX() << ", " << rMinY() << ", " << rMaxX() << ", " << rMaxY()
     << ", " << averagePhi() << ", " << halfPhiSector() << ")";
  sl << std::setprecision(-1);
  return sl;
}
