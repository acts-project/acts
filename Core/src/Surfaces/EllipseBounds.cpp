// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/EllipseBounds.hpp"
#include "Acts/Surfaces/detail/VerticesHelper.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>

using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;

Acts::SurfaceBounds::BoundsType Acts::EllipseBounds::type() const {
  return SurfaceBounds::eEllipse;
}

static inline double square(double x) {
  return x * x;
}

/// @warning This **only** works for tolerance-based checks
bool Acts::EllipseBounds::inside(const Vector2D& lposition,
                                 const BoundaryCheck& bcheck) const {
  double tol0 = bcheck.m_tolerance[0];
  double tol1 = bcheck.m_tolerance[1];
  double phi =
      detail::radian_sym(VectorHelpers::phi(lposition) - get(eAveragePhi));
  double phiHalf = get(eHalfPhiSector) + tol1;

  bool insidePhi = (-phiHalf <= phi) && (phi < phiHalf);
  bool insideInner =
      (get(eInnerRx) <= tol0) || (get(eOuterRx) <= tol0) ||
      (1 < (square(lposition[Acts::eLOC_X] / (get(eInnerRx) - tol0)) +
            square(lposition[Acts::eLOC_Y] / (get(eOuterRx) - tol0))));
  bool insideOuter =
      ((square(lposition[Acts::eLOC_X] / (get(eInnerRy) + tol0)) +
        square(lposition[Acts::eLOC_Y] / (get(eOuterRy) + tol0))) < 1);
  return (insidePhi && insideInner && insideOuter);
}

// For ellipse bound this is only approximation which is valid
// only if m_values.at(EllipseBounds::bv_rMinX) ~=
// m_values.at(EllipseBounds::bv_rMinY)
// and m_values.at(EllipseBounds::bv_rMaxX) ~=
// m_values.at(EllipseBounds::bv_rMaxY)
//
double Acts::EllipseBounds::distanceToBoundary(
    const Vector2D& lposition) const {
  double r = perp(lposition);
  if (r == 0) {
    return std::min(get(eInnerRx), get(eOuterRx));
  }

  double sn = lposition[eLOC_X] / r;
  double cs = lposition[eLOC_Y] / r;
  double dF = detail::radian_sym(phi(lposition) - get(eAveragePhi));
  double sf = 0.;

  if (get(eHalfPhiSector) < M_PI) {
    double df = std::abs(dF) - get(eHalfPhiSector);
    sf = r * std::sin(df);
    if (df > 0.) {
      r *= std::cos(df);
    }
  } else {
    sf = -1.e+10;
  }

  if (sf <= 0.) {
    double a = cs / get(eInnerRy);
    double b = sn / get(eOuterRy);
    double sr0 = r - 1. / std::hypot(a, b);
    if (sr0 >= 0.) {
      return sr0;
    }
    a = cs / get(eInnerRx);
    b = sn / get(eOuterRx);
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
  fb = (dF > 0.) ? (get(eAveragePhi) + get(eHalfPhiSector))
                 : (get(eAveragePhi) - get(eHalfPhiSector));
  sn = sin(fb);
  cs = cos(fb);
  double a = cs / get(eInnerRy);
  double b = sn / get(eOuterRy);
  double sr0 = r - 1. / std::hypot(a, b);
  if (sr0 >= 0.) {
    return std::hypot(sr0, sf);
  }
  a = cs / get(eInnerRx);
  b = sn / get(eOuterRx);
  double sr1 = (1. / std::hypot(a, b)) - r;
  if (sr1 >= 0.) {
    return std::hypot(sr1, sf);
  }
  return sf;
}

std::vector<Acts::Vector2D> Acts::EllipseBounds::vertices(
    unsigned int lseg) const {
  return detail::VerticesHelper::ellispoidVertices(
      get(eInnerRx), get(eInnerRy), get(eOuterRx), get(eOuterRy),
      get(eAveragePhi), get(eHalfPhiSector), lseg);
}

const Acts::RectangleBounds& Acts::EllipseBounds::boundingBox() const {
  return m_boundingBox;
}

// ostream operator overload
std::ostream& Acts::EllipseBounds::toStream(std::ostream& sl) const {
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::EllipseBounds:  (innerRadius0, outerRadius0, innerRadius1, "
        "outerRadius1, hPhiSector, averagePhi) = ";
  sl << "(" << get(eInnerRx) << ", " << get(eInnerRy) << ", " << get(eOuterRx)
     << ", " << get(eOuterRy) << ", " << get(eAveragePhi) << ", "
     << get(eHalfPhiSector) << ", " << get(eAveragePhi) << ")";
  sl << std::setprecision(-1);
  return sl;
}
