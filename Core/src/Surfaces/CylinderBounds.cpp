// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>

using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;

Acts::SurfaceBounds::BoundsType Acts::CylinderBounds::type() const {
  return SurfaceBounds::eCylinder;
}

Acts::Vector2D Acts::CylinderBounds::shifted(
    const Acts::Vector2D& lposition) const {
  return {Acts::detail::radian_sym((lposition[Acts::eLOC_RPHI] / get(eR)) -
                                   get(eAveragePhi)),
          lposition[Acts::eLOC_Z]};
}

Acts::ActsSymMatrixD<2> Acts::CylinderBounds::jacobian() const {
  ActsSymMatrixD<2> j;
  j(0, eLOC_RPHI) = 1 / get(eR);
  j(0, eLOC_Z) = 0;
  j(1, eLOC_RPHI) = 0;
  j(1, eLOC_Z) = 1;
  return j;
}

bool Acts::CylinderBounds::inside(const Vector2D& lposition,
                                  const BoundaryCheck& bcheck) const {
  return bcheck.transformed(jacobian())
      .isInside(shifted(lposition),
                Vector2D(-get(eHalfPhiSector), -get(eHalfLengthZ)),
                Vector2D(get(eHalfPhiSector), get(eHalfLengthZ)));
}

bool Acts::CylinderBounds::inside3D(const Vector3D& position,
                                    const BoundaryCheck& bcheck) const {
  // additional tolerance from the boundary check if configred
  bool checkAbsolute = bcheck.m_type == BoundaryCheck::Type::eAbsolute;

  // this fast check only applies to closed cylindrical bounds
  double addToleranceR =
      (checkAbsolute && m_closed) ? bcheck.m_tolerance[0] : 0.;
  double addToleranceZ = checkAbsolute ? bcheck.m_tolerance[1] : 0.;
  // check if the position compatible with the radius
  if ((s_onSurfaceTolerance + addToleranceR) <=
      std::abs(perp(position) - get(eR))) {
    return false;
  } else if (checkAbsolute && m_closed) {
    return ((s_onSurfaceTolerance + addToleranceZ + get(eHalfLengthZ)) >=
            std::abs(position.z()));
  }
  // detailed, but slower check
  Vector2D lpos(detail::radian_sym(phi(position) - get(eAveragePhi)),
                position.z());
  return bcheck.transformed(jacobian())
      .isInside(lpos, Vector2D(-get(eHalfPhiSector), -get(eHalfLengthZ)),
                Vector2D(get(eHalfPhiSector), get(eHalfLengthZ)));
}

double Acts::CylinderBounds::distanceToBoundary(
    const Acts::Vector2D& lposition) const {
  return BoundaryCheck(true).distance(
      shifted(lposition), Vector2D(-get(eHalfPhiSector), -get(eHalfLengthZ)),
      Vector2D(get(eHalfPhiSector), get(eHalfLengthZ)));
}

std::ostream& Acts::CylinderBounds::toStream(std::ostream& sl) const {
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::CylinderBounds: (radius, halfLengthZ, halfPhiSector, "
        "averagePhi) = ";
  sl << "(" << get(eR) << ", " << get(eHalfLengthZ) << ", ";
  sl << get(eHalfPhiSector) << ", " << get(eAveragePhi) << ")";
  sl << std::setprecision(-1);
  return sl;
}
