// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/CylinderBounds.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>

#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;

Acts::CylinderBounds::CylinderBounds(double radius, double halfZ)
    : CylinderBounds(radius, 0, M_PI, halfZ) {
  m_closed = true;
}

Acts::CylinderBounds::CylinderBounds(double radius, double halfPhi,
                                     double halfZ)
    : CylinderBounds(radius, 0, halfPhi, halfZ) {}

Acts::CylinderBounds::CylinderBounds(double radius, double averagePhi,
                                     double halfPhi, double halfZ)
    : m_radius(std::abs(radius)),
      m_avgPhi(detail::radian_sym(averagePhi)),
      m_halfPhi(std::abs(halfPhi)),
      m_halfZ(std::abs(halfZ)),
      m_closed(false) {
  if (halfPhi == M_PI) {
    m_closed = true;
  }
}

Acts::CylinderBounds* Acts::CylinderBounds::clone() const {
  return new CylinderBounds(*this);
}

Acts::SurfaceBounds::BoundsType Acts::CylinderBounds::type() const {
  return SurfaceBounds::Cylinder;
}

std::vector<TDD_real_t> Acts::CylinderBounds::valueStore() const {
  std::vector<TDD_real_t> values(CylinderBounds::bv_length);
  values[CylinderBounds::bv_radius] = m_radius;
  values[CylinderBounds::bv_averagePhi] = m_avgPhi;
  values[CylinderBounds::bv_halfPhiSector] = m_halfPhi;
  values[CylinderBounds::bv_halfZ] = m_halfZ;
  return values;
}

// Convert from (r*phi,z) to (phi,z) centered around phi0
Acts::Vector2D Acts::CylinderBounds::shifted(
    const Acts::Vector2D& lposition) const {
  return {Acts::detail::radian_sym((lposition[Acts::eLOC_RPHI] / m_radius) -
                                   m_avgPhi),
          lposition[Acts::eLOC_Z]};
}

// Jacobian from (r*phi,z) to (phi,z)
Acts::ActsSymMatrixD<2> Acts::CylinderBounds::jacobian() const {
  ActsSymMatrixD<2> j;
  j(0, eLOC_RPHI) = 1 / m_radius;
  j(0, eLOC_Z) = 0;
  j(1, eLOC_RPHI) = 0;
  j(1, eLOC_Z) = 1;
  return j;
}

bool Acts::CylinderBounds::inside(const Vector2D& lposition,
                                  const BoundaryCheck& bcheck) const {
  return bcheck.transformed(jacobian())
      .isInside(shifted(lposition), Vector2D(-m_halfPhi, -m_halfZ),
                Vector2D(m_halfPhi, m_halfZ));
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
      std::abs(perp(position) - m_radius)) {
    return false;
  } else if (checkAbsolute && m_closed) {
    return ((s_onSurfaceTolerance + addToleranceZ + m_halfZ) >=
            std::abs(position.z()));
  }
  // detailed, but slower check
  Vector2D lpos(detail::radian_sym(phi(position) - m_avgPhi), position.z());
  return bcheck.transformed(jacobian())
      .isInside(lpos, Vector2D(-m_halfPhi, -m_halfZ),
                Vector2D(m_halfPhi, m_halfZ));
}

double Acts::CylinderBounds::distanceToBoundary(
    const Acts::Vector2D& lposition) const {
  return BoundaryCheck(true).distance(shifted(lposition),
                                      Vector2D(-m_halfPhi, -m_halfZ),
                                      Vector2D(m_halfPhi, m_halfZ));
}

// ostream operator overload
std::ostream& Acts::CylinderBounds::toStream(std::ostream& sl) const {
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::CylinderBounds: (radius, averagePhi, halfPhiSector, "
        "halflengthInZ) = ";
  sl << "(" << m_radius << ", " << m_avgPhi << ", ";
  sl << m_halfPhi << ", " << m_halfZ << ")";
  sl << std::setprecision(-1);
  return sl;
}
