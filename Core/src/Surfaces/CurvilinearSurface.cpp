// This file is part of the Acts project.
//
// Copyright (C) 2023-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/CurvilinearSurface.hpp"

#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Utilities/JacobianHelpers.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <iomanip>
#include <ios>
#include <ostream>

namespace Acts {

RotationMatrix3 CurvilinearSurface::referenceFrame() const {
  /// the right-handed coordinate system is defined as
  /// T = normal
  /// U = Z x T
  /// V = T x U
  Vector3 T = m_direction;
  Vector3 U = Vector3::UnitZ().cross(T).normalized();
  Vector3 V = T.cross(U);

  RotationMatrix3 rframe;
  rframe << U, V, T;
  return rframe;
}

Transform3 CurvilinearSurface::transform() const {
  Transform3 transform = Transform3(referenceFrame());
  transform.pretranslate(center());
  return transform;
}

/// Calculate the jacobian from local to global which the surface knows best,
/// hence the calculation is done here.
///
/// @return Jacobian from local to global
BoundToFreeMatrix CurvilinearSurface::boundToFreeJacobian() const {
  // Prepare the jacobian to free
  BoundToFreeMatrix jacobian = BoundToFreeMatrix::Zero();

  auto [cosPhi, sinPhi, cosTheta, sinTheta] =
      VectorHelpers::evaluateTrigonomics(m_direction);
  jacobian(eFreePos0, eBoundLoc0) = -sinPhi;
  jacobian(eFreePos0, eBoundLoc1) = -cosPhi * cosTheta;
  jacobian(eFreePos1, eBoundLoc0) = cosPhi;
  jacobian(eFreePos1, eBoundLoc1) = -sinPhi * cosTheta;
  jacobian(eFreePos2, eBoundLoc1) = sinTheta;
  // Time parameter: stays as is
  jacobian(eFreeTime, eBoundTime) = 1;
  jacobian.block<3, 2>(eFreeDir0, eBoundPhi) =
      sphericalToFreeDirectionJacobian(m_direction);
  // Q/P parameter: stays as is
  jacobian(eFreeQOverP, eBoundQOverP) = 1;

  return jacobian;
}

/// Calculate the jacobian from global to local which the surface knows best,
/// hence the calculation is done here.
///
/// @return Jacobian from global to local
FreeToBoundMatrix CurvilinearSurface::freeToBoundJacobian() const {
  // Prepare the jacobian to curvilinear
  FreeToBoundMatrix jacobian = FreeToBoundMatrix::Zero();

  auto [cosPhi, sinPhi, cosTheta, sinTheta] =
      VectorHelpers::evaluateTrigonomics(m_direction);
  // We operate in curvilinear coordinates defined as follows
  jacobian(eBoundLoc0, eFreePos0) = -sinPhi;
  jacobian(eBoundLoc0, eFreePos1) = cosPhi;
  jacobian(eBoundLoc1, eFreePos0) = -cosPhi * cosTheta;
  jacobian(eBoundLoc1, eFreePos1) = -sinPhi * cosTheta;
  jacobian(eBoundLoc1, eFreePos2) = sinTheta;
  // Time parameter: stays as is
  jacobian(eBoundTime, eFreeTime) = 1.;
  // Directional and momentum parameters for curvilinear
  jacobian.block<2, 3>(eBoundPhi, eFreeDir0) =
      freeToSphericalDirectionJacobian(m_direction);
  // Q/P parameter: stays as is
  jacobian(eBoundQOverP, eFreeQOverP) = 1.;

  return jacobian;
}

FreeToPathMatrix CurvilinearSurface::freeToPathDerivative() const {
  FreeToPathMatrix freeToPath = FreeToPathMatrix::Zero();
  freeToPath.segment<3>(eFreePos0) = -1.0 * m_direction;
  return freeToPath;
}

std::ostream& CurvilinearSurface::toStream(std::ostream& sl) const {
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(4);
  sl << std::boolalpha;
  sl << "Curvilinear Surface" << std::endl;
  sl << "     Center position  (x, y, z) = (" << m_position.x() << ", "
     << m_position.y() << ", " << m_position.z() << ")" << std::endl;
  sl << "     Direction  (x, y, z) = (" << m_direction.x() << ", "
     << m_direction.y() << ", " << m_direction.z() << ")" << std::endl;
  sl << std::setprecision(-1);
  sl << std::noboolalpha;
  return sl;
}

std::string CurvilinearSurface::toString() const {
  std::stringstream ss;
  toStream(ss);
  return ss.str();
}

std::shared_ptr<PlaneSurface> CurvilinearSurface::planeSurface() const {
  return Surface::makeShared<PlaneSurface>(transform());
}

std::shared_ptr<Surface> CurvilinearSurface::surface() const {
  return planeSurface();
}

std::ostream& operator<<(std::ostream& os, const CurvilinearSurface& surface) {
  surface.toStream(os);
  return os;
}

}  // namespace Acts
