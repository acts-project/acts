// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Utilities/JacobianHelpers.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <iomanip>
#include <ios>
#include <ostream>
#include <string>

namespace Acts {

/// @brief Utility class for curvilinear surfaces
///
class CurvilinearSurface {
 public:
  explicit CurvilinearSurface(const Vector3& direction)
      : m_direction{direction} {}

  /// Constructor with direction vector
  CurvilinearSurface(const Vector3& position, const Vector3& direction)
      : m_position{position}, m_direction{direction} {}

  /// Return method for the surface center by reference
  /// @note the center is always recalculated in order to not keep a cache
  /// @return center position by value
  Vector3 center() const { return m_position; }

  /// Return the surface normal at a given @p position and @p direction.
  /// This method is fully generic, and valid for all surface types.
  /// @note For some surface types, the @p direction is ignored, but
  ///       it is **not safe** to pass in a zero vector!
  /// @return The normal vector at the given position and direction
  Vector3 normal() const { return m_direction; }

  /// Return method for the reference frame
  /// This is the frame in which the covariance matrix is defined (specialized
  /// by all surfaces)
  ///
  /// @return RotationMatrix3 which defines the three axes of the measurement
  /// frame
  Acts::RotationMatrix3 referenceFrame() const {
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

  /// Return method for the surface Transform3 by reference
  /// In case a detector element is associated the surface transform
  /// is just forwarded to the detector element in order to keep the
  /// (mis-)alignment cache cetrally handled
  ///
  /// @return the contextual transform
  Transform3 transform() const {
    Transform3 transform = Transform3(referenceFrame());
    transform.pretranslate(center());
    return transform;
  }

  /// Calculate the jacobian from local to global which the surface knows best,
  /// hence the calculation is done here.
  ///
  /// @return Jacobian from local to global
  virtual BoundToFreeMatrix boundToFreeJacobian() const {
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
  FreeToBoundMatrix freeToBoundJacobian() const {
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

  /// Calculate the derivative of path length at the geometry constraint or
  /// point-of-closest-approach w.r.t. free parameters. The calculation is
  /// identical for all surfaces where the reference frame does not depend on
  /// the direction
  ///
  /// @return Derivative of path length w.r.t. free parameters
  FreeToPathMatrix freeToPathDerivative() const {
    FreeToPathMatrix freeToPath = FreeToPathMatrix::Zero();
    freeToPath.segment<3>(eFreePos0) = -1.0 * m_direction;
    return freeToPath;
  }

  /// Output Method for std::ostream
  ///
  /// @param sl is the ostream to be dumped into
  std::ostream& toStream(std::ostream& sl) const {
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

  /// Output into a std::string
  std::string toString() const {
    std::stringstream ss;
    toStream(ss);
    return ss.str();
  }

  std::shared_ptr<PlaneSurface> planeSurface() const {
    return Surface::makeShared<PlaneSurface>(transform());
  }

 private:
  Vector3 m_position = Vector3::Zero();
  Vector3 m_direction = Vector3::UnitZ();
};

inline std::ostream& operator<<(std::ostream& os,
                                const CurvilinearSurface& surface) {
  surface.toStream(os);
  return os;
}

}  // namespace Acts
