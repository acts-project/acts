// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <cmath>

#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

namespace Acts {

/// @class CylinderBounds
///
/// Bounds for a cylindrical Surface.
///
/// These bounds may be used for a CylinderSurface
/// In case of bounds for a StraightLineSurface the radius determines the radius
/// within a localPosition
/// is regarded as inside bounds.
///
/// CylinderBounds also enhance the possibility of a cylinder segment with an
/// opening angle @f$ 2\cdot\phi_{half}@f$
/// around an average @f$ \phi @f$ angle @f$ \phi_{ave} @f$.
///
/// @todo update the documentation picture for cylinder segments
///
/// @image html CylinderBounds.gif

class CylinderBounds : public SurfaceBounds {
 public:
  /// @enum BoundValues for readablility
  /// nested enumeration object
  enum BoundValues : int {
    eRadius = 0,
    eHalfLengthZ = 1,
    eHalfPhiSector = 2,
    eAveragePhi = 3,
    eSize = 4
  };

  CylinderBounds() = delete;

  /// Constructor - full cylinder
  ///
  /// @param radius The radius of the cylinder
  /// @param halfz The half length in z
  /// @param halfphi The half opening angle
  /// @param avphi The phi value from which the opening angle spans (both sides)
  CylinderBounds(double radius, double halfz, double halfphi = M_PI,
                 double avphi = 0.);

  ~CylinderBounds() override = default;

  CylinderBounds* clone() const final;

  /// Type enumeration
  BoundsType type() const final;

  /// The values of the object as dynamically sized vector
  ActsVectorXd values() const final;

  /// Inside check for the bounds object driven by the boundary check directive
  /// Each Bounds has a method inside, which checks if a LocalPosition is inside
  /// the bounds  Inside can be called without/with tolerances.
  ///
  /// @param lposition Local position (assumed to be in right surface frame)
  /// @param bcheck boundary check directive
  /// @return boolean indicator for the success of this operation
  bool inside(const Vector2D& lposition,
              const BoundaryCheck& bcheck) const final;

  /// Specialized method for CylinderBounds that checks if a global position
  /// is within the the cylinder cover
  ///
  /// @param position is the position in the cylinder frame
  /// @param bcheck is the boundary check directive
  /// @return boolean indicator for operation success
  bool inside3D(const Vector3D& position,
                const BoundaryCheck& bcheck = true) const;

  /// Minimal distance to boundary ( > 0 if outside and <=0 if inside)
  ///
  /// @param lposition is the local position to check for the distance
  /// @return is a signed distance parameter
  double distanceToBoundary(const Vector2D& lposition) const final;

  /// Templated access to the bound parameters
  double get(BoundValues bValue) const { return m_parameters[bValue]; }

  /// Returns true for full phi coverage
  bool coversFullAzimuth() const;

  /// Output Method for std::ostream
  std::ostream& toStream(std::ostream& sl) const final;

 private:
  /// The bound radius, half Z, half phi and average phi
  ActsRowVectorD<eSize> m_parameters;
  /// Indicator if the bounds are closed
  bool m_closed;

  Vector2D shifted(const Vector2D& lposition) const;
  ActsSymMatrixD<2> jacobian() const;
};

inline ActsVectorXd CylinderBounds::values() const {
  return m_parameters;
}

inline bool CylinderBounds::coversFullAzimuth() const {
  return m_closed;
}

}  // namespace Acts