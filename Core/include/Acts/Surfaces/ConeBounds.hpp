// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <cfloat>

#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

///  @class ConeBounds
///
///  Bounds for a conical surface,
///  the opening angle is stored in \f$ \tan(\alpha) \f$ and always positively
/// defined.
///  The cone can open to both sides, steered by \f$ z_min \f$ and \f$ z_max
///  \f$.
///
///  @image html ConeBounds.gif
///

class ConeBounds : public SurfaceBounds {
 public:
  /// @enum BoundValues for readablility
  enum BoundValues : int {
    eAlpha = 0,
    eMinZ = 1,
    eMaxZ = 2,
    eHalfPhiSector = 3,
    eAveragePhi = 4,
    eValues = 5
  };

  ConeBounds() = delete;

  /// Constructor - open cone with alpha, by default a full cone
  /// but optionally can make a conical section
  ///
  /// @param alpha is the opening angle of the cone
  /// @param symm is the boolean indicating if the cone is symmetric in +/- z
  /// @param halfphi is the half opening angle (default is pi)
  /// @param avphi is the phi value around which the bounds are opened
  /// (default=0)
  ConeBounds(double alpha, bool symm, double halfphi = M_PI, double avphi = 0.);

  /// Constructor - open cone with alpha, minz and maxz, by
  /// default a full cone but can optionally make it a conical section
  ///
  /// @param alpha is the opening angle of the cone
  /// @param zmin cone expanding from minimal z
  /// @param zmax cone expanding to maximal z
  /// @param halfphi is the half opening angle (default is pi)
  /// @param avphi is the phi value around which the bounds are opened
  /// (default=0)
  ConeBounds(double alpha, double zmin, double zmax, double halfphi = M_PI,
             double avphi = 0.);

  /// Constructor - from parameters vector
  ///
  /// @param parametes The parameter vector
  ConeBounds(const std::vector<double>& parameters);

  ~ConeBounds() override = default;

  ConeBounds* clone() const final;

  /// Return the bounds type
  BoundsType type() const final;

  /// Return the bound values
  ///
  /// @return this returns a copy of the internal parameters
  std::vector<double> boundValues() const final;

  /// inside method for local position
  ///
  /// @param lposition is the local position to be checked
  /// @param bcheck is the boundary check directive
  /// @return is a boolean indicating if the position is inside
  bool inside(const Vector2D& lposition,
              const BoundaryCheck& bcheck = true) const final;

  /// Minimal distance to boundary ( > 0 if outside and <=0 if inside)
  ///
  /// @param lposition is the local position to check for the distance
  /// @return is a signed distance parameter
  double distanceToBoundary(const Vector2D& lposition) const final;

  /// Output Method for std::ostream
  ///
  /// @param sl is the ostrea into which the dump is done
  /// @return is the input obect
  std::ostream& toStream(std::ostream& sl) const final;

  /// Return the radius at a specific z values
  ///
  /// @param z is the z value for which r is requested
  /// @return is the r value associated with z
  double r(double z) const;

  /// Return tangent of alpha (pre-computed)
  double tanAlpha() const;

  /// Templated access to the bound parameters
  template <BoundValues bValue>
  double get() const {
    return m_parameters[bValue];
  }

 private:
  std::vector<double> m_parameters;
  double m_tanAlpha;

  /// Private helper functin to shift a local 2D position
  ///
  /// @param lposition The original local position
  Vector2D shifted(const Vector2D& lposition) const;
};

inline double ConeBounds::r(double z) const {
  return std::abs(z * m_tanAlpha);
}

inline double ConeBounds::tanAlpha() const {
  return m_tanAlpha;
}

inline std::vector<double> ConeBounds::boundValues() const {
  return m_parameters;
}

}  // namespace Acts