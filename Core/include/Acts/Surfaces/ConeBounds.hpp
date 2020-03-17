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
  enum BoundValues {
    bv_alpha = 0,
    bv_minZ = 1,
    bv_maxZ = 2,
    bv_averagePhi = 3,
    bv_halfPhiSector = 4,
    bv_length = 5
  };

  // Deleted default constructor
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

  /// Defaulted destructor
  ~ConeBounds() override = default;

  /// Virtual constructor
  ConeBounds* clone() const final;

  /// The type enumeration
  BoundsType type() const final;

  /// The value store for persistency
  std::vector<TDD_real_t> valueStore() const final;

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

  /// Return the average values for the angles
  double tanAlpha() const;

  /// Return the average values for the angles
  double sinAlpha() const;

  /// Return the average values for the angles
  double cosAlpha() const;

  /// Return the average values for the angles
  double alpha() const;

  /// This method returns the minimum z value in the local
  /// frame for an unbound symmetric cone, it returns -MAXBOUNDVALUE*/
  double minZ() const;

  /// This method returns the maximum z value in the local
  /// frame for an unbound symmetric cone, it returns -MAXBOUNDVALUE*/
  double maxZ() const;

  /// This method returns the average phi value
  /// (i.e. the "middle" phi value for the conical sector we  are describing)
  double averagePhi() const;

  /// This method returns the half-phi width of the sector
  /// (so that averagePhi +/- halfPhiSector gives the phi bounds of the cone)
  double halfPhiSector() const;

 private:
  double m_alpha, m_tanAlpha;
  double m_zMin, m_zMax;
  double m_avgPhi, m_halfPhi;

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

inline double ConeBounds::sinAlpha() const {
  return std::sin(m_alpha);
}

inline double ConeBounds::cosAlpha() const {
  return std::cos(m_alpha);
}

inline double ConeBounds::alpha() const {
  return m_alpha;
}

inline double ConeBounds::minZ() const {
  return m_zMin;
}

inline double ConeBounds::maxZ() const {
  return m_zMax;
}

inline double ConeBounds::averagePhi() const {
  return m_avgPhi;
}

inline double ConeBounds::halfPhiSector() const {
  return m_halfPhi;
}
}  // namespace Acts