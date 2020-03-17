// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Surfaces/DiscBounds.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

namespace Acts {

/// @brief Class that implements a (potentially asymmetric) bounds with
/// difference between surface bound center and surface coordinate center
///
/// These bounds combine two different systems:
///  * module system : radial bounds centred on the moduleOrigin
///  * strip system : phi bounds centred on the stripOrigin
///
/// The measurement will be done in the strip system, with r/phi local
/// coordinates.
///
class AnnulusBounds : public DiscBounds {
 public:
  using Transform2D = Eigen::Affine2d;

  /// enumeration for readability
  enum BoundValues {
    bv_minR = 0,
    bv_maxR = 1,
    bv_phiMin = 2,
    bv_phiMax = 3,
    bv_phiAvg = 4,
    bv_originX = 5,
    bv_originY = 6,
    bv_length = 7
  };

  /// @brief Default constructor from parameters
  /// @param minR inner radius, in module system
  /// @param maxR outer radius, in module system
  /// @param minPhi right angular edge, in strip system
  /// @param maxPhi left angular edge, in strip system
  /// @param moduleOrigin The origin offset between the two systems.
  /// @param avgPhi (Optional) internal rotation of this bounds object's local
  /// frame
  /// @note For @c moduleOrigin you need to actually calculate the cartesian
  /// offset
  AnnulusBounds(double minR, double maxR, double minPhi, double maxPhi,
                const Vector2D& moduleOrigin = {0, 0}, double avgPhi = 0);

  /// @brief copy constructor
  AnnulusBounds(const AnnulusBounds& source) = default;

  /// Virtual constructor
  AnnulusBounds* clone() const final;

  /// Bound type
  SurfaceBounds::BoundsType type() const final;

  /// This returns the stored values for persisitency
  std::vector<TDD_real_t> valueStore() const final;

  /// Inside check for the bounds object driven by the boundary check directive
  /// Each Bounds has a method inside, which checks if a LocalPosition is inside
  /// the bounds  Inside can be called without/with tolerances.
  ///
  /// @param lposition Local position (assumed to be in right surface frame)
  /// @param bcheck boundary check directive
  /// @return boolean indicator for the success of this operation
  virtual bool inside(const Vector2D& lposition,
                      const BoundaryCheck& bcheck) const final;

  /// Minimal distance to boundary ( > 0 if outside and <=0 if inside)
  ///
  /// @param lposition is the local position to check for the distance
  /// @return is a signed distance parameter
  virtual double distanceToBoundary(const Vector2D& lposition) const final;

  /// Outstream operator
  ///
  /// @param sl is the ostream to be dumped into
  std::ostream& toStream(std::ostream& sl) const final;

  /// @brief Returns inner radial bounds (module system)
  /// @return The inner radius
  double rMin() const;

  /// @brief Returns outer radial bounds (module system)
  /// @return The outer radius
  double rMax() const;

  /// @brief Returns the right angular edge of the module
  /// @return The right side angle
  double phiMin() const;

  /// @brief Returns the left angular edge of the module
  /// @return The left side angle
  double phiMax() const;

  /// Returns true for full phi coverage
  bool coversFullAzimuth() const final;

  /// Checks if this is inside the radial coverage
  /// given the a tolerance
  bool insideRadialBounds(double R, double tolerance = 0.) const final;

  /// Return a reference radius for binning
  double binningValueR() const final;

  /// Return a reference radius for binning
  double binningValuePhi() const final;

  /// @brief Returns moduleOrigin, but rotated out, so @c avgPhi is already
  /// considered. The module origin needs to consider the rotation introduced by
  /// @c avgPhi
  /// @return The origin of the local frame
  Vector2D moduleOrigin() const;

  /// This method returns the four corners of the bounds in polar coordinates
  /// Starting from the upper right (max R, pos locX) and proceding clock-wise
  /// i.e. (max R; pos locX), (min R; pos locX), (min R; neg loc X), (max R: neg
  /// locX)
  std::vector<Vector2D> corners() const;

  /// This method returns the xy coordinates of the four corners of the
  /// bounds in module coorindates (in x/y)
  /// Starting from the upper right (max R, pos locX) and proceding clock-wise
  /// i.e. (max R; pos locX), (min R; pos locX), (min R; neg loc X), (max R: neg
  /// locX)
  ///
  ///
  /// @param lseg the number of segments used to approximate
  /// and eventually curved line
  ///
  /// @note that the extremas are given, which may slightly alter the
  /// number of segments returned
  ///
  /// @return vector for vertices in 2D
  std::vector<Vector2D> vertices(unsigned int lseg) const;

 private:
  double m_rMin;
  double m_rMax;
  double m_phiMin;
  double m_phiMax;

  // @TODO: Does this need to be in bound values?
  Vector2D m_moduleOrigin;
  Vector2D m_shiftXY;  // == -m_moduleOrigin
  Vector2D m_shiftPC;
  double m_phiAvg;
  Transform2D m_rotationStripPC;
  Transform2D m_translation;

  // Vectors needed for inside checking
  Vector2D m_outLeftStripPC;
  Vector2D m_inLeftStripPC;
  Vector2D m_outRightStripPC;
  Vector2D m_inRightStripPC;

  Vector2D m_outLeftModulePC;
  Vector2D m_inLeftModulePC;
  Vector2D m_outRightModulePC;
  Vector2D m_inRightModulePC;

  Vector2D m_outLeftStripXY;
  Vector2D m_inLeftStripXY;
  Vector2D m_outRightStripXY;
  Vector2D m_inRightStripXY;

  /// Private helper method:
  /// Inside check for the bounds object driven by the boundary check directive
  /// Each Bounds has a method inside, which checks if a LocalPosition is inside
  /// the bounds  Inside can be called without/with tolerances.
  ///
  /// @param lposition Local position (assumed to be in right surface frame)
  /// @param tolR tolerance on the radius
  /// @param tolPhi tolerance on the polar angle phi
  /// @return boolean indicator for the success of this operation
  virtual bool inside(const Vector2D& lposition, double tolR,
                      double tolPhi) const final;

  /// Private helper mehtod
  Vector2D stripXYToModulePC(const Vector2D& vStripXY) const;

  /// Private helper mehtod
  Vector2D closestOnSegment(const Vector2D& a, const Vector2D& b,
                            const Vector2D& p,
                            const Eigen::Matrix<double, 2, 2>& weight) const;

  /// Private helper mehtod
  double squaredNorm(const Vector2D& v,
                     const Eigen::Matrix<double, 2, 2>& weight) const;
};

inline AnnulusBounds* AnnulusBounds::clone() const {
  return new AnnulusBounds(m_rMin, m_rMax, m_phiMin, m_phiMax, m_moduleOrigin,
                           m_phiAvg);
}

inline SurfaceBounds::BoundsType AnnulusBounds::type() const {
  return SurfaceBounds::Annulus;
}

inline double AnnulusBounds::rMin() const {
  return m_rMin;
}

inline double AnnulusBounds::rMax() const {
  return m_rMax;
}

inline double AnnulusBounds::phiMin() const {
  return m_phiMin + m_phiAvg;
}

inline double AnnulusBounds::phiMax() const {
  return m_phiMax + m_phiAvg;
}

inline bool AnnulusBounds::coversFullAzimuth() const {
  return (std::abs((m_phiMax - m_phiMin) - M_PI) < s_onSurfaceTolerance);
}

inline bool AnnulusBounds::insideRadialBounds(double R,
                                              double tolerance) const {
  return ((R + tolerance) > m_rMin and (R - tolerance) < m_rMax);
}

inline double AnnulusBounds::binningValueR() const {
  return 0.5 * (m_rMin + m_rMax);
}

inline double AnnulusBounds::binningValuePhi() const {
  return m_phiAvg;
}

}  // namespace Acts
