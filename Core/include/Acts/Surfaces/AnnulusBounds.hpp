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

#include <array>
#include <exception>
#include <vector>

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

  enum BoundValues : int {
    eMinR = 0,
    eMaxR = 1,
    eMinPhiRel = 2,
    eMaxPhiRel = 3,
    eAveragePhi = 4,
    eOriginX = 5,
    eOriginY = 6,
    eSize = 7
  };

  AnnulusBounds() = delete;

  /// @brief Default constructor from parameters
  /// @param minR inner radius, in module system
  /// @param maxR outer radius, in module system
  /// @param minPhiRel right angular edge, in strip system, rel to avgOhi
  /// @param maxPhiRel left angular edge, in strip system, rel to avgPhi
  /// @param moduleOrigin The origin offset between the two systems.
  /// @param avgPhi (Optional) internal rotation of this bounds object's local
  /// frame
  /// @note For @c morigin you need to actually calculate the cartesian
  /// offset
  AnnulusBounds(double minR, double maxR, double minPhiRel, double maxPhiRel,
                const Vector2D& moduleOrigin = {0, 0},
                double avgPhi = 0) noexcept(false)
      : AnnulusBounds({minR, maxR, minPhiRel, maxPhiRel, avgPhi,
                       moduleOrigin.x(), moduleOrigin.y()}) {}

  /// Constructor - from parameters array
  ///
  /// @param values The parameter array
  AnnulusBounds(const std::array<double, eSize>& values) noexcept(false);

  AnnulusBounds(const AnnulusBounds& source) = default;

  SurfaceBounds::BoundsType type() const final;

  /// Return the bound values as dynamically sized vector
  ///
  /// @return this returns a copy of the internal values
  std::vector<double> values() const final;

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

  /// Access to the bound values
  /// @param bValue the class nested enum for the array access
  double get(BoundValues bValue) const { return m_values[bValue]; }

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

  /// @brief Returns moduleOrigin, but rotated out, so @c averagePhi is already
  /// considered. The module origin needs to consider the rotation introduced by
  /// @c averagePhi
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
  /// @param lseg the number of segments used to approximate
  /// and eventually curved line
  ///
  /// @note that the extremas are given, which may slightly alter the
  /// number of segments returned
  ///
  /// @return vector for vertices in 2D
  std::vector<Vector2D> vertices(unsigned int lseg) const;

  /// This method returns inner radius
  double rMin() const final;

  /// This method returns outer radius
  double rMax() const final;

 private:
  std::array<double, eSize> m_values;

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

  /// Check the input values for consistency, will throw a logic_exception
  /// if consistency is not given
  void checkConsistency() noexcept(false);

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

  /// Transform the strip cartesien
  /// into the module polar system
  ///
  /// @param vStripXy the position in the cartesian strip system
  /// @return the poistion in the module polar coordiante system
  Vector2D stripXYToModulePC(const Vector2D& vStripXY) const;

  /// Private helper method
  Vector2D closestOnSegment(const Vector2D& a, const Vector2D& b,
                            const Vector2D& p,
                            const Eigen::Matrix<double, 2, 2>& weight) const;

  /// Private helper mehtod
  double squaredNorm(const Vector2D& v,
                     const Eigen::Matrix<double, 2, 2>& weight) const;
};

inline SurfaceBounds::BoundsType AnnulusBounds::type() const {
  return SurfaceBounds::eAnnulus;
}

inline double AnnulusBounds::rMin() const {
  return get(eMinR);
}

inline double AnnulusBounds::rMax() const {
  return get(eMaxR);
}

inline double AnnulusBounds::phiMin() const {
  return get(eMinPhiRel) + get(eAveragePhi);
}

inline double AnnulusBounds::phiMax() const {
  return get(eMaxPhiRel) + get(eAveragePhi);
}

inline bool AnnulusBounds::coversFullAzimuth() const {
  return (std::abs((get(eMinPhiRel) - get(eMaxPhiRel)) - M_PI) <
          s_onSurfaceTolerance);
}

inline bool AnnulusBounds::insideRadialBounds(double R,
                                              double tolerance) const {
  return ((R + tolerance) > get(eMinR) and (R - tolerance) < get(eMaxR));
}

inline double AnnulusBounds::binningValueR() const {
  return 0.5 * (get(eMinR) + get(eMaxR));
}

inline double AnnulusBounds::binningValuePhi() const {
  return get(eAveragePhi);
}

inline std::vector<double> AnnulusBounds::values() const {
  std::vector<double> valvector;
  valvector.insert(valvector.begin(), m_values.begin(), m_values.end());
  return valvector;
}

inline void AnnulusBounds::checkConsistency() noexcept(false) {
  if (get(eMinR) < 0. or get(eMaxR) < 0. or get(eMinR) > get(eMaxR) or
      std::abs(get(eMinR) - get(eMaxR)) < s_epsilon) {
    throw std::invalid_argument("AnnulusBounds: invalid radial setup.");
  }
  if (get(eMinPhiRel) != detail::radian_sym(get(eMinPhiRel)) or
      get(eMaxPhiRel) != detail::radian_sym(get(eMaxPhiRel)) or
      get(eMinPhiRel) > get(eMaxPhiRel)) {
    throw std::invalid_argument("AnnulusBounds: invalid phi boundary setup.");
  }
  if (get(eAveragePhi) != detail::radian_sym(get(eAveragePhi))) {
    throw std::invalid_argument("AnnulusBounds: invalid phi positioning.");
  }
}

}  // namespace Acts
