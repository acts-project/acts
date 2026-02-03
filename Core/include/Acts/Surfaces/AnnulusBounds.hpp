// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Surfaces/DiscBounds.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"

#include <array>
#include <cmath>
#include <iosfwd>
#include <numbers>
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
  /// Enumeration for the different bound values
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

  /// @brief Default constructor from parameters
  /// @param minR The inner radius of the annulus
  /// @param maxR The outer radius of the annulus
  /// @param minPhiRel The minimum phi relative to average phi
  /// @param maxPhiRel The maximum phi relative to average phi
  /// @param moduleOrigin The origin of the module in the strip frame
  /// @param avgPhi The average phi value
  /// @note For @c morigin you need to actually calculate the cartesian
  /// offset
  explicit AnnulusBounds(double minR, double maxR, double minPhiRel,
                         double maxPhiRel, const Vector2& moduleOrigin = {0, 0},
                         double avgPhi = 0) noexcept(false)
      : AnnulusBounds({minR, maxR, minPhiRel, maxPhiRel, avgPhi,
                       moduleOrigin.x(), moduleOrigin.y()}) {}

  /// Constructor - from fixed size array
  ///
  /// @param values The bound values stored in a fixed size array
  explicit AnnulusBounds(const std::array<double, eSize>& values) noexcept(
      false);

  BoundsType type() const final { return Annulus; }

  /// @copydoc SurfaceBounds::isCartesian
  bool isCartesian() const final { return false; }

  /// @copydoc SurfaceBounds::boundToCartesianJacobian
  SquareMatrix2 boundToCartesianJacobian(const Vector2& lposition) const final;

  /// @copydoc SurfaceBounds::boundToCartesianMetric
  SquareMatrix2 boundToCartesianMetric(const Vector2& lposition) const final;

  /// Return the bound values as dynamically sized vector
  /// @return this returns a copy of the internal values
  std::vector<double> values() const final;

  /// @copydoc SurfaceBounds::inside
  bool inside(const Vector2& lposition) const final;

  /// @copydoc SurfaceBounds::closestPoint
  Vector2 closestPoint(const Vector2& lposition,
                       const SquareMatrix2& metric) const final;

  using SurfaceBounds::inside;

  /// @copydoc SurfaceBounds::center
  /// @note For AnnulusBounds: returns pre-calculated center from corner vertices in strip polar coordinates (r, phi), accounting for average phi rotation
  Vector2 center() const final;

  /// Outstream operator
  /// @param sl is the ostream to be dumped into
  /// @return Reference to the output stream
  std::ostream& toStream(std::ostream& sl) const final;

  /// Access to the bound values
  /// @param bValue the class nested enum for the array access
  /// @return The value of the specified bound parameter
  double get(BoundValues bValue) const { return m_values[bValue]; }

  /// @brief Returns the right angular edge of the module
  /// @return The right side angle
  double phiMin() const { return get(eMinPhiRel) + get(eAveragePhi); }

  /// @brief Returns the left angular edge of the module
  /// @return The left side angle
  double phiMax() const { return get(eMaxPhiRel) + get(eAveragePhi); }

  /// Returns true for full phi coverage
  /// @return True if the annulus covers the full azimuthal range, false otherwise
  bool coversFullAzimuth() const final {
    return (std::abs((get(eMinPhiRel) - get(eMaxPhiRel)) - std::numbers::pi) <
            s_onSurfaceTolerance);
  }

  /// Checks if this is inside the radial coverage
  /// given the a tolerance
  /// @param R The radius value to check
  /// @param tolerance The tolerance for the check
  /// @return True if the radius is within bounds (plus tolerance), false otherwise
  bool insideRadialBounds(double R, double tolerance = 0.) const final {
    return ((R + tolerance) > get(eMinR) && (R - tolerance) < get(eMaxR));
  }

  /// Return a reference radius for binning
  /// @return Average radius for binning purposes
  double binningValueR() const final { return 0.5 * (get(eMinR) + get(eMaxR)); }

  /// Return a reference phi for binning
  /// @return Average phi angle for binning purposes
  double binningValuePhi() const final { return get(eAveragePhi); }

  /// @brief Returns moduleOrigin, but rotated out, so @c averagePhi is already
  /// considered. The module origin needs to consider the rotation introduced by
  /// @c averagePhi
  /// @return The origin of the local frame
  Vector2 moduleOrigin() const;

  /// This method returns the four corners of the bounds in polar coordinates
  /// Starting from the upper right (max R, pos locX) and proceeding clock-wise
  /// i.e. (max R; pos locX), (min R; pos locX), (min R; neg loc X), (max R: neg
  /// locX)
  /// @return Vector of corner points in polar coordinates
  std::vector<Vector2> corners() const;

  /// This method returns the xy coordinates of the four corners of the
  /// bounds in module coordinates (in x/y), and if quarterSegments is bigger or
  /// equal to 0, the curved part of the segment is included and approximated
  /// by the corresponding number of segments.
  ///
  /// Starting from the upper right (max R, pos locX) and proceeding clock-wise
  /// i.e. (max R; pos locX), (min R; pos locX), (min R; neg loc X), (max R: neg
  /// locX)
  ///
  /// @param quarterSegments the number of segments used to approximate
  /// a quarter of a circle
  ///
  /// @return vector for vertices in 2D
  std::vector<Vector2> vertices(
      unsigned int quarterSegments = 2u) const override;

  /// This method returns inner radius
  /// @return Minimum radius of the annulus
  double rMin() const final { return get(eMinR); }

  /// This method returns outer radius
  /// @return Maximum radius of the annulus
  double rMax() const final { return get(eMaxR); }

 private:
  std::array<double, eSize> m_values;

  // @TODO: Does this need to be in bound values?
  Vector2 m_moduleOrigin{
      Vector2::Zero()};  ///< The origin of the module in the strip frame
  Vector2 m_shiftXY{Vector2::Zero()};  // == -m_moduleOrigin
  Vector2 m_shiftPC{Vector2::Zero()};
  Transform2 m_rotationStripPC{Transform2::Identity()};  ///< Rotation to strip
  Transform2 m_translation{Transform2::Identity()};  ///< Translation to strip

  // Vectors needed for inside checking
  Vector2 m_outLeftStripPC{Vector2::Zero()};
  Vector2 m_inLeftStripPC{Vector2::Zero()};
  Vector2 m_outRightStripPC{Vector2::Zero()};
  Vector2 m_inRightStripPC{Vector2::Zero()};

  Vector2 m_outLeftModulePC{Vector2::Zero()};
  Vector2 m_inLeftModulePC{Vector2::Zero()};
  Vector2 m_outRightModulePC{Vector2::Zero()};
  Vector2 m_inRightModulePC{Vector2::Zero()};

  Vector2 m_outLeftStripXY{Vector2::Zero()};
  Vector2 m_inLeftStripXY{Vector2::Zero()};
  Vector2 m_outRightStripXY{Vector2::Zero()};
  Vector2 m_inRightStripXY{Vector2::Zero()};

  /// Pre-calculated center point (average of vertices)
  Vector2 m_center{Vector2::Zero()};

  /// Check the input values for consistency, will throw a logic_exception
  /// if consistency is not given
  void checkConsistency() noexcept(false);

  /// Transform the strip cartesian into the module polar system
  ///
  /// @param vStripXY the position in the cartesian strip system
  /// @return the position in the module polar coordinate system
  Vector2 stripXYToModulePC(const Vector2& vStripXY) const;

  Vector2 stripPCToModulePC(const Vector2& vStripPC) const;

  Vector2 modulePCToStripPC(const Vector2& vModulePC) const;

  SquareMatrix2 stripPCToModulePCJacobian(
      const Vector2& lpositionRotated) const;
};

}  // namespace Acts
