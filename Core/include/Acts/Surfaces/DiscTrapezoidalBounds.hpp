// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// DiscTrapezoidalBounds.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#include <cmath>

#include "Acts/Surfaces/DiscBounds.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"
#include "Acts/Utilities/VariantDataFwd.hpp"

namespace Acts {

///
/// @class DiscTrapezoidalBounds
///
/// Class to describe the bounds for a planar DiscSurface.
/// By providing an argument for hphisec, the bounds can
/// be restricted to a phi-range around the center position.
///

class DiscTrapezoidalBounds : public DiscBounds
{
public:
  /// @enum BoundValues
  /// enumeration for readability
  enum BoundValues {
    bv_rMin       = 0,
    bv_rMax       = 1,
    bv_minHalfX   = 2,
    bv_maxHalfX   = 3,
    bv_averagePhi = 4,
    bv_stereo     = 5,
    bv_length     = 6
  };

  DiscTrapezoidalBounds() = delete;

  /// Constructor for a symmetric Trapezoid giving min X length, max X length,
  /// Rmin and R max
  /// @param minhalfx half length in X at min radius
  /// @param maxhalfx half length in X at maximum radius
  /// @param maxR outer radius
  /// @param minR inner radius
  /// @param avephi average phi value
  /// @param stereo optional stero angle applied
  DiscTrapezoidalBounds(double minhalfx,
                        double maxhalfx,
                        double maxR,
                        double minR,
                        double avephi = M_PI_2,
                        double stereo = 0.);

  /// Constructor which accepts @c variant_data
  ///
  /// @param vardata the @c variant_data to build from
  DiscTrapezoidalBounds(const variant_data& vardata);

  ~DiscTrapezoidalBounds() override;

  DiscTrapezoidalBounds*
  clone() const final;

  SurfaceBounds::BoundsType
  type() const final;

  std::vector<TDD_real_t>
  valueStore() const final;

  ///  This method cheks if the radius given in the LocalPosition is inside
  ///  [rMin,rMax]
  /// if only tol0 is given and additional in the phi sector is tol1 is given
  /// @param lpos is the local position to be checked (in polar coordinates)
  /// @param bcheck is the boundary check directive
  bool
  inside(const Vector2D& lpos, const BoundaryCheck& bcheck = true) const final;

  /// Minimal distance to boundary
  /// @param lpos is the local position to be checked (in polar coordinates)
  /// @return is the minimal distance ( > 0 if outside and <=0 if inside)
  double
  distanceToBoundary(const Vector2D& lpos) const final;

  /// Output Method for std::ostream
  std::ostream&
  dump(std::ostream& sl) const final;

  /// This method returns inner radius
  double
  rMin() const;

  /// This method returns outer radius
  double
  rMax() const;

  /// This method returns the average phi
  double
  averagePhi() const;

  /// This method returns the center radius
  double
  rCenter() const;

  /// This method returns the stereo angle
  double
  stereo() const;

  /// This method returns the halfPhiSector which is covered by the disc
  double
  halfPhiSector() const;

  /// This method returns the minimal halflength in X
  double
  minHalflengthX() const;

  /// This method returns the maximal halflength in X
  double
  maxHalflengthX() const;

  /// This method returns the halflength in Y (this is Rmax -Rmin)
  double
  halflengthY() const;

  /// Produce a @c variant_data representation of this object
  /// @return The representation
  variant_data
  toVariantData() const override;

private:
  double m_rMin, m_rMax, m_minHalfX, m_maxHalfX, m_avgPhi;
  double m_stereo;  // TODO 2017-04-09 msmk: what is this good for?

  Vector2D
  toLocalXY(const Vector2D& lpos) const;
  ActsMatrixD<2, 2>
  jacobianToLocalXY(const Vector2D& lpos) const;
};

inline double
DiscTrapezoidalBounds::rMin() const
{
  return m_rMin;
}

inline double
DiscTrapezoidalBounds::rMax() const
{
  return m_rMax;
}

inline double
DiscTrapezoidalBounds::minHalflengthX() const
{
  return m_minHalfX;
}

inline double
DiscTrapezoidalBounds::maxHalflengthX() const
{
  return m_maxHalfX;
}

inline double
DiscTrapezoidalBounds::averagePhi() const
{
  return m_avgPhi;
}

inline double
DiscTrapezoidalBounds::stereo() const
{
  return m_stereo;
}

inline double
DiscTrapezoidalBounds::halfPhiSector() const
{
  auto minHalfPhi = std::asin(m_minHalfX / m_rMin);
  auto maxHalfPhi = std::asin(m_maxHalfX / m_rMax);
  return std::max(minHalfPhi, maxHalfPhi);
}

inline double
DiscTrapezoidalBounds::rCenter() const
{
  auto hmin = std::sqrt(m_rMin * m_rMin - m_minHalfX * m_minHalfX);
  auto hmax = std::sqrt(m_rMax * m_rMax - m_maxHalfX * m_maxHalfX);
  return (hmin + hmax) / 2.0;
}

inline double
DiscTrapezoidalBounds::halflengthY() const
{
  auto hmin = std::sqrt(m_rMin * m_rMin - m_minHalfX * m_minHalfX);
  auto hmax = std::sqrt(m_rMax * m_rMax - m_maxHalfX * m_maxHalfX);
  return (hmax - hmin) / 2.0;
}

}  // namespace Acts