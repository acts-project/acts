// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// RadialBounds.h, c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

#pragma once
#include <cmath>

#include "Acts/Surfaces/DiscBounds.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"
#include "Acts/Utilities/VariantDataFwd.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

namespace Acts {

/// @class RadialBounds
///
/// Class to describe the bounds for a planar DiscSurface.
/// By providing an argument for hphisec, the bounds can
/// be restricted to a phi-range around the center position.
///
/// @image html RadialBounds.gif

class RadialBounds : public DiscBounds
{
public:
  /// enumeration for readability
  enum BoundValues {
    bv_rMin          = 0,
    bv_rMax          = 1,
    bv_averagePhi    = 2,
    bv_halfPhiSector = 3,
    bv_length        = 4
  };

  RadialBounds();

  /// Constructor for full disc of symmetric disc around phi=0
  ///
  /// @param minrad is the inner radius of the disc (0 for full disc)
  /// @param maxrad is the outer radius of the disc
  /// @param hphisec is the half opening angle of the disc (Pi for full angular
  /// coverage)
  RadialBounds(double minrad, double maxrad, double hphisec = M_PI);

  /// Constructor for full disc of symmetric disc around phi!=0
  ///
  /// @param minrad is the inner radius of the disc (0 for full disc)
  /// @param maxrad is the outer radius of the disc
  /// @param avephi is the phi value of the local x-axis in the local 3D frame
  /// @param hphisec is the half opening angle of the disc (Pi for full angular
  /// coverage)
  RadialBounds(double minrad, double maxrad, double avephi, double hphisec);

  /// Constructor which accepts @c variant_data
  ///
  /// @param vardata the @c variant_data to build from
  RadialBounds(const variant_data& vardata);

  ~RadialBounds() override;

  RadialBounds*
  clone() const final;

  SurfaceBounds::BoundsType
  type() const final;

  std::vector<TDD_real_t>
  valueStore() const final;

  /// For disc surfaces the local position in (r,phi) is checked
  ///
  /// @param lpos local position to be checked
  /// @param bcheck boundary check directive
  ///
  /// @return is a boolean indicating the operation success
  bool
  inside(const Vector2D& lpos, const BoundaryCheck& bcheck) const final;

  /// Minimal distance to boundary calculation
  ///
  /// @param lpos local 2D position in surface coordinate frame
  ///
  /// @return distance to boundary ( > 0 if outside and <=0 if inside)
  double
  distanceToBoundary(const Vector2D& lpos) const final;

  /// Outstream operator
  ///
  /// @param sl is the ostream to be dumped into
  std::ostream&
  dump(std::ostream& sl) const final;

  /// Return method for inner Radius
  double
  rMin() const;

  /// Return method for outer Radius
  double
  rMax() const;

  /// Return method for the central phi value
  ///(i.e. phi value of x-axis of local 3D frame)
  double
  averagePhi() const;

  /// Return method for the half phi span
  double
  halfPhiSector() const;

  /// Produce a @c variant_data representation of this object
  /// @return The representation
  variant_data
  toVariantData() const override;

private:
  double m_rMin, m_rMax, m_avgPhi, m_halfPhi;

  Vector2D
  shifted(const Vector2D& lpos) const;
};

inline double
RadialBounds::rMin() const
{
  return m_rMin;
}

inline double
RadialBounds::rMax() const
{
  return m_rMax;
}

inline double
RadialBounds::averagePhi() const
{
  return m_avgPhi;
}

inline double
RadialBounds::halfPhiSector() const
{
  return m_halfPhi;
}

}  // namespace