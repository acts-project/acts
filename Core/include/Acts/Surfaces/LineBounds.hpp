// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// LineBounds.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/VariantDataFwd.hpp"

namespace Acts {

/// @class LineBounds
///
/// Bounds for a LineSurface.
///

class LineBounds : public SurfaceBounds
{
public:
  /// @enum BoundValues for readablility
  /// nested enumeration object
  enum BoundValues { bv_radius = 0, bv_halfZ = 1, bv_length = 2 };

  /// Constructor
  ///
  /// @param radius is the radius of the cylinder, default = 0.
  /// @param halez is the half length in z, defualt = 0.
  LineBounds(double radius = 0., double halez = 0.);

  /// Constructor which accepts @c variant_data
  ///
  /// @param vardata the @c variant_data to build from
  LineBounds(const variant_data& vardata);

  ~LineBounds() override;

  LineBounds*
  clone() const final;

  BoundsType
  type() const final;

  std::vector<TDD_real_t>
  valueStore() const final;

  /// Inside check for the bounds object driven by the boundary check directive
  /// Each Bounds has a method inside, which checks if a LocalPosition is inside
  /// the bounds  Inside can be called without/with tolerances.
  ///
  /// @param lpos Local position (assumed to be in right surface frame)
  /// @param bcheck boundary check directive
  ///
  /// @return boolean indicator for the success of this operation
  bool
  inside(const Vector2D& lpos, const BoundaryCheck& bcheck) const final;

  /// Minimal distance to boundary ( > 0 if outside and <=0 if inside)
  ///
  /// @param lpos is the local position to check for the distance
  ///
  /// @return is a signed distance parameter
  double
  distanceToBoundary(const Vector2D& lpos) const final;

  /// This method returns the radius
  virtual double
  r() const;

  /// This method returns the halflengthZ
  double
  halflengthZ() const;

  /// Output Method for std::ostream
  ///
  /// @param sl is the ostream to be dumped into
  std::ostream&
  dump(std::ostream& sl) const final;

  /// Produce a @c variant_data representation of this object
  /// @return The representation
  variant_data
  toVariantData() const override;

private:
  double m_radius, m_halfZ;
};

inline double
LineBounds::r() const
{
  return m_radius;
}

inline double
LineBounds::halflengthZ() const
{
  return m_halfZ;
}

}  // namespace