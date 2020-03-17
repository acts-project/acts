// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

/// @class RectangleBounds
///
/// Bounds for a rectangular, planar surface.
/// The two local coordinates Acts::eLOC_X, Acts::eLOC_Y are for legacy reasons
/// also called @f$ phi @f$ respectively @f$ eta @f$. The orientation
/// with respect to the local surface framce can be seen in the attached
/// illustration.
///
/// @image html RectangularBounds.gif

class RectangleBounds : public PlanarBounds {
 public:
  /// @enum BoundValues for readability
  enum BoundValues { bv_halfX = 0, bv_halfY = 1, bv_length = 2 };

  RectangleBounds() = delete;

  /// Constructor with halflength in x and y
  ///
  /// @param halex halflength in X
  /// @param haley halflength in Y
  RectangleBounds(double halex, double haley);

  /// Constructor with explicit min and max vertex
  ///
  /// @param vmin Minimum vertex
  /// @param vmax Maximum vertex
  RectangleBounds(const Vector2D& vmin, const Vector2D& vmax);

  ~RectangleBounds() override = default;

  RectangleBounds* clone() const final;

  BoundsType type() const final;

  std::vector<TDD_real_t> valueStore() const final;

  /// Inside check for the bounds object driven by the boundary check directive
  /// Each Bounds has a method inside, which checks if a LocalPosition is inside
  /// the bounds  Inside can be called without/with tolerances.
  ///
  /// @param lposition Local position (assumed to be in right surface frame)
  /// @param bcheck boundary check directive
  /// @return boolean indicator for the success of this operation
  bool inside(const Vector2D& lposition,
              const BoundaryCheck& bcheck) const final;

  /// Minimal distance to boundary ( > 0 if outside and <=0 if inside)
  ///
  /// @param lposition is the local position to check for the distance
  /// @return is a signed distance parameter
  double distanceToBoundary(const Vector2D& lposition) const final;

  /// Return the vertices
  ///
  /// @param lseg the number of segments used to approximate
  /// and eventually curved line
  ///
  /// @note the number of segements is ignored in this representation
  ///
  /// @return vector for vertices in 2D
  std::vector<Vector2D> vertices(unsigned int lseg = 1) const final;

  // Bounding box representation
  const RectangleBounds& boundingBox() const final;

  /// Output Method for std::ostream
  ///
  /// @param sl is the ostream for the dump
  std::ostream& toStream(std::ostream& sl) const final;

  /// Return method for the half length in X
  double halflengthX() const;

  /// Return method for the half length in Y
  double halflengthY() const;

  /// Get the min vertex defining the bounds
  /// @return The min vertex
  const Vector2D& min() const;

  /// Get the max vertex defining the bounds
  /// @return The max vertex
  const Vector2D& max() const;

 private:
  Vector2D m_min;
  Vector2D m_max;
};

inline double RectangleBounds::halflengthX() const {
  return std::abs(m_max.x() - m_min.x()) * 0.5;
}

inline double RectangleBounds::halflengthY() const {
  return std::abs(m_max.y() - m_min.y()) * 0.5;
}

inline SurfaceBounds::BoundsType RectangleBounds::type() const {
  return SurfaceBounds::Rectangle;
}

inline const Vector2D& RectangleBounds::min() const {
  return m_min;
}

inline const Vector2D& RectangleBounds::max() const {
  return m_max;
}

}  // namespace Acts
