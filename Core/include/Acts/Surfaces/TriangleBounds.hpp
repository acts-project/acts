// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// TriangleBounds.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#include <array>
#include <utility>

#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

namespace Acts {

/// @class TriangleBounds
///
/// Bounds for a triangular, planar surface.
///
/// @image html TriangularBounds.gif

class TriangleBounds : public PlanarBounds {
 public:
  // @enum BoundValues for readability
  enum BoundValues {
    bv_x1 = 0,
    bv_y1 = 1,
    bv_x2 = 2,
    bv_y2 = 3,
    bv_x3 = 4,
    bv_y3 = 5,
    bv_length = 6
  };

  TriangleBounds() = delete;

  /// Constructor with coordinates of vertices
  ///
  /// @param vertices is the vector of vertices
  TriangleBounds(const std::array<Vector2D, 3>& vertices);

  ~TriangleBounds() override;

  TriangleBounds* clone() const final;

  BoundsType type() const final;

  std::vector<TDD_real_t> valueStore() const final;

  /// This method checks if the provided local coordinates are inside the
  /// surface bounds
  ///
  /// @param lposition local position in 2D local Cartesian frame
  /// @param bcheck is the boundary check directive
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
  /// @param sl is the ostream to be dumped into
  std::ostream& toStream(std::ostream& sl) const final;

 private:
  std::array<Vector2D, 3> m_vertices;
  RectangleBounds m_boundingBox;  ///< internal bounding box cache
};

}  // namespace Acts