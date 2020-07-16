// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

/// @brief Boolean bounds of planar surfaces
///
/// It combines two bounds on one planar surface as either
/// a union or a intersection.
///
/// As they are on the same surface, they have to have the same
/// type of local coordinate system (local cartesian), and hence
/// only a 2-dimensional shift on the local position is
/// to be applied in for inside/outside checks
///
/// Boolean operations are defined as:
/// https://en.wikipedia.org/wiki/Boolean_operations_on_polygons
class PlanarBooleanBounds : public PlanarBounds {
 public:
  PlanarBooleanBounds() = delete;

  /// Constructor with arguments
  ///
  /// @param left The first part of the union, left hand side
  /// @param leftShift The shift of the frist part wrt base surface
  /// @param right The second part of the union
  /// @param rightShift The shift of the second part wrt base surface
  /// @param bOperation The boolean operation to be applied
  PlanarBooleanBounds(std::shared_ptr<const PlanarBounds> left,
                      const Vector2D& leftShift,
                      std::shared_ptr<const PlanarBounds> right,
                      const Vector2D& rightShift,
                      BooleanOperation bOperation = eUnion) noexcept(false);

  virtual ~PlanarBooleanBounds() = default;

  /// Return the bounds type - for persistency optimization
  ///
  /// @return is a BoundsType enum
  BoundsType type() const override { return SurfaceBounds::ePlanarBoolean; }

  /// Access method for bound values, this is a dynamically sized
  /// vector containing the parameters needed to describe these bounds
  ///
  /// @note For this union it is a composite of one and the other
  ///
  /// @return of the stored values for this SurfaceBounds object
  std::vector<double> values() const override;

  /// Inside check for the bounds object driven by the boundary check directive
  ///
  /// @note For this union both components will be checked
  ///
  /// @param lposition Local position (assumed to be in right surface frame)
  /// @param bcheck boundary check directive
  ///
  /// @return boolean indicator for the success of this operation
  bool inside(const Vector2D& lposition,
              const BoundaryCheck& bcheck) const override;

  /// Minimal distance to boundary ( > 0 if outside and <=0 if inside)
  ///
  /// @note For this union the smaller distance will be provided
  ///
  /// @param lposition is the local position to check for the distance
  /// @return is a signed distance parameter
  double distanceToBoundary(const Vector2D& lposition) const override;

  /// Return the vertices
  ///
  /// @param lseg the number of segments used to approximate
  /// and eventually curved line
  ///
  /// @note that the extremas are given, which may slightly alter the
  /// number of segments returned
  ///
  /// @return vector for vertices in 2D
  std::vector<Vector2D> vertices(unsigned int lseg = 1) const override;

  /// The left side planar bounds
  ///
  /// @return a PlanarBounds reference of the boolean left side
  const PlanarBounds& leftBounds() const;

  /// The left side shift wrt base surface
  ///
  /// @return a 2-dimensional vector for left side bounds shift
  const Vector2D& leftShift() const;

  /// The right side planar bounds
  ///
  /// @return a PlanarBounds reference of the boolean right side
  const PlanarBounds& rightBounds() const;

  /// The left side shift wrt base surface
  ///
  /// @return a 2-dimensional vector for right side bounds shift
  const Vector2D& rightShift() const;

  /// The applied boolean operation
  ///
  /// @return a boolean operation indicator
  BooleanOperation bOperation() const;

  /// Indicate if there is an overlap between the bounds
  ///
  /// @return a boolean indicating if there's overlap or not
  bool boundsOverlap() const;

  /// Bounding box parameters
  ///
  /// @return rectangle bounds for a bounding box
  const RectangleBounds& boundingBox() const override;

  /// Output Method for std::ostream, to be overloaded by child classes
  ///
  /// @param sl is the outstream in which the string dump is done
  std::ostream& toStream(std::ostream& os) const override;

 private:
  std::shared_ptr<const PlanarBounds> m_left;
  Vector2D m_leftShift;

  std::shared_ptr<const PlanarBounds> m_right;
  Vector2D m_rightShift;

  BooleanOperation m_bOperation = eUnion;
  bool m_overlap = false;

  RectangleBounds m_boundingBox;
};

inline const PlanarBounds& PlanarBooleanBounds::leftBounds() const {
  return (*m_left.get());
}

inline const Vector2D& PlanarBooleanBounds::leftShift() const {
  return m_leftShift;
}

inline const PlanarBounds& PlanarBooleanBounds::rightBounds() const {
  return (*m_right.get());
}

inline const Vector2D& PlanarBooleanBounds::rightShift() const {
  return m_rightShift;
}

inline BooleanOperation PlanarBooleanBounds::bOperation() const {
  return m_bOperation;
}

inline bool PlanarBooleanBounds::boundsOverlap() const {
  return m_overlap;
}

inline const RectangleBounds& PlanarBooleanBounds::boundingBox() const {
  return m_boundingBox;
}

}  // namespace Acts