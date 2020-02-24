// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <cmath>

#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/Definitions.hpp"

#include <boost/container/small_vector.hpp>

namespace Acts {

/// @brief base class for convex polygon bounds
///
/// This class serves as a base class for the actual bounds class.
/// The only deriving type is the templated `ConvexPolygonBounds`.
///
class ConvexPolygonBoundsBase : public PlanarBounds {
 public:
  /// Output Method for std::ostream
  /// @param sl is the ostream to be written into
  std::ostream& toStream(std::ostream& sl) const final;

  /// Return vector containing defining parameters
  /// @return the parameters
  std::vector<TDD_real_t> valueStore() const final;

 protected:
  /// Return a rectangle bounds instance that encloses a set of vertices.
  /// @param vertices A collection of vertices to enclose.
  /// @return Enclosing rectangle.
  template <typename coll_t>
  static RectangleBounds makeBoundingBox(const coll_t& vertices);

  /// Calculates whether a set of vertices forms a convex polygon. This is
  /// generic over the number of vertices, so it's factored out of the concrete
  /// classes and into this base class.
  /// @param vertices A collection of vertices.
  /// @return Whether the vertices form a convex polygon.
  template <typename coll_t>
  static bool convex_impl(const coll_t& vertices);
};

/// This is the actual implementation of the bounds.
/// It is templated on the number of vertices, but there is a specialization for
/// *dynamic* number of vertices, where the underlying storage is then a vector.
///
/// @tparam N Number of vertices
template <int N>
class ConvexPolygonBounds : public ConvexPolygonBoundsBase {
  /// Expose number of vertices given as template parameter.
  ///
  static constexpr size_t num_vertices = N;
  /// Type that's used to store the vertices, in this case a fixed size array.
  ///
  using vertex_array = std::array<Vector2D, num_vertices>;

 public:
  static_assert(N >= 3, "ConvexPolygonBounds needs at least 3 sides.");

  /// Default constructor, deleted
  ConvexPolygonBounds() = delete;

  /// Constructor from a vector of vertices, to facilitate construction.
  /// This will throw if the vector size does not match `num_vertices`.
  /// This will throw if the vertices do not form a convex polygon.
  /// @param vertices The list of vertices.
  ConvexPolygonBounds(const std::vector<Vector2D>& vertices);

  /// Constructor from a fixed size array of vertices.
  /// This will throw if the vertices do not form a convex polygon.
  /// @param vertices The vertices
  ConvexPolygonBounds(const vertex_array& vertices);

  /// Defaulted destructor
  ~ConvexPolygonBounds() override = default;

  /// Return a copy of this bounds object.
  /// @return The cloned instance
  ConvexPolygonBounds<N>* clone() const final;

  /// Return the bounds type of this bounds object.
  /// @return The bounds type
  BoundsType type() const final;

  /// Return whether a local 2D point lies inside of the bounds defined by this
  /// object.
  /// @param lposition The local position to check
  /// @param bcheck The `BoundaryCheck` object handling tolerances.
  /// @return Whether the points is inside
  bool inside(const Vector2D& lposition,
              const BoundaryCheck& bcheck) const final;

  /// Return the smallest distance to any point on the boundary of this bounds
  /// object.
  /// @param lposition The local position to get the distance to
  /// @return The smallest distance to the boundary.
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

  /// Return a rectangle bounds object that encloses this polygon.
  /// @return The rectangular bounds
  const RectangleBounds& boundingBox() const final;

  /// Return whether this bounds class is in fact convex
  /// @return Whether the bounds are convex.
  bool convex() const;

 private:
  vertex_array m_vertices;
  RectangleBounds m_boundingBox;
};

/// Tag to trigger specialization of a dynamic polygon
constexpr int PolygonDynamic = -1;

/// This is the specialization handling a polygon with a dynamic number of
/// points. It can accept any number of points.
///
template <>
class ConvexPolygonBounds<PolygonDynamic> : public ConvexPolygonBoundsBase {
 public:
  /// Default constructor, deleted
  ConvexPolygonBounds() = delete;

  /// Defaulted destructor
  ~ConvexPolygonBounds() override = default;

  /// Constructor from a vector of vertices, to facilitate construction.
  /// This will throw if the vertices do not form a convex polygon.
  /// @param vertices The list of vertices.
  ConvexPolygonBounds(const std::vector<Vector2D>& vertices);

  /// Return a copy of this bounds object.
  /// @return The cloned instance
  ConvexPolygonBounds<PolygonDynamic>* clone() const final;

  /// Return the bounds type of this bounds object.
  /// @return The bounds type
  BoundsType type() const final;

  /// Return whether a local 2D point lies inside of the bounds defined by this
  /// object.
  /// @param lposition The local position to check
  /// @param bcheck The `BoundaryCheck` object handling tolerances.
  /// @return Whether the points is inside
  bool inside(const Vector2D& lposition,
              const BoundaryCheck& bcheck) const final;

  /// Return the smallest distance to any point on the boundary of this bounds
  /// object.
  /// @param lpos The lposition position to get the distance to
  /// @return The smallest distance to the boundary.
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

  ///
  /// Return a rectangle bounds object that encloses this polygon.
  /// @return The rectangular bounds
  ///
  const RectangleBounds& boundingBox() const final;

  ///
  /// Return whether this bounds class is in fact convex
  /// @return Whether the bounds are convex.
  ///
  bool convex() const;

 private:
  boost::container::small_vector<Vector2D, 10> m_vertices;
  RectangleBounds m_boundingBox;
};

}  // namespace Acts

#include "Acts/Surfaces/ConvexPolygonBounds.ipp"
