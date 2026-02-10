// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"

#include <array>
#include <cstddef>
#include <iosfwd>
#include <span>
#include <vector>

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
  /// @return Reference to the output stream after writing
  std::ostream& toStream(std::ostream& sl) const final;

  /// Return the bounds type of this bounds object.
  /// @return The bounds type
  Type type() const final { return ConvexPolygon; }

  /// Return the bound values as dynamically sized vector
  /// @return this returns a copy of the internal values
  std::vector<double> values() const final;

  /// @copydoc SurfaceBounds::center
  /// @note For ConvexPolygonBounds: returns average of all vertices (vertex
  ///       centroid)
  Vector2 center() const final;

  /// Return a rectangle bounds object that encloses this polygon.
  /// @return The rectangular bounds
  ///
  const RectangleBounds& boundingBox() const final;

 protected:
  /// Creates a rectangle bounds instance that encloses a set of vertices.
  /// @param vertices A collection of vertices to enclose.
  void makeBoundingBox(std::span<const Vector2> vertices);

  /// Calculates whether a set of vertices forms a convex polygon. This is
  /// generic over the number of vertices, so it's factored out of the concrete
  /// classes and into this base class.
  /// @param vertices A collection of vertices.
  /// throws a logic error if this is not the case
  template <typename coll_t>
    requires std::same_as<typename coll_t::value_type, Acts::Vector2>
  static void convex_impl(const coll_t& vertices) noexcept(false);

  /// Calculate and cache the center point from vertices
  /// @param vertices The vertices to calculate center from
  void calculateCenter(std::span<const Vector2> vertices);

  /// Return whether this bounds class is in fact convex
  /// thorws a logic error if not
  /// @param vertices The vertices to check for consistency
  static void checkConsistency(std::span<const Vector2> vertices) noexcept(
      false);

 private:
  /// Cached center position
  Vector2 m_center{Vector2::Zero()};
  RectangleBounds m_boundingBox{0, 0};
};

template <int N>
concept isValidConvexPolygonSize = requires { requires(N >= 3 || N == -1); };

/// This is the actual implementation of the bounds.
/// It is templated on the number of vertices, but there is a specialization
/// for *dynamic* number of vertices, where the underlying storage is then a
/// vector.
///
/// @tparam N Number of vertices
template <int N>
  requires isValidConvexPolygonSize<N>
class ConvexPolygonBounds : public ConvexPolygonBoundsBase {
 public:
  /// Expose number of vertices given as template parameter.
  static constexpr std::size_t nVertices = N;

  /// Expose number of parameters as a template parameter
  /// @note The `eSize` name here emulates the size of the *bound values* in other
  ///       bounds classes.
  static constexpr std::size_t eSize = 2 * N;

  /// Constructor from a vector of vertices, to facilitate construction.
  /// This will throw if the vector size does not match `num_vertices`.
  /// This will throw if the vertices do not form a convex polygon.
  /// @param vertices The list of vertices.
  explicit ConvexPolygonBounds(std::span<const Vector2> vertices) noexcept(
      false);

  /// Constructor from a fixed size array of parameters
  /// This will throw if the vertices do not form a convex polygon.
  /// @param values The values to build up the vertices
  explicit ConvexPolygonBounds(std::span<const double> values) noexcept(false);

  /// @copydoc SurfaceBounds::inside
  bool inside(const Vector2& lposition) const final;

  /// @copydoc SurfaceBounds::closestPoint
  Vector2 closestPoint(const Vector2& lposition,
                       const SquareMatrix2& metric) const final;

  using SurfaceBounds::inside;

  /// Return the vertices
  ///
  /// @param ignoredSegments the number of segments used to approximate
  /// and eventually curved line
  ///
  /// @note the number of segments is ignored in this representation
  ///
  /// @return vector for vertices in 2D
  std::vector<Vector2> vertices(unsigned int ignoredSegments = 0u) const final;

 private:
  std::array<Vector2, nVertices> m_vertices{};
};

/// Tag to trigger specialization of a dynamic polygon
constexpr int PolygonDynamic = -1;

/// This is the specialization handling a polygon with a dynamic number of
/// points. It can accept any number of points.
///
template <>
class ConvexPolygonBounds<PolygonDynamic> : public ConvexPolygonBoundsBase {
 public:
  /// Expose number of vertices given as template parameter.
  constexpr static int nVertices = PolygonDynamic;
  /// Expose number of parameters as a template parameter
  /// @note The `eSize` name here emulates the size of the *bound values* in other
  ///       bounds classes.
  constexpr static int eSize = -1;

  /// Constructor from a vector of vertices, to facilitate construction.
  /// This will throw if the vertices do not form a convex polygon.
  /// @param vertices The list of vertices.
  explicit ConvexPolygonBounds(std::span<const Vector2> vertices);

  /// Constructor from a vector of vertices, to facilitate construction.
  /// This will throw if the vertices do not form a convex polygon.
  /// @param vertices The list of vertices.
  explicit ConvexPolygonBounds(const std::vector<Vector2>& vertices);

  /// @copydoc SurfaceBounds::inside
  bool inside(const Vector2& lposition) const final;

  /// @copydoc SurfaceBounds::closestPoint
  Vector2 closestPoint(const Vector2& lposition,
                       const SquareMatrix2& metric) const final;

  using SurfaceBounds::inside;

  /// Return the vertices
  ///
  /// @param lseg the number of segments used to approximate
  /// and eventually curved line
  ///
  /// @note the number of segments is ignored in this representation
  ///
  /// @return vector for vertices in 2D
  std::vector<Vector2> vertices(unsigned int lseg = 1) const final;

 private:
  boost::container::small_vector<Vector2, 10> m_vertices{};
};

}  // namespace Acts

#include "Acts/Surfaces/ConvexPolygonBounds.ipp"
