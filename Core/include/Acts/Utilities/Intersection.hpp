// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <type_traits>

#include <boost/container/static_vector.hpp>

namespace Acts {

class Surface;

/// Status enum
enum class IntersectionStatus : int {
  unreachable = 0,
  reachable = 1,
  onSurface = 2
};

/// Ostream-operator for the IntersectionStatus enum
inline std::ostream& operator<<(std::ostream& os, IntersectionStatus status) {
  constexpr static std::array<const char*, 3> names = {
      {"missed/unreachable", "reachable", "onSurface"}};

  os << names[static_cast<std::size_t>(status)];
  return os;
}

///  @struct Intersection
///
///  Intersection struct used for position
template <unsigned int DIM>
class Intersection {
 public:
  /// Position type
  using Position = Eigen::Map<const ActsVector<DIM>>;

  /// Constructor with arguments
  ///
  /// @param position is the position of the intersection
  /// @param pathLength is the path length to the intersection
  /// @param status is an enum indicating the status of the intersection
  constexpr Intersection(const ActsVector<DIM>& position, double pathLength,
                         IntersectionStatus status) noexcept
      : Intersection(std::span<const double, DIM>{position.data(), DIM},
                     pathLength, status) {}

  constexpr Intersection(const Position& position, double pathLength,
                         IntersectionStatus status) noexcept
      : Intersection(std::span<const double, DIM>{position.data(), DIM},
                     pathLength, status) {}

  constexpr Intersection(std::span<const double, DIM> position,
                         double pathLength, IntersectionStatus status) noexcept
      : m_pathLength(pathLength), m_status(status) {
    std::ranges::copy(position, m_position.begin());
  }

  Intersection(const Intersection&) noexcept = default;
  Intersection(Intersection&&) noexcept = default;
  Intersection& operator=(const Intersection&) noexcept = default;
  Intersection& operator=(Intersection&&) noexcept = default;

  /// Returns whether the intersection was successful or not
  constexpr bool isValid() const {
    return m_status != IntersectionStatus::unreachable;
  }

  /// Returns the position of the interseciton
  constexpr Position position() const { return Position{m_position.data()}; }

  /// Returns the path length to the interseciton
  constexpr double pathLength() const { return m_pathLength; }

  /// Returns the intersection status enum
  constexpr IntersectionStatus status() const { return m_status; }

  /// Static factory to creae an invalid instesection
  constexpr static Intersection invalid() { return Intersection(); }

  /// Comparison function for path length order i.e. intersection closest to
  /// -inf will be first.
  constexpr static bool pathLengthOrder(const Intersection& aIntersection,
                                        const Intersection& bIntersection) {
    auto a = aIntersection.pathLength();
    auto b = bIntersection.pathLength();
    return a < b;
  }

  /// Comparison function for closest order i.e. intersection closest to 0 will
  /// be first.
  constexpr static bool closestOrder(const Intersection& aIntersection,
                                     const Intersection& bIntersection) {
    if ((aIntersection.status() == IntersectionStatus::unreachable) &&
        (bIntersection.status() != IntersectionStatus::unreachable)) {
      return false;
    }
    if ((aIntersection.status() != IntersectionStatus::unreachable) &&
        (bIntersection.status() == IntersectionStatus::unreachable)) {
      return true;
    }
    // both are reachable or onSurface now
    auto a = aIntersection.pathLength();
    auto b = bIntersection.pathLength();
    return std::abs(a) < std::abs(b);
  }

  /// Comparison function for closest forward order i.e. intersection closest to
  /// 0 with positive path length will be first.
  constexpr static bool closestForwardOrder(const Intersection& aIntersection,
                                            const Intersection& bIntersection) {
    auto a = aIntersection.pathLength();
    auto b = bIntersection.pathLength();
    return std::signbit(a) == std::signbit(b) ? std::abs(a) < std::abs(b)
                                              : a > b;
  }

 private:
  /// Position of the intersection
  std::array<double, DIM> m_position{};
  /// Signed path length to the intersection (if valid)
  double m_pathLength = std::numeric_limits<double>::infinity();
  /// The Status of the intersection
  IntersectionStatus m_status = IntersectionStatus::unreachable;

  constexpr Intersection() = default;
};

using Intersection2D = Intersection<2>;
using Intersection3D = Intersection<3>;

static_assert(std::is_trivially_move_constructible_v<Intersection2D>);
static_assert(std::is_trivially_copy_constructible_v<Intersection2D>);
static_assert(std::is_trivially_move_assignable_v<Intersection2D>);

static constexpr std::uint8_t s_maximumNumberOfIntersections = 2;
using MultiIntersection3D =
    boost::container::static_vector<Intersection3D,
                                    s_maximumNumberOfIntersections>;

class SurfaceIntersection {
 public:
  /// @param intersection is the intersection
  /// @param surface is the surface that has been intersected
  /// @param index is the intersection index
  /// @param boundaryTolerance is the boundary tolerance for the intersection
  constexpr SurfaceIntersection(
      const Intersection3D& intersection, const Surface& surface,
      std::uint8_t index = 0,
      BoundaryTolerance boundaryTolerance = BoundaryTolerance::None()) noexcept
      : m_intersection(intersection),
        m_surface(&surface),
        m_index(index),
        m_boundaryTolerance(boundaryTolerance) {}

  SurfaceIntersection(const SurfaceIntersection&) noexcept = default;
  SurfaceIntersection(SurfaceIntersection&&) noexcept = default;
  SurfaceIntersection& operator=(const SurfaceIntersection&) noexcept = default;
  SurfaceIntersection& operator=(SurfaceIntersection&&) noexcept = default;

  /// Returns whether the intersection was successful or not
  constexpr bool isValid() const { return m_intersection.isValid(); }

  /// Returns the intersection
  constexpr const Intersection3D& intersection() const {
    return m_intersection;
  }

  /// Returns the position of the interseciton
  Intersection3D::Position position() const {
    return m_intersection.position();
  }

  /// Returns the path length to the interseciton
  constexpr double pathLength() const { return m_intersection.pathLength(); }

  /// Returns the status of the interseciton
  constexpr IntersectionStatus status() const {
    return m_intersection.status();
  }

  /// Returns the surface that has been intersected
  [[deprecated(
      "Use surface() instead. This will be removed in a future release.")]]
  constexpr const Surface* object() const {
    return m_surface;
  }

  /// Returns the surface that has been intersected
  constexpr const Surface& surface() const { return *m_surface; }

  constexpr std::uint8_t index() const { return m_index; }
  constexpr BoundaryTolerance boundaryTolerance() const {
    return m_boundaryTolerance;
  }

  constexpr static SurfaceIntersection invalid() {
    return SurfaceIntersection(Intersection3D::invalid());
  }
  constexpr static SurfaceIntersection invalid(const Surface& surface) {
    return SurfaceIntersection(Intersection3D::invalid(), surface);
  }

  constexpr static bool pathLengthOrder(
      const SurfaceIntersection& aIntersection,
      const SurfaceIntersection& bIntersection) {
    return Intersection3D::pathLengthOrder(aIntersection.intersection(),
                                           bIntersection.intersection());
  }

  constexpr static bool closestOrder(const SurfaceIntersection& aIntersection,
                                     const SurfaceIntersection& bIntersection) {
    return Intersection3D::closestOrder(aIntersection.intersection(),
                                        bIntersection.intersection());
  }

  constexpr static bool closestForwardOrder(
      const SurfaceIntersection& aIntersection,
      const SurfaceIntersection& bIntersection) {
    return Intersection3D::closestForwardOrder(aIntersection.intersection(),
                                               bIntersection.intersection());
  }

 private:
  /// The intersection itself
  Intersection3D m_intersection = Intersection3D::invalid();
  /// The surface that was (tried to be) intersected
  const Surface* m_surface = nullptr;
  /// The intersection index
  std::uint8_t m_index = 0;
  /// The boundary tolerance for the intersection
  BoundaryTolerance m_boundaryTolerance = BoundaryTolerance::None();

  constexpr SurfaceIntersection() = default;
  explicit constexpr SurfaceIntersection(const Intersection3D& intersection)
      : m_intersection(intersection) {}
};

static_assert(std::is_trivially_move_constructible_v<SurfaceIntersection>);
static_assert(std::is_trivially_move_assignable_v<SurfaceIntersection>);

class SurfaceMultiIntersection {
 public:
  using SplitIntersections =
      boost::container::static_vector<SurfaceIntersection,
                                      s_maximumNumberOfIntersections>;

  /// @param intersections are the intersections
  /// @param surface is the surface that has been intersected
  /// @param boundaryTolerance is the boundary tolerance for the intersection
  SurfaceMultiIntersection(
      const MultiIntersection3D& intersections, const Surface& surface,
      BoundaryTolerance boundaryTolerance = BoundaryTolerance::None())
      : m_intersections(intersections),
        m_surface(&surface),
        m_boundaryTolerance(boundaryTolerance) {}

  SurfaceIntersection operator[](std::uint8_t index) const {
    return {m_intersections[index], *m_surface, index, m_boundaryTolerance};
  }

  const MultiIntersection3D& intersections() const { return m_intersections; }

  std::size_t size() const { return m_intersections.size(); }

  /// Returns the surface that has been intersected
  [[deprecated(
      "Use surface() instead. This will be removed in a future release.")]]
  constexpr const Surface* object() const {
    return m_surface;
  }

  /// Returns the surface that has been intersected
  constexpr const Surface& surface() const { return *m_surface; }

  constexpr BoundaryTolerance boundaryTolerance() const {
    return m_boundaryTolerance;
  }

  SplitIntersections split() const {
    SplitIntersections result;
    for (std::size_t i = 0; i < size(); ++i) {
      result.push_back(operator[](i));
    }
    return result;
  }

  SurfaceIntersection closest() const {
    auto splitIntersections = split();
    return *std::min_element(splitIntersections.begin(),
                             splitIntersections.end(),
                             SurfaceIntersection::closestOrder);
  }

  SurfaceIntersection closestForward() const {
    auto splitIntersections = split();
    return *std::min_element(splitIntersections.begin(),
                             splitIntersections.end(),
                             SurfaceIntersection::closestForwardOrder);
  }

 private:
  /// The intersections
  MultiIntersection3D m_intersections;
  /// The surface that was (tried to be) intersected
  const Surface* m_surface = nullptr;
  /// The boundary tolerance for the intersection
  BoundaryTolerance m_boundaryTolerance = BoundaryTolerance::None();
};

namespace detail {

/// This function checks if an intersection path length is valid for the
/// specified near-limit and far-limit
///
/// @param pathLength The path length of the intersection
/// @param nearLimit The minimum path length for an intersection to be considered
/// @param farLimit The maximum path length for an intersection to be considered
/// @param logger A optionally supplied logger which prints out a lot of infos
///               at VERBOSE level
bool checkPathLength(double pathLength, double nearLimit, double farLimit,
                     const Logger& logger = getDummyLogger());

}  // namespace detail

}  // namespace Acts
