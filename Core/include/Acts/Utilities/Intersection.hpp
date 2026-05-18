// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <type_traits>

namespace Acts {

/// Status enum
enum class IntersectionStatus : int {
  unreachable = 0,
  reachable = 1,
  onSurface = 2
};

/// Ostream-operator for the IntersectionStatus enum
/// @param os Output stream
/// @param status IntersectionStatus to output
/// @return Reference to output stream
inline std::ostream& operator<<(std::ostream& os, IntersectionStatus status) {
  constexpr static std::array<const char*, 3> names = {
      {"missed/unreachable", "reachable", "onSurface"}};

  os << names[static_cast<std::size_t>(status)];
  return os;
}

/// Intersection struct containing the position, path length and status of an
/// intersection.
template <unsigned int DIM>
class Intersection {
 public:
  /// Position type
  using Position = Eigen::Map<const Vector<DIM>>;

  /// Constructor with arguments
  ///
  /// @param position is the position of the intersection
  /// @param pathLength is the path length to the intersection
  /// @param status is an enum indicating the status of the intersection
  constexpr Intersection(const Vector<DIM>& position, double pathLength,
                         IntersectionStatus status) noexcept
      : Intersection(std::span<const double, DIM>{position.data(), DIM},
                     pathLength, status) {}

  /// Constructor from position vector, path length, and status
  /// @param position The intersection position
  /// @param pathLength The path length to the intersection
  /// @param status The intersection status
  constexpr Intersection(const Position& position, double pathLength,
                         IntersectionStatus status) noexcept
      : Intersection(std::span<const double, DIM>{position.data(), DIM},
                     pathLength, status) {}

  /// Constructor from position span, path length, and status
  /// @param position Span of position coordinates
  /// @param pathLength The path length to the intersection
  /// @param status The intersection status
  constexpr Intersection(std::span<const double, DIM> position,
                         double pathLength, IntersectionStatus status) noexcept
      : m_pathLength(pathLength), m_status(status) {
    std::ranges::copy(position, m_position.begin());
  }

  /// Copy constructor
  constexpr Intersection(const Intersection&) noexcept = default;
  /// Move constructor
  constexpr Intersection(Intersection&&) noexcept = default;
  /// Copy assignment operator
  /// @return Reference to this intersection for chaining
  constexpr Intersection& operator=(const Intersection&) noexcept = default;
  /// Move assignment operator
  /// @return Reference to this intersection for chaining
  constexpr Intersection& operator=(Intersection&&) noexcept = default;

  /// Returns whether the intersection was successful or not
  /// @return True if intersection is reachable or on surface, false if unreachable
  constexpr bool isValid() const noexcept {
    return m_status != IntersectionStatus::unreachable;
  }

  /// Returns the position of the interseciton
  /// @return Position vector of the intersection point
  Position position() const noexcept { return Position{m_position.data()}; }

  /// Returns the path length to the intersection
  /// @return Signed path length from origin to intersection point
  constexpr double pathLength() const noexcept { return m_pathLength; }

  /// Returns the intersection status enum
  /// @return Status indicating if intersection is unreachable, reachable, or on surface
  constexpr IntersectionStatus status() const noexcept { return m_status; }

  /// Static factory to create an invalid intersection
  /// @return Invalid intersection with unreachable status
  constexpr static Intersection Invalid() noexcept { return Intersection(); }

  /// Comparison function for path length order i.e. intersection closest to
  /// -inf will be first.
  /// @param aIntersection First intersection to compare
  /// @param bIntersection Second intersection to compare
  /// @return True if first intersection has smaller path length than second
  constexpr static bool pathLengthOrder(
      const Intersection& aIntersection,
      const Intersection& bIntersection) noexcept {
    auto a = aIntersection.pathLength();
    auto b = bIntersection.pathLength();
    return a < b;
  }

  /// Comparison function for closest order i.e. intersection closest to 0 will
  /// be first.
  /// @param aIntersection First intersection to compare
  /// @param bIntersection Second intersection to compare
  /// @return True if first intersection is closer to zero path length than second
  constexpr static bool closestOrder(
      const Intersection& aIntersection,
      const Intersection& bIntersection) noexcept {
    using enum IntersectionStatus;

    if ((aIntersection.status() == unreachable) &&
        (bIntersection.status() != unreachable)) {
      return false;
    }
    if ((aIntersection.status() != unreachable) &&
        (bIntersection.status() == unreachable)) {
      return true;
    }
    // both are reachable or onSurface now
    auto a = aIntersection.pathLength();
    auto b = bIntersection.pathLength();
    return std::abs(a) < std::abs(b);
  }

  /// Comparison function for closest forward order i.e. intersection closest to
  /// 0 with positive path length will be first.
  /// @param aIntersection First intersection to compare
  /// @param bIntersection Second intersection to compare
  /// @return True if first intersection is closer to zero with preference for forward direction
  constexpr static bool closestForwardOrder(
      const Intersection& aIntersection,
      const Intersection& bIntersection) noexcept {
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

  constexpr Intersection() noexcept = default;
};

/// Type alias for 2D intersection
using Intersection2D = Intersection<2>;
/// Type alias for 3D intersection
using Intersection3D = Intersection<3>;

static_assert(std::is_trivially_copy_constructible_v<Intersection2D>);
static_assert(std::is_trivially_move_constructible_v<Intersection2D>);
static_assert(std::is_trivially_move_assignable_v<Intersection2D>);

/// Index type for intersections
using IntersectionIndex = std::uint8_t;
/// Maximum number of intersections that can be stored
static constexpr IntersectionIndex s_maximumNumberOfIntersections = 2;

/// Container for up to two intersections in a given dimension.
template <unsigned int DIM>
class MultiIntersection {
 public:
  /// Intersection type for this dimension
  using IntersectionType = Intersection<DIM>;
  /// Pair of intersection and its index
  using IndexedIntersection = std::pair<IntersectionType, IntersectionIndex>;

  /// Container type for storing intersections
  using Container =
      std::array<IntersectionType, s_maximumNumberOfIntersections>;

  /// Size type for indexing
  using size_type = IntersectionIndex;

  /// Construct from single intersection
  /// @param intersection The intersection
  constexpr explicit MultiIntersection(
      const IntersectionType& intersection) noexcept
      : m_intersections{intersection, IntersectionType::Invalid()}, m_size{1} {}
  /// Construct from two intersections
  /// @param intersection1 The first intersection
  /// @param intersection2 The second intersection
  constexpr MultiIntersection(const IntersectionType& intersection1,
                              const IntersectionType& intersection2) noexcept
      : m_intersections{intersection1, intersection2}, m_size{2} {}

  /// Copy constructor
  constexpr MultiIntersection(const MultiIntersection&) noexcept = default;
  /// Move constructor
  constexpr MultiIntersection(MultiIntersection&&) noexcept = default;
  /// Copy assignment operator
  /// @return Reference to this object
  constexpr MultiIntersection& operator=(const MultiIntersection&) noexcept =
      default;
  /// Move assignment operator
  /// @return Reference to this object
  constexpr MultiIntersection& operator=(MultiIntersection&&) noexcept =
      default;

  /// Access intersection by index
  /// @param index The index of the intersection
  /// @return Reference to the intersection
  constexpr const IntersectionType& operator[](IntersectionIndex index) const {
    return m_intersections[index];
  }

  /// Access intersection at index with bounds checking
  /// @param index The index of the intersection
  /// @return Reference to the intersection
  constexpr const IntersectionType& at(IntersectionIndex index) const {
    return m_intersections.at(index);
  }

  /// Get the number of intersections
  /// @return The number of intersections
  constexpr IntersectionIndex size() const noexcept { return m_size; }

  /// Get begin iterator
  /// @return Iterator to the beginning
  constexpr auto begin() const noexcept {
    return std::span(m_intersections.data(), m_size).begin();
  }
  /// Get end iterator
  /// @return Iterator to the end
  constexpr auto end() const noexcept {
    return std::span(m_intersections.data(), m_size).end();
  }

  /// Get closest intersection
  /// @return The closest intersection
  constexpr IntersectionType closest() const noexcept {
    return closestWithIndex().first;
  }
  /// Get closest intersection with its index
  /// @return Pair of intersection and its index
  constexpr IndexedIntersection closestWithIndex() const noexcept {
    auto min = std::ranges::min_element(m_intersections,
                                        IntersectionType::closestOrder);
    return {*min, static_cast<IntersectionIndex>(
                      std::distance(m_intersections.begin(), min))};
  }

  /// Get closest forward intersection
  /// @return The closest forward intersection
  constexpr IntersectionType closestForward() const noexcept {
    return closestForwardWithIndex().first;
  }
  /// Get closest forward intersection with its index
  /// @return Pair of intersection and its index
  constexpr IndexedIntersection closestForwardWithIndex() const noexcept {
    auto min = std::ranges::min_element(m_intersections,
                                        IntersectionType::closestForwardOrder);
    return {*min, static_cast<IntersectionIndex>(
                      std::distance(m_intersections.begin(), min))};
  }

 private:
  Container m_intersections{};
  IntersectionIndex m_size{};
};

/// Container for up to two 2D intersections
using MultiIntersection2D = MultiIntersection<2>;
/// Container for up to two 3D intersections
using MultiIntersection3D = MultiIntersection<3>;

static_assert(std::is_trivially_copy_constructible_v<MultiIntersection2D>);
static_assert(std::is_trivially_move_constructible_v<MultiIntersection2D>);
static_assert(std::is_trivially_move_assignable_v<MultiIntersection2D>);

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
