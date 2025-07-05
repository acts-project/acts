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

#include <boost/container/static_vector.hpp>

namespace Acts {

class Surface;

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

/// Intersection struct containing the position, path length and status of an
/// intersection.
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

  constexpr Intersection(const Intersection&) noexcept = default;
  constexpr Intersection(Intersection&&) noexcept = default;
  constexpr Intersection& operator=(const Intersection&) noexcept = default;
  constexpr Intersection& operator=(Intersection&&) noexcept = default;

  /// Returns whether the intersection was successful or not
  constexpr bool isValid() const noexcept {
    return m_status != IntersectionStatus::unreachable;
  }

  /// Returns the position of the interseciton
  constexpr Position position() const noexcept {
    return Position{m_position.data()};
  }

  /// Returns the path length to the interseciton
  constexpr double pathLength() const noexcept { return m_pathLength; }

  /// Returns the intersection status enum
  constexpr IntersectionStatus status() const noexcept { return m_status; }

  /// Static factory to creae an invalid instesection
  constexpr static Intersection invalid() noexcept { return Intersection(); }

  /// Comparison function for path length order i.e. intersection closest to
  /// -inf will be first.
  constexpr static bool pathLengthOrder(
      const Intersection& aIntersection,
      const Intersection& bIntersection) noexcept {
    auto a = aIntersection.pathLength();
    auto b = bIntersection.pathLength();
    return a < b;
  }

  /// Comparison function for closest order i.e. intersection closest to 0 will
  /// be first.
  constexpr static bool closestOrder(
      const Intersection& aIntersection,
      const Intersection& bIntersection) noexcept {
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

using Intersection2D = Intersection<2>;
using Intersection3D = Intersection<3>;

static_assert(std::is_trivially_move_constructible_v<Intersection2D>);
static_assert(std::is_trivially_copy_constructible_v<Intersection2D>);
static_assert(std::is_trivially_move_assignable_v<Intersection2D>);

using IntersectionIndex = std::uint8_t;

template <unsigned int DIM>
class IndexedIntersection {
 public:
  using Intersection = Intersection<DIM>;
  using Position = Intersection::Position;

  constexpr IndexedIntersection(const Intersection& intersection,
                                IntersectionIndex index) noexcept
      : m_intersection(intersection), m_index(index) {}

  constexpr IndexedIntersection(const IndexedIntersection&) noexcept = default;
  constexpr IndexedIntersection(IndexedIntersection&&) noexcept = default;
  constexpr IndexedIntersection& operator=(
      const IndexedIntersection&) noexcept = default;
  constexpr IndexedIntersection& operator=(IndexedIntersection&&) noexcept =
      default;

  constexpr const Intersection& intersection() const { return m_intersection; }

  constexpr IntersectionIndex index() const { return m_index; }

  /// Returns whether the intersection was successful or not
  constexpr bool isValid() const { return m_intersection.isValid(); }

  /// Returns the position of the interseciton
  constexpr Position position() const { return m_intersection.position(); }

  /// Returns the path length to the interseciton
  constexpr double pathLength() const { return m_intersection.pathLength(); }

  /// Returns the intersection status enum
  constexpr IntersectionStatus status() const {
    return m_intersection.status();
  }

 private:
  Intersection m_intersection = Intersection::invalid();
  IntersectionIndex m_index = 0;
};

using IndexedIntersection2D = IndexedIntersection<2>;
using IndexedIntersection3D = IndexedIntersection<3>;

static constexpr IntersectionIndex s_maximumNumberOfIntersections = 2;

template <unsigned int DIM>
class MultiIntersection {
 public:
  using Intersection = Intersection<DIM>;

  using Container = std::array<Intersection, s_maximumNumberOfIntersections>;

  constexpr explicit MultiIntersection(
      const Intersection& intersection) noexcept
      : m_intersections{intersection, Intersection::invalid()}, m_size{1} {}
  constexpr MultiIntersection(const Intersection& intersection1,
                              const Intersection& intersection2) noexcept
      : m_intersections{intersection1, intersection2}, m_size{2} {}

  constexpr MultiIntersection(const MultiIntersection&) noexcept = default;
  constexpr MultiIntersection(MultiIntersection&&) noexcept = default;
  constexpr MultiIntersection& operator=(const MultiIntersection&) noexcept =
      default;
  constexpr MultiIntersection& operator=(MultiIntersection&&) noexcept =
      default;

  constexpr Intersection& at(IntersectionIndex index) {
    return m_intersections.at(index);
  }
  constexpr const Intersection& at(IntersectionIndex index) const {
    return m_intersections.at(index);
  }

  constexpr std::uint8_t size() const noexcept { return m_size; }

  class Iterator {
   public:
    using container_iterator = Container::const_iterator;

    using value_type = IndexedIntersection<DIM>;
    using difference_type = std::ptrdiff_t;
    using pointer = const value_type*;
    using reference = const value_type&;
    using iterator_category = std::forward_iterator_tag;

    constexpr Iterator(container_iterator it, IntersectionIndex index) noexcept
        : m_it(it), m_index(index) {}

    constexpr Iterator& operator++() noexcept {
      ++m_it;
      ++m_index;
      return *this;
    }

    constexpr value_type operator*() const noexcept { return {*m_it, m_index}; }

   private:
    container_iterator m_it;
    IntersectionIndex m_index;

    friend bool operator==(const Iterator& lhs, const Iterator& rhs) noexcept {
      return lhs.m_it == rhs.m_it;
    }
  };

  constexpr auto begin() const noexcept {
    return Iterator(m_intersections.begin(), 0);
  }
  constexpr auto end() const noexcept {
    return Iterator(m_intersections.begin() + m_size, m_size);
  }

  constexpr IndexedIntersection<DIM> closest() const noexcept {
    auto min =
        std::ranges::min_element(m_intersections, Intersection::closestOrder);
    return {*min, static_cast<IntersectionIndex>(
                      std::distance(m_intersections.begin(), min))};
  }

  constexpr IndexedIntersection<DIM> closestForward() const noexcept {
    auto min = std::ranges::min_element(m_intersections,
                                        Intersection::closestForwardOrder);
    return {*min, static_cast<IntersectionIndex>(
                      std::distance(m_intersections.begin(), min))};
  }

 private:
  Container m_intersections;
  std::uint8_t m_size = 0;
};

using MultiIntersection2D = MultiIntersection<2>;
using MultiIntersection3D = MultiIntersection<3>;

class SurfaceIntersection {
 public:
  using Intersection = Intersection3D;
  using Position = Intersection::Position;

  /// Create a surface intersection from a 3D intersection, intersection index,
  /// and a surface
  ///
  /// @param intersection is the intersection
  /// @param index is the intersection index
  /// @param object is the object to be instersected
  constexpr SurfaceIntersection(const Intersection3D& intersection,
                                IntersectionIndex index,
                                const Surface& surface) noexcept
      : m_intersection(intersection), m_index(index), m_surface(&surface) {}
  constexpr SurfaceIntersection(
      const IndexedIntersection3D& indexedIntersection,
      const Surface& surface) noexcept
      : m_intersection(indexedIntersection.intersection()),
        m_index(indexedIntersection.index()),
        m_surface(&surface) {}

  constexpr SurfaceIntersection(const SurfaceIntersection&) noexcept = default;
  constexpr SurfaceIntersection(SurfaceIntersection&&) noexcept = default;
  constexpr SurfaceIntersection& operator=(
      const SurfaceIntersection&) noexcept = default;
  constexpr SurfaceIntersection& operator=(SurfaceIntersection&&) noexcept =
      default;

  /// Returns the intersection
  constexpr const Intersection3D& intersection() const {
    return m_intersection;
  }

  constexpr IntersectionIndex index() const noexcept { return m_index; }

  /// Returns the surface that has been intersected
  constexpr const Surface& surface() const noexcept { return *m_surface; }

  /// Returns whether the intersection was successful or not
  constexpr bool isValid() const { return m_intersection.isValid(); }

  /// Returns the position of the interseciton
  Position position() const noexcept { return m_intersection.position(); }

  /// Returns the path length to the interseciton
  constexpr double pathLength() const noexcept {
    return m_intersection.pathLength();
  }

  /// Returns the status of the interseciton
  constexpr IntersectionStatus status() const noexcept {
    return m_intersection.status();
  }

  constexpr static SurfaceIntersection invalid() noexcept {
    return SurfaceIntersection();
  }
  constexpr static SurfaceIntersection invalid(
      const Surface& surface) noexcept {
    return SurfaceIntersection(Intersection3D::invalid(), 0, surface);
  }

  constexpr static bool pathLengthOrder(
      const SurfaceIntersection& aIntersection,
      const SurfaceIntersection& bIntersection) noexcept {
    return Intersection3D::pathLengthOrder(aIntersection.intersection(),
                                           bIntersection.intersection());
  }

  constexpr static bool closestOrder(
      const SurfaceIntersection& aIntersection,
      const SurfaceIntersection& bIntersection) noexcept {
    return Intersection3D::closestOrder(aIntersection.intersection(),
                                        bIntersection.intersection());
  }

  constexpr static bool closestForwardOrder(
      const SurfaceIntersection& aIntersection,
      const SurfaceIntersection& bIntersection) noexcept {
    return Intersection3D::closestForwardOrder(aIntersection.intersection(),
                                               bIntersection.intersection());
  }

 private:
  /// The intersection itself
  Intersection3D m_intersection = Intersection3D::invalid();
  /// The intersection index
  IntersectionIndex m_index = 0;
  /// The surface that was (tried to be) intersected
  const Surface* m_surface = nullptr;

  constexpr SurfaceIntersection() = default;
};

static_assert(std::is_trivially_move_constructible_v<SurfaceIntersection>);
static_assert(std::is_trivially_move_assignable_v<SurfaceIntersection>);

}  // namespace Acts
