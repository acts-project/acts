// This file is part of the Acts project.
//
// Copyright (C) 2016-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>

#include <boost/container/static_vector.hpp>

namespace Acts {

/// Status enum
enum class IntersectionStatus : int {
  missed = 0,
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
  using Position = ActsVector<DIM>;
  /// Status enum
  using Status = IntersectionStatus;

  /// Constructor with arguments
  ///
  /// @param position is the position of the intersection
  /// @param pathLength is the path length to the intersection
  /// @param status is an enum indicating the status of the intersection
  constexpr Intersection(const Position& position, double pathLength,
                         Status status)
      : m_position(position), m_pathLength(pathLength), m_status(status) {}

  /// Returns whether the intersection was successful or not
  constexpr explicit operator bool() const {
    return m_status != Status::missed;
  }

  constexpr const Position& position() const { return m_position; }

  constexpr ActsScalar pathLength() const { return m_pathLength; }

  constexpr Status status() const { return m_status; }

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
    if ((aIntersection.status() == Status::unreachable) &&
        (bIntersection.status() != Status::unreachable)) {
      return false;
    }
    if ((aIntersection.status() != Status::unreachable) &&
        (bIntersection.status() == Status::unreachable)) {
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
  Position m_position = Position::Zero();
  /// Signed path length to the intersection (if valid)
  ActsScalar m_pathLength = std::numeric_limits<double>::infinity();
  /// The Status of the intersection
  Status m_status = Status::unreachable;

  constexpr Intersection() = default;
};

using Intersection2D = Intersection<2>;
using Intersection3D = Intersection<3>;

static constexpr std::uint8_t s_maximumNumberOfIntersections = 2;
using MultiIntersection3D =
    boost::container::static_vector<Intersection3D,
                                    s_maximumNumberOfIntersections>;

template <typename object_t>
class ObjectIntersection {
 public:
  /// Object intersection
  ///
  /// @param intersection is the intersection
  /// @param object is the object to be instersected
  /// @param index is the intersection index
  constexpr ObjectIntersection(const Intersection3D& intersection,
                               const object_t* object, std::uint8_t index = 0)
      : m_intersection(intersection), m_object(object), m_index(index) {}

  /// Returns whether the intersection was successful or not
  constexpr explicit operator bool() const {
    return m_intersection.operator bool();
  }

  constexpr const Intersection3D& intersection() const {
    return m_intersection;
  }

  constexpr const Intersection3D::Position& position() const {
    return m_intersection.position();
  }

  constexpr ActsScalar pathLength() const {
    return m_intersection.pathLength();
  }

  constexpr Intersection3D::Status status() const {
    return m_intersection.status();
  }

  constexpr const object_t* object() const { return m_object; }

  constexpr std::uint8_t index() const { return m_index; }

  constexpr static ObjectIntersection invalid() { return ObjectIntersection(); }

  constexpr static bool pathLengthOrder(
      const ObjectIntersection& aIntersection,
      const ObjectIntersection& bIntersection) {
    return Intersection3D::pathLengthOrder(aIntersection.intersection(),
                                           bIntersection.intersection());
  }

  constexpr static bool closestOrder(const ObjectIntersection& aIntersection,
                                     const ObjectIntersection& bIntersection) {
    return Intersection3D::closestOrder(aIntersection.intersection(),
                                        bIntersection.intersection());
  }

  constexpr static bool closestForwardOrder(
      const ObjectIntersection& aIntersection,
      const ObjectIntersection& bIntersection) {
    return Intersection3D::closestForwardOrder(aIntersection.intersection(),
                                               bIntersection.intersection());
  }

 private:
  /// The intersection itself
  Intersection3D m_intersection = Intersection3D::invalid();
  /// The object that was (tried to be) intersected
  const object_t* m_object = nullptr;
  /// The intersection index
  std::uint8_t m_index = 0;

  constexpr ObjectIntersection() = default;
};

template <typename object_t>
class ObjectMultiIntersection {
 public:
  using SplitIntersections =
      boost::container::static_vector<ObjectIntersection<object_t>,
                                      s_maximumNumberOfIntersections>;

  /// Object intersection
  ///
  /// @param intersections are the intersections
  /// @param object is the object to be instersected
  constexpr ObjectMultiIntersection(const MultiIntersection3D& intersections,
                                    const object_t* object)
      : m_intersections(intersections), m_object(object) {}

  constexpr ObjectIntersection<object_t> operator[](std::uint8_t index) const {
    return {m_intersections[index], m_object, index};
  }

  constexpr const MultiIntersection3D& intersections() const {
    return m_intersections;
  }

  constexpr std::size_t size() const { return m_intersections.size(); }

  constexpr const object_t* object() const { return m_object; }

  constexpr SplitIntersections split() const {
    SplitIntersections result;
    for (std::size_t i = 0; i < size(); ++i) {
      result.push_back(operator[](i));
    }
    return result;
  }

  constexpr ObjectIntersection<object_t> closest() const {
    auto splitIntersections = split();
    return *std::min_element(splitIntersections.begin(),
                             splitIntersections.end(),
                             ObjectIntersection<object_t>::closestOrder);
  }

  constexpr ObjectIntersection<object_t> closestForward() const {
    auto splitIntersections = split();
    return *std::min_element(splitIntersections.begin(),
                             splitIntersections.end(),
                             ObjectIntersection<object_t>::closestForwardOrder);
  }

 private:
  /// The intersections
  MultiIntersection3D m_intersections;
  /// The object that was (tried to be) intersected
  const object_t* m_object = nullptr;
};

namespace detail {

/// This function checks if an intersection is valid for the specified
/// path-limit and overstep-limit
///
/// @tparam intersection_t Type of the intersection object
///
/// @param intersection The intersection to check
/// @param nearLimit The minimum distance for an intersection to be considered
/// @param farLimit The maximum distance for an intersection to be considered
/// @param logger A optionally supplied logger which prints out a lot of infos
///               at VERBOSE level
template <typename intersection_t>
bool checkIntersection(const intersection_t& intersection, double nearLimit,
                       double farLimit,
                       const Logger& logger = getDummyLogger()) {
  const double distance = intersection.pathLength();
  // TODO why?
  const double tolerance = s_onSurfaceTolerance;

  ACTS_VERBOSE(" -> near limit, far limit, distance: "
               << nearLimit << ", " << farLimit << ", " << distance);

  const bool coCriterion = distance > nearLimit;
  const bool cpCriterion = distance < farLimit + tolerance;

  const bool accept = coCriterion && cpCriterion;

  if (accept) {
    ACTS_VERBOSE("Intersection is WITHIN limit");
  } else {
    ACTS_VERBOSE("Intersection is OUTSIDE limit because: ");
    if (!coCriterion) {
      ACTS_VERBOSE("- intersection path length "
                   << distance << " <= near limit " << nearLimit);
    }
    if (!cpCriterion) {
      ACTS_VERBOSE("- intersection path length "
                   << distance << " is over the far limit "
                   << (farLimit + tolerance) << " (including tolerance of "
                   << tolerance << ")");
    }
  }

  return accept;
}

}  // namespace detail

}  // namespace Acts
