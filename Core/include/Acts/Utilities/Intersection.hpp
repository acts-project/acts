// This file is part of the Acts project.
//
// Copyright (C) 2016-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <array>
#include <limits>

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
struct Intersection {
  /// Status enum
  using Status = IntersectionStatus;

  /// Position of the intersection
  ActsVector<DIM> position = ActsVector<DIM>::Zero();
  /// Signed path length to the intersection (if valid)
  ActsScalar pathLength{std::numeric_limits<double>::infinity()};
  /// The Status of the intersection
  Status status{Status::unreachable};
  /// The intersection index
  std::uint8_t index{0};

  /// Default constructor
  Intersection() = default;

  /// Constructor with arguments
  ///
  /// @param sinter is the position of the intersection
  /// @param slength is the path length to the intersection
  /// @param sstatus is an enum indicating the status of the intersection
  Intersection(const ActsVector<DIM>& sinter, double slength, Status sstatus,
               std::uint8_t sindex = 0)
      : position(sinter), pathLength(slength), status(sstatus), index(sindex) {}

  /// Returns wheter the intersection was successful or not
  explicit operator bool() const { return status != Status::missed; }
};

using Intersection2D = Intersection<2>;

using Intersection3D = Intersection<3>;

/// @brief class extensions to return also the object and a representation
template <typename object_t, typename representation_t = object_t>
class ObjectIntersection {
 public:
  /// The intersection itself
  Intersection3D intersection{};
  /// The object that was (tried to be) intersected
  const object_t* object{nullptr};
  /// The representation of this object
  const representation_t* representation{nullptr};

  /// Default constructor
  ObjectIntersection() = default;

  /// Object intersection - symmetric setup
  ///
  /// @param sIntersection is the intersection
  /// @param sObject is the object to be instersected
  template <typename T = representation_t,
            std::enable_if_t<std::is_same<T, object_t>::value, int> = 0>
  ObjectIntersection(const Intersection3D& sIntersection,
                     const object_t* sObject)
      : intersection(sIntersection), object(sObject), representation(sObject) {}

  /// Object intersection
  ///
  /// @param sIntersection is the intersection
  /// @param sObject is the object to be instersected
  /// @param sRepresentation is the object representation
  ObjectIntersection(const Intersection3D& sIntersection,
                     const object_t* sObject,
                     const representation_t* sRepresentation)
      : intersection(sIntersection),
        object(sObject),
        representation(sRepresentation) {}

  /// Returns wheter the intersection was successful or not
  explicit operator bool() const { return intersection.operator bool(); }
};

/// @brief class extensions to return also the object and a representation
template <typename object_t, typename representation_t = object_t>
class ObjectMultiIntersection {
 public:
  /// The intersections
  std::array<Intersection3D, 2> intersections;
  /// The object that was (tried to be) intersected
  const object_t* object{nullptr};
  /// The representation of this object
  const representation_t* representation{nullptr};

  /// Default constructor
  ObjectMultiIntersection() = default;

  /// Object intersection - symmetric setup
  ///
  /// @param sIntersections are the intersections
  /// @param sObject is the object to be instersected
  template <typename T = representation_t,
            std::enable_if_t<std::is_same<T, object_t>::value, int> = 0>
  ObjectMultiIntersection(const std::array<Intersection3D, 2>& sIntersections,
                          const object_t* sObject)
      : intersections(sIntersections),
        object(sObject),
        representation(sObject) {
    assert(checkOrder());
  }

  /// Object intersection
  ///
  /// @param sIntersections are the intersections
  /// @param sObject is the object to be instersected
  /// @param sRepresentation is the object representation
  ObjectMultiIntersection(const std::array<Intersection3D, 2>& sIntersections,
                          const object_t* sObject,
                          const representation_t* sRepresentation)
      : intersections(sIntersections),
        object(sObject),
        representation(sRepresentation) {
    assert(checkOrder());
  }

  bool checkOrder() const {
    return intersections[0].pathLength <= intersections[1].pathLength;
  }

  ObjectIntersection<object_t, representation_t> first() const {
    return {intersections[0], object, representation};
  }

  ObjectIntersection<object_t, representation_t> second() const {
    return {intersections[1], object, representation};
  }
};

namespace detail {

/// This function checks if an intersection is valid for the specified
/// path-limit and overstep-limit
///
/// @tparam intersection_t Type of the intersection object
/// @tparam logger_t The logger type, which defaults to std::false_type to
/// prevent the generation of logging code
///
/// @param intersection The intersection to check
/// @param pLimit The path-limit
/// @param oLimit The overstep-limit
/// @param tolerance The tolerance that is applied to the path-limit criterion
/// @param logger A optionally supplied logger which prints out a lot of infos
/// at VERBOSE level
template <typename intersection_t, typename logger_t = std::false_type>
bool checkIntersection(const intersection_t& intersection, double pLimit,
                       double oLimit, double tolerance,
                       const Logger& logger = getDummyLogger()) {
  const double cLimit = intersection.pathLength;

  ACTS_VERBOSE(" -> pLimit, oLimit, cLimit: " << pLimit << ", " << oLimit
                                              << ", " << cLimit);

  const bool coCriterion = cLimit > oLimit;
  const bool cpCriterion = std::abs(cLimit) < std::abs(pLimit) + tolerance;

  const bool accept = coCriterion and cpCriterion;

  if (accept) {
    ACTS_VERBOSE("Intersection is WITHIN limit");
  } else {
    ACTS_VERBOSE("Intersection is OUTSIDE limit because: ");
    if (not coCriterion) {
      ACTS_VERBOSE("- intersection path length "
                   << cLimit << " <= overstep limit " << oLimit);
    }
    if (not cpCriterion) {
      ACTS_VERBOSE("- intersection path length "
                   << std::abs(cLimit) << " is over the path limit "
                   << (std::abs(pLimit) + tolerance)
                   << " (including tolerance of " << tolerance << ")");
    }
  }

  return accept;
}

}  // namespace detail

}  // namespace Acts
