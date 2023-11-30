// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <variant>

namespace Acts {
namespace detail {

/// @brief Type erased intersection
///
/// TODO this would be obsolete if the intersection could carry a `std::any`
class AnyIntersection {
 public:
  using Any = std::variant<SurfaceIntersection, LayerIntersection,
                           BoundaryIntersection>;

  AnyIntersection(const Any& any) : m_any(any) {}
  AnyIntersection(Any&& any) : m_any(std::move(any)) {}

  explicit operator bool() const {
    return std::visit(
        [](const auto& intersection) { return intersection.operator bool(); },
        m_any);
  }

  template <typename intersection_t>
  bool checkType() const {
    return std::holds_alternative<intersection_t>(m_any);
  }

  template <typename intersection_t>
  const typename intersection_t::Object* object() const {
    return std::get<intersection_t>(m_any).object();
  }

  const Surface* representation() const {
    return std::visit(
        [](const auto& intersection) -> const Surface* {
          return intersection.representation();
        },
        m_any);
  }

  const Intersection3D& intersection() const {
    return std::visit(
        [](const auto& intersection) -> const Intersection3D& {
          return intersection.intersection();
        },
        m_any);
  }

  const Intersection3D::Position& position() const {
    return std::visit(
        [](const auto& intersection) -> const Intersection3D::Position& {
          return intersection.position();
        },
        m_any);
  }

  ActsScalar pathLength() const {
    return std::visit(
        [](const auto& intersection) { return intersection.pathLength(); },
        m_any);
  }

  Intersection3D::Status status() const {
    return std::visit(
        [](const auto& intersection) { return intersection.status(); }, m_any);
  }

  std::uint8_t index() const {
    return std::visit(
        [](const auto& intersection) { return intersection.index(); }, m_any);
  }

  static bool forwardOrder(const AnyIntersection& aIntersection,
                           const AnyIntersection& bIntersection) {
    return Intersection3D::forwardOrder(aIntersection.intersection(),
                                        bIntersection.intersection());
  }

  static bool closestOrder(const AnyIntersection& aIntersection,
                           const AnyIntersection& bIntersection) {
    return Intersection3D::closestOrder(aIntersection.intersection(),
                                        bIntersection.intersection());
  }

 private:
  Any m_any;
};

/// @brief Type erased any intersection
///
/// TODO this would be obsolete if the intersection could carry a `std::any`
class AnyMultiIntersection {
 public:
  using Any = std::variant<SurfaceMultiIntersection, LayerMultiIntersection,
                           BoundaryMultiIntersection>;
  using SplitIntersections =
      boost::container::static_vector<AnyIntersection,
                                      s_maximumNumberOfIntersections>;

  AnyMultiIntersection(const Any& any) : m_any(any) {}
  AnyMultiIntersection(Any&& any) : m_any(std::move(any)) {}

  AnyIntersection operator[](std::uint8_t index) const {
    return std::visit(
        [&](const auto& intersection) {
          return AnyIntersection(intersection[index]);
        },
        m_any);
  }

  template <typename intersection_t>
  bool checkType() const {
    return std::holds_alternative<intersection_t>(m_any);
  }

  std::size_t size() const {
    return std::visit(
        [](const auto& intersection) { return intersection.size(); }, m_any);
  }

  template <typename intersection_t>
  const typename intersection_t::Object* object() const {
    return std::get<intersection_t>(m_any).object();
  }

  const Surface* representation() const {
    return std::visit(
        [](const auto& intersection) -> const Surface* {
          return intersection.representation();
        },
        m_any);
  }

  SplitIntersections split() const {
    SplitIntersections result;
    for (std::size_t i = 0; i < size(); ++i) {
      result.push_back(operator[](i));
    }
    return result;
  }

 private:
  Any m_any;
};

}  // namespace detail
}  // namespace Acts
