// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Propagator/detail/AnyIntersection.hpp"

#include <variant>

namespace Acts {
namespace detail {

struct NavigationCandidate {
  using SurfaceObject = const Surface*;
  using LayerObject = const Layer*;
  using BoundaryObject = const BoundarySurfaceT<TrackingVolume>*;
  using AnyObject = std::variant<SurfaceObject, LayerObject, BoundaryObject>;

  AnyObject object;
  const Surface* representation = nullptr;
  BoundaryCheck boundaryCheck;

  NavigationCandidate(AnyObject _object, const Surface* _representation,
                      BoundaryCheck _boundaryCheck)
      : object(_object),
        representation(_representation),
        boundaryCheck(std::move(_boundaryCheck)) {}

  AnyMultiIntersection intersect(const GeometryContext& gctx,
                                 const Vector3& position,
                                 const Vector3& direction,
                                 ActsScalar tolerance) const {
    if (std::holds_alternative<SurfaceObject>(object)) {
      const auto& surface = std::get<SurfaceObject>(object);
      auto intersection = representation->intersect(gctx, position, direction,
                                                    boundaryCheck, tolerance);
      return AnyMultiIntersection(SurfaceMultiIntersection(
          intersection.intersections(), surface, representation));
    }
    if (std::holds_alternative<LayerObject>(object)) {
      const auto& layer = std::get<LayerObject>(object);
      auto intersection = representation->intersect(gctx, position, direction,
                                                    boundaryCheck, tolerance);
      return AnyMultiIntersection(LayerMultiIntersection(
          intersection.intersections(), layer, representation));
    }
    if (std::holds_alternative<BoundaryObject>(object)) {
      const auto& boundary = std::get<BoundaryObject>(object);
      auto intersection = boundary->surfaceRepresentation().intersect(
          gctx, position, direction, boundaryCheck, tolerance);
      return AnyMultiIntersection(BoundaryMultiIntersection(
          intersection.intersections(), boundary, representation));
    }
    throw std::runtime_error("unknown type");
  }
};

}  // namespace detail
}  // namespace Acts
