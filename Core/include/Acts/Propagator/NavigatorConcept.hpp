// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <concepts>

namespace Acts {

class IVolumeMaterial;

/// @brief Concept that is satisfied by navigators.
template <typename Navigator, typename Options = typename Navigator::Options,
          typename State = typename Navigator::State>
concept NavigatorConcept = requires {
  typename Navigator::Config;
  typename Navigator::Options;
  typename Navigator::State;

  requires requires(const Navigator& n, const Options& o, State& s,
                    const Surface& sf, const Vector3& position,
                    const Vector3& direction, Direction propagationDirection) {
    { n.makeState(o) } -> std::same_as<State>;
    { n.currentSurface(s) } -> std::same_as<const Surface*>;
    { n.currentVolume(s) } -> std::same_as<const TrackingVolume*>;
    { n.currentVolumeMaterial(s) } -> std::same_as<const IVolumeMaterial*>;
    { n.startSurface(s) } -> std::same_as<const Surface*>;
    { n.targetSurface(s) } -> std::same_as<const Surface*>;
    { n.endOfWorldReached(s) } -> std::same_as<bool>;
    { n.navigationBreak(s) } -> std::same_as<bool>;
    {
      n.initialize(s, position, direction, propagationDirection)
    } -> std::same_as<Result<void>>;
    { n.nextTarget(s, position, direction) } -> std::same_as<NavigationTarget>;
    { n.checkTargetValid(s, position, direction) } -> std::same_as<bool>;
    {
      n.handleSurfaceReached(s, position, direction, sf)
    } -> std::same_as<void>;
  };
};

}  // namespace Acts
