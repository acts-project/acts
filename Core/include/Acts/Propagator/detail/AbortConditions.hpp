// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cmath>
#include <cstdint>

namespace Acts::detail {

/// Check if the propagation reached a path limit
///
/// If the limit is not reached yet, this constrains the next step to it.
///
/// @param [in,out] state The propagation state
/// @param [in] stepper Stepper used for the propagation
/// @param [in] pathLimit Signed path limit
/// @param [in] logger A logger instance
/// @return True if the path limit is reached
template <typename propagator_state_t, typename stepper_t>
bool pathLimitReached(propagator_state_t& state, const stepper_t& stepper,
                      double pathLimit, const Logger& logger) {
  const double distance =
      std::abs(pathLimit) - std::abs(stepper.pathLength(state.stepping));
  const double tolerance = state.options.surfaceTolerance;
  if (std::abs(distance) < std::abs(tolerance)) {
    ACTS_VERBOSE("Path limit reached at distance " << distance);
    return true;
  }
  stepper.updateStepSize(state.stepping, distance,
                         ConstrainedStep::Type::Actor);
  ACTS_VERBOSE("Step size updated to the path limit: "
               << stepper.outputStepSize(state.stepping));
  return false;
}

/// Check if the propagation reached a surface
///
/// @param [in] state The propagation state
/// @param [in] navigator Navigator used for the propagation
/// @param [in] surface The surface to reach, or nullptr for no surface
/// @param [in] logger A logger instance
/// @return True if the surface is the current surface of the navigator
template <typename propagator_state_t, typename navigator_t>
bool surfaceReached(propagator_state_t& state, const navigator_t& navigator,
                    const Surface* surface, const Logger& logger) {
  if (surface == nullptr ||
      navigator.currentSurface(state.navigation) != surface) {
    return false;
  }
  ACTS_VERBOSE("Target surface reached.");
  return true;
}

/// Check if the propagation left the world
///
/// @param [in] state The propagation state
/// @param [in] navigator Navigator used for the propagation
/// @return True if the navigator reached the end of the world
template <typename propagator_state_t, typename navigator_t>
bool endOfWorldReached(propagator_state_t& state,
                       const navigator_t& navigator) {
  return navigator.endOfWorldReached(state.navigation);
}

/// Check if the propagation left the allowed volumes
///
/// This uses `constrainToVolumeIds` and `endOfWorldVolumeIds` of the
/// propagator options.
///
/// @param [in] state The propagation state
/// @param [in] navigator Navigator used for the propagation
/// @param [in] logger A logger instance
/// @return True if the current volume is not allowed
template <typename propagator_state_t, typename navigator_t>
bool volumeConstraintReached(propagator_state_t& state,
                             const navigator_t& navigator,
                             const Logger& logger) {
  const auto& constrainToVolumeIds = state.options.constrainToVolumeIds;
  const auto& endOfWorldVolumeIds = state.options.endOfWorldVolumeIds;

  if (constrainToVolumeIds.empty() && endOfWorldVolumeIds.empty()) {
    return false;
  }
  const auto* currentVolume = navigator.currentVolume(state.navigation);

  // We need a volume to check its ID
  if (currentVolume == nullptr) {
    return false;
  }

  const auto currentVolumeId =
      static_cast<std::uint32_t>(currentVolume->geometryId().volume());

  if (!constrainToVolumeIds.empty() &&
      !rangeContainsValue(constrainToVolumeIds, currentVolumeId)) {
    ACTS_VERBOSE("Abort with volume constraint " << currentVolumeId);
    return true;
  }

  if (!endOfWorldVolumeIds.empty() &&
      rangeContainsValue(endOfWorldVolumeIds, currentVolumeId)) {
    ACTS_VERBOSE("Abort with additional end of world volume "
                 << currentVolumeId);
    return true;
  }

  return false;
}

/// Check if the propagation reached a surface other than the start and the
/// target surface
///
/// @param [in] state The propagation state
/// @param [in] navigator Navigator used for the propagation
/// @return True if the current surface is neither the start nor the target
template <typename propagator_state_t, typename navigator_t>
bool anySurfaceReached(propagator_state_t& state,
                       const navigator_t& navigator) {
  const Surface* startSurface = navigator.startSurface(state.navigation);
  const Surface* targetSurface = navigator.targetSurface(state.navigation);
  const Surface* currentSurface = navigator.currentSurface(state.navigation);

  // `startSurface` is excluded because we want to reach a new surface
  // `targetSurface` is excluded because another condition should handle it
  return currentSurface != nullptr && currentSurface != startSurface &&
         currentSurface != targetSurface;
}

}  // namespace Acts::detail
