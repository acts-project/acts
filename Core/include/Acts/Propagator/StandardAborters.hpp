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
#include <limits>

namespace Acts {

/// This is the condition that the pathLimit has been reached
struct PathLimitReached {
  /// Internal path limit for loop protection
  double internalLimit = std::numeric_limits<double>::max();

  /// boolean operator for abort condition without using the result
  ///
  /// @tparam propagator_state_t Type of the propagator state
  /// @tparam stepper_t Type of the stepper
  /// @tparam navigator_t Type of the navigator
  ///
  /// @param [in,out] state The propagation state object
  /// @param [in] stepper Stepper used for propagation
  /// @param [in] navigator Navigator used for propagation
  /// @param logger a logger instance
  /// @return True if path limit exceeded and propagation should abort
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  bool checkAbort(propagator_state_t& state, const stepper_t& stepper,
                  const navigator_t& navigator, const Logger& logger) const {
    static_cast<void>(navigator);

    // Check if the maximum allowed step size has to be updated
    double distance =
        std::abs(internalLimit) - std::abs(state.stepping.pathAccumulated);
    double tolerance = state.options.surfaceTolerance;
    bool limitReached = (std::abs(distance) < std::abs(tolerance));
    if (limitReached) {
      ACTS_VERBOSE("PathLimit aborter | " << "Path limit reached at distance "
                                          << distance);
      return true;
    }
    stepper.updateStepSize(state.stepping, distance,
                           ConstrainedStep::Type::Actor);
    ACTS_VERBOSE("PathLimit aborter | "
                 << "Target stepSize (path limit) updated to "
                 << stepper.outputStepSize(state.stepping));
    return false;
  }
};

/// Tag used in place of a target aborter type to build a propagator state for
/// a propagation without a target surface
struct NoTargetAborter {};

/// This is the condition that the target surface has been reached. It aborts
/// the propagation once the navigator reports the target as the current
/// surface, and it does not steer the propagation towards it.
///
/// @note The navigator steers the propagation onto the target surface. The
///       propagator hands its target surface to the navigator. An actor that
///       stops on a surface of its own registers it as an additional surface
///       of the navigator.
struct SurfaceReached {
  /// Target surface to reach for propagation termination
  const Surface* surface = nullptr;

  /// boolean operator for abort condition without using the result
  ///
  /// @tparam propagator_state_t Type of the propagator state
  /// @tparam stepper_t Type of the stepper
  /// @tparam navigator_t Type of the navigator
  ///
  /// @param [in,out] state The propagation state object
  /// @param [in] stepper Stepper used for propagation
  /// @param [in] navigator Navigator used for propagation
  /// @param logger a logger instance
  /// @return true if abort condition is met (surface reached)
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  bool checkAbort(propagator_state_t& state, const stepper_t& stepper,
                  const navigator_t& navigator, const Logger& logger) const {
    static_cast<void>(stepper);

    if (surface == nullptr ||
        navigator.currentSurface(state.navigation) != surface) {
      return false;
    }

    ACTS_VERBOSE("SurfaceReached aborter | Target surface reached.");
    return true;
  }
};

/// This is the condition that the end of world has been reached
/// it then triggers an propagation abort
struct EndOfWorldReached {
  /// boolean operator for abort condition without using the result
  ///
  /// @tparam propagator_state_t Type of the propagator state
  /// @tparam navigator_t Type of the navigator
  ///
  /// @param [in,out] state The propagation state object
  /// @param [in] navigator The navigator object
  /// @return True if end of world reached and propagation should abort
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  bool checkAbort(propagator_state_t& state, const stepper_t& /*stepper*/,
                  const navigator_t& navigator,
                  const Logger& /*logger*/) const {
    bool endOfWorld = navigator.endOfWorldReached(state.navigation);
    return endOfWorld;
  }
};

/// This is the condition that the end of world has been reached
/// it then triggers a propagation abort
struct VolumeConstraintAborter {
  /// boolean operator for abort condition without using the result
  ///
  /// @tparam propagator_state_t Type of the propagator state
  /// @tparam navigator_t Type of the navigator
  ///
  /// @param [in,out] state The propagation state object
  /// @param [in] navigator The navigator object
  /// @param logger a logger instance
  /// @return True if volume constraints violated and propagation should abort
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  bool checkAbort(propagator_state_t& state, const stepper_t& /*stepper*/,
                  const navigator_t& navigator, const Logger& logger) const {
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
      ACTS_VERBOSE(
          "VolumeConstraintAborter aborter | Abort with volume constrain "
          << currentVolumeId);
      return true;
    }

    if (!endOfWorldVolumeIds.empty() &&
        rangeContainsValue(endOfWorldVolumeIds, currentVolumeId)) {
      ACTS_VERBOSE(
          "VolumeConstraintAborter aborter | Abort with additional end of "
          "world volume "
          << currentVolumeId);
      return true;
    }

    return false;
  }
};

/// Aborter that checks if the propagation has reached any surface
struct AnySurfaceReached {
  /// Check if any surface has been reached during propagation
  /// @tparam propagator_state_t Type of the propagator state
  /// @tparam stepper_t Type of the stepper
  /// @tparam navigator_t Type of the navigator
  /// @param state The propagation state object
  /// @param stepper Stepper used for propagation (unused)
  /// @param navigator Navigator used for propagation
  /// @param logger Logger instance (unused)
  /// @return true if any surface has been reached, false otherwise
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  bool checkAbort(propagator_state_t& state, const stepper_t& stepper,
                  const navigator_t& navigator, const Logger& logger) const {
    static_cast<void>(stepper);
    static_cast<void>(logger);

    const Surface* startSurface = navigator.startSurface(state.navigation);
    const Surface* targetSurface = navigator.targetSurface(state.navigation);
    const Surface* currentSurface = navigator.currentSurface(state.navigation);

    // `startSurface` is excluded because we want to reach a new surface
    // `targetSurface` is excluded because another aborter should handle it
    if (currentSurface != nullptr && currentSurface != startSurface &&
        currentSurface != targetSurface) {
      return true;
    }

    return false;
  }
};

}  // namespace Acts
