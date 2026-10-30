// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Propagator/detail/AbortConditions.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"

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
    return detail::pathLimitReached(state, stepper, internalLimit, logger);
  }
};

/// Tag used in place of a target aborter type to build a propagator state for
/// a propagation without a target surface
struct NoTargetAborter {};

/// Aborts once the navigator reports the target surface as the current surface.
///
/// @note It does not steer. The navigator does, so an actor with its own target
///       registers it as an additional surface of the navigator.
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
  /// @param [in] navigator Navigator used for propagation
  /// @param logger a logger instance
  /// @return true if abort condition is met (surface reached)
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  bool checkAbort(propagator_state_t& state, const stepper_t& /*stepper*/,
                  const navigator_t& navigator, const Logger& logger) const {
    return detail::surfaceReached(state, navigator, surface, logger);
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
    return detail::endOfWorldReached(state, navigator);
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
    return detail::volumeConstraintReached(state, navigator, logger);
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
    return detail::anySurfaceReached(state, navigator);
  }
};

}  // namespace Acts
