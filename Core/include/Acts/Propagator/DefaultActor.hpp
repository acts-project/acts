// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Propagator/PropagatorState.hpp"
#include "Acts/Propagator/detail/AbortConditions.hpp"
#include "Acts/Propagator/detail/LoopProtection.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <limits>

namespace Acts {

/// @brief Actor that stops the propagation where the configuration says
///
/// The propagator runs this actor if the options give no other actor. It stops
/// the propagation
/// - at the target surface of the propagation, if there is one,
/// - at the path limit of the options, lowered by the loop protection,
/// - outside the volumes that the options allow,
/// - at the end of the world, if the options say so.
///
/// An actor of your own replaces this actor. To keep these conditions, put
/// this actor in an @c ActorList with your actor. Put it first, so the loop
/// protection sees the start momentum before other actors change it.
struct DefaultActor {
  /// State of the actor
  struct State {
    /// Signed path limit, set at the start of the propagation
    double pathLimit = std::numeric_limits<double>::max();
  };

  /// Type alias for the actor state
  using state_type = State;

  /// Set the path limit at the start of the propagation
  ///
  /// @param [in,out] state The propagation state
  /// @param [in] stepper Stepper used for the propagation
  /// @param [in] navigator Navigator used for the propagation
  /// @param [in,out] actorState The state of the actor
  /// @param [in] logger A logger instance
  /// @return Always success
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  Result<void> act(propagator_state_t& state, const stepper_t& stepper,
                   const navigator_t& navigator, State& actorState,
                   const Logger& logger) const {
    static_cast<void>(navigator);
    if (state.stage == PropagatorStage::prePropagation) {
      actorState.pathLimit = state.options.pathLimit;
      detail::setupLoopProtection(state, stepper, actorState.pathLimit, false,
                                  logger);
    }
    return Result<void>::success();
  }

  /// Check the abort conditions
  ///
  /// @param [in,out] state The propagation state
  /// @param [in] stepper Stepper used for the propagation
  /// @param [in] navigator Navigator used for the propagation
  /// @param [in] actorState The state of the actor
  /// @param [in] logger A logger instance
  /// @return True if the propagation should stop
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  bool checkAbort(propagator_state_t& state, const stepper_t& stepper,
                  const navigator_t& navigator, const State& actorState,
                  const Logger& logger) const {
    return (state.options.stopAtEndOfWorld &&
            detail::endOfWorldReached(state, navigator)) ||
           detail::volumeConstraintReached(state, navigator, logger) ||
           detail::surfaceReached(state, navigator,
                                  navigator.targetSurface(state.navigation),
                                  logger) ||
           detail::pathLimitReached(state, stepper, actorState.pathLimit,
                                    logger);
  }
};

}  // namespace Acts
