// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Surfaces/Surface.hpp"

namespace Acts {

/// This
struct MultiStepperSurfaceReached {
  MultiStepperSurfaceReached() = default;

  /// boolean operator for abort condition without using the result
  ///
  /// @tparam propagator_state_t Type of the propagator state
  /// @tparam stepper_t Type of the stepper
  ///
  /// @param [in,out] state The propagation state object
  /// @param [in] stepper Stepper used for propagation
  template <typename propagator_state_t, typename stepper_t>
  bool operator()(propagator_state_t& state, const stepper_t& stepper) const {
    return (*this)(state, stepper, *state.navigation.targetSurface);
  }

  /// boolean operator for abort condition without using the result
  ///
  /// @tparam propagator_state_t Type of the propagator state
  /// @tparam stepper_t Type of the stepper
  ///
  /// @param [in,out] state The propagation state object
  /// @param [in] stepper Stepper used for the progation
  /// @param [in] targetSurface The target surface
  template <typename propagator_state_t, typename stepper_t>
  bool operator()(propagator_state_t& state, const stepper_t& stepper,
                  const Surface& targetSurface) const {
    bool reached = true;
    const auto oldCurrentSurface = state.navigation.currentSurface;

    for (auto cmp : stepper.componentIterable(state.stepping)) {
      auto singleState = cmp.singleState(state);
      const auto& singleStepper = cmp.singleStepper(stepper);

      if (!SurfaceReached{}(singleState, singleStepper, targetSurface)) {
        reached = false;
      }
    }

    // These values are changed by the single component aborters but must be
    // reset if not all components are on the target
    if (!reached) {
      state.navigation.currentSurface = oldCurrentSurface;
      state.navigation.targetReached = false;
    }

    return reached;
  }
};
}  // namespace Acts
