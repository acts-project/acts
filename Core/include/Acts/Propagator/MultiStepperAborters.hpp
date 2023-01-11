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

  /// If this is set, we are also happy if the mean of the components is on the
  /// surface. How the averaging is performed depends on the stepper
  /// implementation
  bool averageOnSurface = true;

  /// A configurable tolerance within which distance to the intersection we
  /// consider the surface as reached. Has no effect if averageOnSurface is
  /// false
  double averageOnSurfaceTolerance = 0.2;

  /// boolean operator for abort condition without using the result
  ///
  /// @tparam propagator_state_t Type of the propagator state
  /// @tparam stepper_t Type of the stepper
  ///
  /// @param [in,out] state The propagation state object
  /// @param [in] stepper Stepper used for propagation
  /// @param logger a logger instance
  template <typename propagator_state_t, typename stepper_t>
  bool operator()(propagator_state_t& state, const stepper_t& stepper,
                  const Logger& logger) const {
    return (*this)(state, stepper, *state.navigation.targetSurface, logger);
  }

  /// boolean operator for abort condition without using the result
  ///
  /// @tparam propagator_state_t Type of the propagator state
  /// @tparam stepper_t Type of the stepper
  ///
  /// @param [in,out] state The propagation state object
  /// @param [in] stepper Stepper used for the progation
  /// @param [in] targetSurface The target surface
  /// @param logger a logger instance
  template <typename propagator_state_t, typename stepper_t>
  bool operator()(propagator_state_t& state, const stepper_t& stepper,
                  const Surface& targetSurface, const Logger& logger) const {
    bool reached = true;
    const auto oldCurrentSurface = state.navigation.currentSurface;

    for (auto cmp : stepper.componentIterable(state.stepping)) {
      auto singleState = cmp.singleState(state);
      const auto& singleStepper = cmp.singleStepper(stepper);

      if (!SurfaceReached{}(singleState, singleStepper, targetSurface,
                            logger)) {
        reached = false;
      }
    }

    // However, if mean of all is on surface, we are happy as well
    if (averageOnSurface) {
      const auto sIntersection = targetSurface.intersect(
          state.geoContext, stepper.position(state.stepping),
          state.stepping.navDir * stepper.direction(state.stepping), true);

      if (sIntersection.intersection.status ==
              Intersection3D::Status::onSurface or
          sIntersection.intersection.pathLength < averageOnSurfaceTolerance) {
        ACTS_VERBOSE("Reached target in average mode");
        state.navigation.currentSurface = &targetSurface;
        state.navigation.targetReached = true;
        return true;
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
