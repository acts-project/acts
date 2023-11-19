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
  /// @tparam navigator_t Type of the navigator
  ///
  /// @param [in,out] state The propagation state object
  /// @param [in] stepper Stepper used for propagation
  /// @param [in] navigator Navigator used for the propagation
  /// @param logger a logger instance
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  bool operator()(propagator_state_t& state, const stepper_t& stepper,
                  const navigator_t& navigator, const Logger& logger) const {
    return (*this)(state, stepper, navigator,
                   *navigator.targetSurface(state.navigation), logger);
  }

  /// boolean operator for abort condition without using the result
  ///
  /// @tparam propagator_state_t Type of the propagator state
  /// @tparam stepper_t Type of the stepper
  /// @tparam navigator_t Type of the navigator
  ///
  /// @param [in,out] state The propagation state object
  /// @param [in] stepper Stepper used for the propagation
  /// @param [in] navigator Navigator used for the propagation
  /// @param [in] targetSurface The target surface
  /// @param logger a logger instance
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  bool operator()(propagator_state_t& state, const stepper_t& stepper,
                  const navigator_t& navigator, const Surface& targetSurface,
                  const Logger& logger) const {
    bool reached = true;
    const auto oldCurrentSurface = navigator.currentSurface(state.navigation);

    for (auto cmp : stepper.componentIterable(state.stepping)) {
      auto singleState = cmp.singleState(state);
      const auto& singleStepper = cmp.singleStepper(stepper);

      if (!SurfaceReached{}(singleState, singleStepper, navigator,
                            targetSurface, logger)) {
        reached = false;
      } else {
        cmp.status() = Acts::Intersection3D::Status::onSurface;
      }
    }

    // However, if mean of all is on surface, we are happy as well
    if (averageOnSurface) {
      const auto sIntersection =
          targetSurface
              .intersect(
                  state.geoContext, stepper.position(state.stepping),
                  state.options.direction * stepper.direction(state.stepping),
                  BoundaryCheck(true), averageOnSurfaceTolerance)
              .closest();

      if (sIntersection.status() == Intersection3D::Status::onSurface) {
        ACTS_VERBOSE("Reached target in average mode");
        navigator.currentSurface(state.navigation, &targetSurface);
        navigator.targetReached(state.navigation, true);

        for (auto cmp : stepper.componentIterable(state.stepping)) {
          cmp.status() = Intersection3D::Status::onSurface;
        }

        return true;
      }
    }

    // These values are changed by the single component aborters but must be
    // reset if not all components are on the target
    if (!reached) {
      navigator.currentSurface(state.navigation, oldCurrentSurface);
      navigator.targetReached(state.navigation, false);
    }

    return reached;
  }
};
}  // namespace Acts
