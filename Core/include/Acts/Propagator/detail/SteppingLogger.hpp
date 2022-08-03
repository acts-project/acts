// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"

#include <memory>
#include <vector>

namespace Acts {

class Surface;
class Layer;
class TrackingVolume;

namespace detail {

/// @brief the step information for
struct Step {
  ConstrainedStep stepSize;
  Vector3 position = Vector3(0., 0., 0.);
  Vector3 momentum = Vector3(0., 0., 0.);
  std::shared_ptr<const Surface> surface = nullptr;
  const TrackingVolume* volume = nullptr;
};

/// @brief a step length logger for debugging the stepping
///
/// It simply logs the constrained step length per step
struct SteppingLogger {
  /// Simple result struct to be returned
  struct this_result {
    std::vector<Step> steps;
  };

  using result_type = this_result;

  /// Set the Logger to sterile
  bool sterile = false;

  /// SteppingLogger action for the ActionList of the Propagator
  ///
  /// @tparam stepper_t is the type of the Stepper
  /// @tparam propagator_state_t is the type of Propagator state
  ///
  /// @param [in,out] state is the mutable stepper state object
  /// @param [in,out] result is the mutable result object
  template <typename propagator_state_t, typename stepper_t>
  void operator()(propagator_state_t& state, const stepper_t& stepper,
                  result_type& result) const {
    // don't log if you have reached the target
    if (sterile or state.navigation.targetReached) {
      return;
    }
    // record the propagation state
    Step step;
    step.stepSize = state.stepping.stepSize;
    step.position = stepper.position(state.stepping);
    step.momentum =
        stepper.momentum(state.stepping) * stepper.direction(state.stepping);

    if (state.navigation.currentSurface != nullptr) {
      // hang on to surface
      step.surface = state.navigation.currentSurface->getSharedPtr();
    }

    step.volume = state.navigation.currentVolume;
    result.steps.push_back(std::move(step));
  }

  /// Pure observer interface
  /// - this does not apply to the logger
  template <typename propagator_state_t, typename stepper_t>
  void operator()(propagator_state_t& /*unused*/,
                  const stepper_t& /*unused*/) const {}
};

}  // namespace detail
}  // namespace Acts
