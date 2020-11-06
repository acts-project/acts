// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Propagator/ConstrainedStep.hpp"

namespace Acts {

/// @brief Helper struct for constrained step handling.
///
/// It handles the update, resetting of the step size,
/// and step size components in the Stepper.
///
template <typename stepper_t>
struct ConstrainedStepControl {
  /// Retrieve the step size
  ///
  /// @param state [in,out] The stepping state (thread-local cache)
  double size(const typename stepper_t::State& state,
              ConstrainedStep::Type stype = ConstrainedStep::actor) const {
    return state.stepSize.value(stype);
  }

  /// Update step size
  ///
  /// It checks the status to the reference surface & updates
  /// the step size accordingly
  ///
  /// @param state [in,out] The stepping state (thread-local cache)
  /// @param oIntersection [in] The ObjectIntersection to layer, boundary, etc
  /// @param release [in] boolean to trigger step size release
  template <typename object_intersection_t>
  void update(typename stepper_t::State& state,
              const object_intersection_t& oIntersection,
              bool release = true) const {
    state.stepSize.update(oIntersection.intersection.pathLength,
                          ConstrainedStep::actor, release);
  }

  /// Update step size component
  ///
  /// This method intersects the provided surface and update the navigation
  /// step estimation accordingly (hence it changes the state). It also
  /// returns the status of the intersection to trigger onSurface in case
  /// the surface is reached.
  ///
  /// @param state [in,out] The stepping state (thread-local cache)
  /// @param stepSize [in] The new step size to be set
  /// @param stype [in] The type of the step size component to be changed
  void update(typename stepper_t::State& state, double stepSize,
              ConstrainedStep::Type stype = ConstrainedStep::actor) const {
    state.stepSize.update(stepSize, stype);
  }

  /// Set Step size - explicitely with a double
  ///
  /// @param state [in,out] The stepping state (thread-local cache)
  /// @param stepSize [in] The step size value
  /// @param stype [in] The step size type to be set
  void set(typename stepper_t::State& state, double stepSize,
           ConstrainedStep::Type stype = ConstrainedStep::actor) const {
    state.previousStepSize = state.stepSize;
    state.stepSize.update(stepSize, stype, true);
  }

  /// Release the Step size
  ///
  /// @param state [in,out] The stepping state (thread-local cache)
  /// @param stype [in] The step size type to be set
  void release(typename stepper_t::State& state,
               ConstrainedStep::Type stype = ConstrainedStep::actor) const {
    state.stepSize.release(stype);
  }

  /// Output the Step Size - single component
  ///
  /// @param state [in,out] The stepping state (thread-local cache)
  std::string output(const typename stepper_t::State& state) const {
    return state.stepSize.toString();
  }
};

}  // namespace Acts
