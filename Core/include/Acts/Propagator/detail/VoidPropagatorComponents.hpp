// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

namespace Acts {
namespace detail {

/// @brief The void navigator struct as a default navigator
///
/// It does not provide any navigation action, the compiler
/// should eventually optimise that the function call is not done
///
struct VoidNavigator {
  /// @brief Nested State struct, minimal requirement
  struct State {
    /// Navigation state - external state: the start surface
    const Surface* startSurface = nullptr;

    /// Navigation state - external state: the current surface
    const Surface* currentSurface = nullptr;

    /// Navigation state - external state: the target surface
    const Surface* targetSurface = nullptr;

    /// Indicator if the target is reached
    bool targetReached = false;

    /// Navigation state : a break has been detected
    bool navigationBreak = false;
  };

  /// Unique typedef to publish to the Propagator
  using state_type = State;

  State makeState(const Surface* startSurface,
                  const Surface* targetSurface) const {
    State result;
    result.startSurface = startSurface;
    result.targetSurface = targetSurface;
    return result;
  }

  const Surface* currentSurface(const State& state) const {
    return state.currentSurface;
  }

  const Surface* startSurface(const State& state) const {
    return state.startSurface;
  }

  // TODO sounds like an aborters job?
  const Surface* targetSurface(const State& state) const {
    return state.targetSurface;
  }

  // TODO sounds like an aborters job?
  bool targetReached(const State& state) const { return state.targetReached; }

  // TODO not sure why the navigation has to break - propergator has control
  // over that
  bool navigationBreak(const State& state) const {
    return state.navigationBreak;
  }

  // TODO why would anybody set this? sounds dangerous
  void currentSurface(State& state, const Surface* surface) const {
    state.currentSurface = surface;
  }

  // TODO sounds like an aborters job? shouldnt we tell the propagator that we
  // are done?
  void targetReached(State& state, bool targetReached) const {
    state.targetReached = targetReached;
  }

  // TODO not sure why the navigation has to break - propergator has control
  // over that
  void navigationBreak(State& state, bool navigationBreak) const {
    state.navigationBreak = navigationBreak;
  }

  /// Navigation call - void
  ///
  /// @tparam propagator_state_t is the type of Propagatgor state
  /// @tparam stepper_t Type of the Stepper
  ///
  /// Empty call, compiler should optimise that
  template <typename propagator_state_t, typename stepper_t>
  void status(propagator_state_t& /*state*/,
              const stepper_t& /*stepper*/) const {}

  /// Navigation call - void
  ///
  /// @tparam propagator_state_t is the type of Propagatgor state
  /// @tparam stepper_t Type of the Stepper
  ///
  /// Empty call, compiler should optimise that
  template <typename propagator_state_t, typename stepper_t>
  void target(propagator_state_t& /*state*/,
              const stepper_t& /*stepper*/) const {}
};

}  // namespace detail
}  // namespace Acts
