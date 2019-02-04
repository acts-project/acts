// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
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
  struct VoidNavigator
  {

    /// @brief Nested State struct, minimal requirement
    struct State
    {
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

    /// Navigation call - void
    ///
    /// @tparam propagator_state_t is the type of Propagatgor state
    /// @tparam stepper_t Type of the Stepper
    ///
    /// Empty call, compiler should optimise that
    template <typename propagator_state_t, typename stepper_t>
    void
    status(propagator_state_t& /*state*/, const stepper_t& /*unused*/) const
    {
    }

    /// Navigation call - void
    ///
    /// @tparam propagator_state_t is the type of Propagatgor state
    /// @tparam stepper_t Type of the Stepper
    ///
    /// Empty call, compiler should optimise that
    template <typename propagator_state_t, typename stepper_t>
    void
    target(propagator_state_t& /*state*/, const stepper_t& /*unused*/) const
    {
    }

    /// Navigation call - void
    ///
    /// @tparam propagator_state_t is the type of Propagatgor state
    ///
    /// Empty call, compiler should optimise that
    template <typename propagator_state_t>
    void
    operator()(propagator_state_t& /*state*/) const
    {
    }
  };

}  // namespace detail
}  // namespace Acts
