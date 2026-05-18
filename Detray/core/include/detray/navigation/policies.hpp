// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// detray definitions
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/units.hpp"
#include "detray/propagator/base_actor.hpp"

// system includes
#include <climits>

namespace detray {

/// Struct that represents the most conservative navigation policy: always re-
/// initialize the current volume
struct always_init : public base_actor {
  struct state {};

  /// Sets the navigation trust level to 'no trust'
  ///
  /// @param pol_state not used
  /// @param propagation state of the propagation
  template <typename propagator_state_t>
  DETRAY_HOST_DEVICE inline void operator()(
      const state & /*pol_state*/, propagator_state_t &propagation) const {
    propagation.navigation().set_no_trust();
  }
};

/// During guided navigation only the next surface should be re-evaluated. This
/// maps to the 'high trust' level in the navigator
struct guided_navigation : public base_actor {
  struct state {};

  /// Sets the navigation trust level to 'no trust'
  ///
  /// @param pol_state not used
  /// @param propagation state of the propagation
  template <typename propagator_state_t>
  DETRAY_HOST_DEVICE inline void operator()(
      const state & /*pol_state*/, propagator_state_t &propagation) const {
    propagation.navigation().set_high_trust();
  }
};

/// Default navigation update policy for the steppers: If a constraint has been
/// hit, lower the trustlevel to 'fair trust', otherwise stay in 'high trust'.
/// The reasoning is, that the track state might have changed much when a
/// constraint was triggered.
template <concepts::scalar scalar_t>
struct stepper_default_policy : public base_actor {
  struct state {
    scalar_t tol{std::numeric_limits<scalar_t>::epsilon()};
  };

  /// Sets the navigation trust level depending on the step size limit
  ///
  /// @param pol_state not used
  /// @param propagation state of the propagation
  template <typename propagator_state_t>
  DETRAY_HOST_DEVICE inline void operator()(
      state &pol_state, propagator_state_t &propagation) const {
    const auto &stepping = propagation.stepping();
    auto &navigation = propagation.navigation();

    // Not a severe change to track state expected
    if (math::fabs(stepping.step_size()) <
        math::fabs(
            stepping.constraints().template size<>(stepping.direction())) -
            pol_state.tol) {
      // Re-evaluate only next candidate
      navigation.set_high_trust();
    }
    // Step size hit a constraint - the track state was probably changed a
    // lot
    else {
      // Re-evaluate all candidates
      navigation.set_fair_trust();
    }
  }
};

/// Specific navigation policy for the Runge-Kutta stepper: Use the relative
/// amount of step size correction as a measure for the change in direction of
/// the track state.
template <concepts::scalar scalar_t>
struct stepper_rk_policy : public base_actor {
  struct state {
    scalar_t m_threshold_fair_trust{0.05f};
    scalar_t m_threshold_no_trust{0.1f};
  };

  /// Sets the navigation trust level depending on the step size correction
  ///
  /// @param pol_state contains the thresholds
  /// @param propagation state of the propagation
  template <typename propagator_state_t>
  DETRAY_HOST_DEVICE inline void operator()(
      const state &pol_state, propagator_state_t &propagation) const {
    const auto &stepping = propagation.stepping();
    auto &navigation = propagation.navigation();

    // In case of an overlap, have navigator re-evaluate the next candidate
    if (math::fabs(navigation()) <= 1.f * unit<float>::um) {
      navigation.set_high_trust();
      return;
    }

    // How strongly did the RKN algorithm reduce the step size?
    const scalar_t rel_correction{(stepping.step_size() - navigation()) /
                                  navigation()};

    // Large correction to the stepsize - re-initialize the volume
    if (rel_correction > pol_state.m_threshold_no_trust) {
      navigation.set_no_trust();
    }
    // Medium correction - re-evaluate all current candidates
    else if (rel_correction > pol_state.m_threshold_fair_trust) {
      navigation.set_fair_trust();
    } else {
      // Small correction - re-evaluate only the next candidate
      navigation.set_high_trust();
    }
  }
};

}  // namespace detray
