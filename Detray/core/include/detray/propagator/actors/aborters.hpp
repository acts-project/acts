// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/propagator/base_stepper.hpp"
#include "detray/utils/logging.hpp"

// System include(s)
#include <limits>

namespace detray::actor {

/// Aborter that checks whether the track has exceeded its pathlimit
template <concepts::scalar scalar_t>
struct pathlimit_aborter : public base_actor {
  /// Pathlimit for a single propagation workflow
  struct state {
    /// Absolute path limit
    scalar_t _path_limit = std::numeric_limits<scalar_t>::max();

    /// Set the path limit to a scalar @param pl
    DETRAY_HOST_DEVICE
    inline void set_path_limit(const scalar_t pl) { _path_limit = pl; }

    /// @returns this states remaining path length.
    DETRAY_HOST_DEVICE
    inline scalar_t path_limit() const { return _path_limit; }
  };

  /// Enforces the path limit on a stepper state
  ///
  /// @param abrt_state contains the path limit
  /// @param prop_state state of the propagation
  template <typename propagator_state_t>
  DETRAY_HOST_DEVICE void operator()(state &abrt_state,
                                     propagator_state_t &prop_state) const {
    DETRAY_VERBOSE_HOST_DEVICE("Aborter: Check path limits");

    auto &step_state = prop_state.stepping();
    auto &nav_state = prop_state.navigation();

    // Nothing left to do. Propagation will exit successfully
    if (nav_state.finished()) {
      return;
    }

    const scalar_t step_limit =
        abrt_state.path_limit() -
        math::fabs(prop_state.stepping().abs_path_length());

    // Check the path limit
    if (step_limit <= 0.f) {
      DETRAY_VERBOSE_HOST_DEVICE(
          "Path lengths: %f mm",
          math::fabs(prop_state.stepping().abs_path_length()));

      // Stop navigation
      nav_state.abort("Aborter: Maximal path length reached");
      prop_state.heartbeat(false);
    }

    // Don't go over the path limit in the next step
    step_state.template set_constraint<step::constraint::e_aborter>(step_limit);
  }
};

/// Aborter that checks whether the track fell below a minimum momentum
template <concepts::scalar scalar_t>
struct momentum_aborter : public base_actor {
  struct state {
    /// @returns the momentum limit.
    DETRAY_HOST_DEVICE
    constexpr scalar_t p_limit() const { return m_min_p; }

    /// @returns the momentum limit.
    DETRAY_HOST_DEVICE
    constexpr scalar_t pT_limit() const { return m_min_pT; }

    /// Set the momentum limit to @param p
    DETRAY_HOST_DEVICE
    constexpr void min_p(const scalar_t p) { m_min_p = p; }

    /// Set the transverse momentum limit to @param p
    DETRAY_HOST_DEVICE
    constexpr void min_pT(const scalar_t pT) { m_min_pT = pT; }

   private:
    /// Absolute momentum magnitude limit
    scalar_t m_min_p = 10.f * unit<scalar_t>::MeV;
    scalar_t m_min_pT = 10.f * unit<scalar_t>::MeV;
  };

  /// Actor interface: Enforces a minimum momentum magnitude
  ///
  /// @param abrt_state contains the momentum limit
  /// @param prop_state state of the propagation
  template <typename propagator_state_t>
  DETRAY_HOST_DEVICE void operator()(state &abrt_state,
                                     propagator_state_t &prop_state) const {
    DETRAY_VERBOSE_HOST_DEVICE("Aborter: Check momentum");

    auto &step_state = prop_state.stepping();
    auto &nav_state = prop_state.navigation();

    // Nothing left to do. Propagation will exit successfully
    if (nav_state.finished()) {
      return;
    }

    const auto &track = step_state();
    const scalar_t q{step_state.particle_hypothesis().charge()};

    if (track.pT(q) <= abrt_state.pT_limit()) {
      DETRAY_VERBOSE_HOST_DEVICE("Track |pT| = %f MeV", track.pT(q));

      // Stop navigation
      nav_state.abort("Aborter: Minimum transverse momentum (pT) reached");
      prop_state.heartbeat(false);
      return;
    }

    if (track.p(q) <= abrt_state.p_limit()) {
      DETRAY_VERBOSE_HOST_DEVICE("Track |p| = %f MeV", track.p(q));

      // Stop navigation
      nav_state.abort("Aborter: Minimum momentum (p) reached");
      prop_state.heartbeat(false);
    }
  }
};

/// Aborter checks whether a specific surface was reached
struct target_aborter : public base_actor {
  /// Keeps the index for the target surface
  struct state {
    /// Unique surface id of the target
    geometry::identifier _target_surface;
  };

  /// Exits the navigation as soon as the target surface has been found.
  ///
  /// @param abrt_state contains the target surface index
  /// @param prop_state state of the propagation
  template <typename propagator_state_t>
  DETRAY_HOST_DEVICE void operator()(const state &abrt_state,
                                     propagator_state_t &prop_state) const {
    DETRAY_VERBOSE_HOST_DEVICE("Aborter: Check target surface");

    auto &navigation = prop_state.navigation();
    const auto &stepping = prop_state.stepping();

    // In case the propagation starts on a module, make sure to not abort
    // directly
    if (navigation.is_on_surface() &&
        (navigation.geometry_identifier() == abrt_state._target_surface) &&
        (stepping.path_length() > 0.f)) {
      navigation.abort("Aborter: Reached target surface");
      prop_state.heartbeat(false);
    }
  }
};

}  // namespace detray::actor
