// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/propagator/actors/parameter_updater.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/tracks/bound_track_parameters.hpp"
#include "detray/tracks/free_track_parameters.hpp"

// Test include(s)
#include "detray/test/utils/data_record.hpp"

namespace detray {

namespace detail {
/// Data for a single step
template <concepts::algebra algebra_t>
struct step_data {
  using scalar_type = dscalar<algebra_t>;
  using vector3_type = dvector3D<algebra_t>;
  using track_param_type = free_track_parameters<algebra_t>;
  using bound_param_type = bound_track_parameters<algebra_t>;
  using free_matrix_type = free_matrix<algebra_t>;

  scalar_type step_size{0.f};
  scalar_type path_length{0.f};
  std::size_t n_total_trials{0u};
  navigation::direction nav_dir = navigation::direction::e_forward;
  geometry::identifier identifier{};
  track_param_type track_params{};
  bound_param_type bound_params{};
  free_matrix_type jacobian{};
};
}  // namespace detail

namespace actor {

/// Collect information at every step
template <concepts::algebra algebra_t, template <typename...> class vector_t>
struct step_tracer : public base_actor {
  using step_data_t = detail::step_data<algebra_t>;

  /// Actor state that collects the data
  struct state {
    friend struct step_tracer;

    state() = delete;

    /// Construct the vector containers with a given resource
    /// @param resource
    DETRAY_HOST
    explicit state(vecmem::memory_resource& resource) : m_steps(&resource) {}

    /// Construct from externally provided vector for the @param steps
    DETRAY_HOST_DEVICE
    explicit state(vector_t<step_data_t>&& steps) : m_steps(std::move(steps)) {}

    /// Access to the recorded step data of every step along the track -
    /// const
    DETRAY_HOST_DEVICE
    const auto& get_step_data() const { return m_steps; }

    /// Move the recorded step data out of the actor
    DETRAY_HOST
    auto&& release_step_data() && { return std::move(m_steps); }

   private:
    /// The collected data for the steps
    vector_t<step_data_t> m_steps;
  };

  /// Collect data at every step
  /// @note Primary actor call
  template <typename propagator_state_t>
  DETRAY_HOST_DEVICE void operator()(state& tracer_state,
                                     propagator_state_t& prop_state) const {
    tracer_state.m_steps.push_back(collect_data(prop_state));
  }

  /// Collect only when transport to bound track parameters happens
  /// @note Observer to the parameter updater
  template <typename propagator_state_t>
  DETRAY_HOST_DEVICE void operator()(
      state& tracer_state, propagator_state_t& prop_state,
      const parameter_transporter_result<algebra_t>& res) const {
    const auto& navigation = prop_state.navigation();
    assert(navigation.is_on_surface());

    // If the state has already been collected by the other call operator,
    // update the bound track parameters
    if (!tracer_state.m_steps.empty() &&
        navigation.geometry_identifier() ==
            tracer_state.m_steps.back().identifier) {
      tracer_state.m_steps.back().bound_params = res.destination_params();
    } else {
      tracer_state.m_steps.push_back(
          collect_data(prop_state, res.destination_params()));
    }
  }

 private:
  /// Collect step data
  template <typename propagator_state_t>
  DETRAY_HOST_DEVICE step_data_t collect_data(
      propagator_state_t& prop_state,
      const bound_track_parameters<algebra_t>& bound_param = {}) const {
    const auto& navigation = prop_state.navigation();
    const auto& stepping = prop_state.stepping();

    const auto geo_id{navigation.is_on_surface()
                          ? navigation.geometry_identifier()
                          : geometry::identifier{}};

    return {stepping.step_size(),
            stepping.path_length(),
            stepping.n_total_trials(),
            navigation.direction(),
            geo_id,
            stepping(),
            bound_param,
            stepping.transport_jacobian()};
  }
};

}  // namespace actor

}  // namespace detray
