// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/definitions/detail/macros.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/base_stepper.hpp"
#include "detray/propagator/concepts.hpp"
#include "detray/propagator/detail/noise_estimation.hpp"
#include "detray/propagator/propagation_config.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/utils/logging.hpp"

// System include(s).
#include <iomanip>

namespace detray {

/// Templated propagator class, using a stepper and a navigator object in
/// succession.
///
/// @tparam stepper_t for the transport
/// @tparam navigator_t for the navigation
template <typename stepper_t, typename navigator_t,
          concepts::actor_chain actor_chain_t>
struct propagator {
  /// Stepper implementation used to advance the track.
  using stepper_type = stepper_t;
  /// Navigator implementation used to find geometry intersections.
  using navigator_type = navigator_t;
  /// Chain of actors invoked during propagation.
  using actor_chain_type = actor_chain_t;

  /// Detector described by the navigator.
  using detector_type = typename navigator_type::detector_type;
  /// Algebra backend used by the detector.
  using algebra_type = typename detector_type::algebra_type;
  /// Scalar type associated with the detector algebra.
  using scalar_type = dscalar<algebra_type>;
  /// Intersection representation produced by the navigator.
  using intersection_type = typename navigator_type::intersection_type;
  /// Free track parameter representation advanced by the stepper.
  using free_track_parameters_type =
      typename stepper_t::free_track_parameters_type;
  /// Bound track parameter representation accepted by the stepper.
  using bound_track_parameters_type =
      typename stepper_t::bound_track_parameters_type;

  /// Propagation and stepping configuration; lifetime must outlive this object.
  const propagation::config &m_cfg;

  /// Track integration algorithm.
  stepper_t m_stepper{};
  /// Geometry navigation algorithm.
  navigator_t m_navigator{};

  /// Register the actor types
  const actor_chain_type run_actors{};

  /// Construct from a propagator configuration.
  /// @param cfg Propagation configuration, retained by reference.
  DETRAY_HOST_DEVICE
  explicit constexpr propagator(const propagation::config &cfg) : m_cfg{cfg} {}

  /// Propagation that state aggregates a stepping and a navigation state.
  /// It also keeps references to the actor states.
  template <bool states_as_reference = false>
  struct state_base {
    /// Detector associated with the navigator.
    using detector_type = typename navigator_t::detector_type;
    /// Geometry context type for the detector.
    using context_type = typename detector_type::geometry_context;
    /// Detector algebra backend.
    using algebra_type = typename detector_type::algebra_type;
    /// Scalar type used by propagation state.
    using scalar_type = typename detector_type::scalar_type;

    /// Stepper-owned state type.
    using stepper_state_type = typename stepper_t::state;
    /// Navigator-owned state type.
    using navigator_state_type = typename navigator_t::state;
    /// Actor chain associated with this propagator.
    using actor_chain_type = actor_chain_t;

    /// Whether the configured stepper transports field gradients.
    static constexpr bool stepper_uses_gradient = stepper_type::uses_gradient;

    /// Construct an owning propagation state from free parameters.
    /// @param free_params Initial free track parameters.
    /// @param det Detector used for navigation.
    /// @param ctx Geometry context for detector transforms.
    DETRAY_HOST_DEVICE state_base(const free_track_parameters_type &free_params,
                                  const detector_type &det,
                                  const context_type &ctx)
      requires(!states_as_reference)
        : m_stepping(free_params), m_navigation(det), m_context(ctx) {}

    /// Construct an owning state from free parameters and a magnetic field.
    /// @param free_params Initial free track parameters.
    /// @param magnetic_field Magnetic field used by the stepper.
    /// @param det Detector used for navigation.
    /// @param ctx Geometry context for detector transforms.
    template <typename field_t>
      requires(!states_as_reference &&
               concepts::is_field_of<stepper_type, field_t>)
    DETRAY_HOST_DEVICE state_base(const free_track_parameters_type &free_params,
                                  const field_t &magnetic_field,
                                  const detector_type &det,
                                  const context_type &ctx = {})
        : m_stepping(free_params, magnetic_field),
          m_navigation(det),
          m_context(ctx) {}

    /// Construct an owning state from free parameters and a navigator view.
    /// @param free_params Initial free track parameters.
    /// @param det Detector used for navigation.
    /// @param nav_view Existing navigator state view to initialize from.
    /// @param ctx Geometry context for detector transforms.
    DETRAY_HOST_DEVICE state_base(
        const free_track_parameters_type &free_params, const detector_type &det,
        typename navigator_type::state::view_type nav_view,
        const context_type &ctx = {})
      requires(!states_as_reference)
        : m_stepping(free_params),
          m_navigation(det, nav_view),
          m_context(ctx) {}

    /// Construct an owning state from free parameters, field, and navigator
    /// view.
    /// @param free_params Initial free track parameters.
    /// @param magnetic_field Magnetic field used by the stepper.
    /// @param det Detector used for navigation.
    /// @param nav_view Existing navigator state view to initialize from.
    /// @param ctx Geometry context for detector transforms.
    template <typename field_t>
      requires(!states_as_reference &&
               concepts::is_field_of<stepper_type, field_t>)
    DETRAY_HOST_DEVICE state_base(
        const free_track_parameters_type &free_params,
        const field_t &magnetic_field, const detector_type &det,
        typename navigator_type::state::view_type nav_view,
        const context_type &ctx = {})
        : m_stepping(free_params, magnetic_field),
          m_navigation(det, nav_view),
          m_context(ctx) {}

    /// Construct an owning state from bound parameters.
    /// @param param Initial bound track parameters.
    /// @param det Detector used for navigation.
    /// @param ctx Geometry context for detector transforms.
    DETRAY_HOST_DEVICE state_base(const bound_track_parameters_type &param,
                                  const detector_type &det,
                                  const context_type &ctx = {})
      requires(!states_as_reference)
        : m_stepping(param, det, ctx), m_navigation(det), m_context(ctx) {
      m_navigation.set_volume(param.surface_link().volume());
    }

    /// Construct a non-owning view of existing stepper and navigator states.
    /// @param stepping Existing stepper state.
    /// @param navigation Existing navigator state.
    /// @param ctx Geometry context for detector transforms.
    DETRAY_HOST_DEVICE state_base(stepper_state_type &stepping,
                                  navigator_state_type &navigation,
                                  const context_type &ctx = {})
      requires(states_as_reference)
        : m_stepping(stepping), m_navigation(navigation), m_context(ctx) {}

    /// Construct an owning state from bound parameters and a magnetic field.
    /// @param param Initial bound track parameters.
    /// @param magnetic_field Magnetic field used by the stepper.
    /// @param det Detector used for navigation.
    /// @param ctx Geometry context for detector transforms.
    template <typename field_t>
      requires concepts::is_field_of<stepper_type, field_t>
    DETRAY_HOST_DEVICE state_base(const bound_track_parameters_type &param,
                                  const field_t &magnetic_field,
                                  const detector_type &det,
                                  const context_type &ctx = {})
        : m_stepping(param, magnetic_field, det, ctx),
          m_navigation(det),
          m_context(ctx) {
      m_navigation.set_volume(param.surface_link().volume());
    }

    /// Construct an owning state from bound parameters and navigator view.
    /// @param param Initial bound track parameters.
    /// @param magnetic_field Magnetic field used by the stepper.
    /// @param det Detector used for navigation.
    /// @param nav_view Existing navigator state view to initialize from.
    /// @param ctx Geometry context for detector transforms.
    template <typename field_t>
      requires concepts::is_field_of<stepper_type, field_t>
    DETRAY_HOST_DEVICE state_base(
        const bound_track_parameters_type &param, const field_t &magnetic_field,
        const detector_type &det,
        typename navigator_type::state::view_type nav_view,
        const context_type &ctx = {})
        : m_stepping(param, magnetic_field, det, ctx),
          m_navigation(det, nav_view),
          m_context(ctx) {
      m_navigation.set_volume(param.surface_link().volume());
    }

    /// Set the particle hypothesis used by the stepper.
    /// @param ptc Particle properties, including mass and charge.
    DETRAY_HOST_DEVICE
    void set_particle(const pdg_particle<scalar_type> &ptc) {
      m_stepping.set_particle(ptc);
    }

    /// @returns Whether propagation is currently active.
    DETRAY_HOST_DEVICE
    bool is_alive() const { return m_heartbeat; }

    /// @returns Whether propagation is currently active.
    DETRAY_HOST_DEVICE
    bool heartbeat() const { return m_heartbeat; }

    /// Set whether propagation is currently active.
    /// @param heartbeat True to continue propagation.
    DETRAY_HOST_DEVICE
    void heartbeat(bool heartbeat) { m_heartbeat = heartbeat; }

    /// @returns The stepper state (read-only).
    DETRAY_HOST_DEVICE
    const stepper_state_type &stepping() const { return m_stepping; }

    /// @returns The mutable stepper state.
    DETRAY_HOST_DEVICE
    stepper_state_type &stepping() { return m_stepping; }

    /// @returns The navigator state (read-only).
    DETRAY_HOST_DEVICE
    const navigator_state_type &navigation() const { return m_navigation; }

    /// @returns The mutable navigator state.
    DETRAY_HOST_DEVICE
    navigator_state_type &navigation() { return m_navigation; }

    /// @returns The geometry context used during propagation (read-only).
    DETRAY_HOST_DEVICE
    const context_type &context() const { return m_context; }

    /// @returns The mutable geometry context used during propagation.
    DETRAY_HOST_DEVICE
    context_type &context() { return m_context; }

    /// @returns Whether verbose propagation diagnostics are enabled.
    DETRAY_HOST_DEVICE
    bool debug() const { return m_do_debug; }

    /// Enable or disable verbose propagation diagnostics.
    /// @param b True to enable diagnostics.
    DETRAY_HOST_DEVICE
    void debug(bool b) { m_do_debug = b; }

   private:
    std::conditional_t<states_as_reference, stepper_state_type &,
                       stepper_state_type>
        m_stepping;
    std::conditional_t<states_as_reference, navigator_state_type &,
                       navigator_state_type>
        m_navigation;

    context_type m_context;

    // Is the propagation still alive?
    bool m_heartbeat = false;
    bool m_do_debug = false;
  };

  /// Owning propagation state type.
  using state = state_base<false>;

  /// Propagate method finale: Return whether or not the propagation
  /// completed successfully.
  ///
  /// @param propagation the state of a propagation flow
  ///
  /// @return propagation success.
  template <bool is_owning>
  DETRAY_HOST_DEVICE bool finished(
      const state_base<is_owning> &propagation) const {
    return propagation.navigation().finished();
  }

  /// @returns true if the @param propagation is suspended
  template <bool is_owning>
  DETRAY_HOST_DEVICE inline auto is_paused(
      const state_base<is_owning> &propagation) const -> bool {
    return !propagation.is_alive() && propagation.navigation().is_alive();
  }

  /// Revive a paused propagation.
  /// @param propagation State to resume.
  template <bool is_owning>
  DETRAY_HOST_DEVICE inline void resume(
      state_base<is_owning> &propagation) const {
    assert(propagation.navigation().is_alive());
    propagation.heartbeat(true);
  }

  /// Propagate method: Coordinates the calls of the stepper, navigator
  /// and all registered actors.
  ///
  /// @param propagation the state of a propagation flow
  /// @param actor_state_refs tuple containing references to the actor
  /// states
  ///
  /// @return True if propagation completed or was paused; false if aborted.
  template <typename actor_states_t, bool is_owning>
    requires(concepts::is_state_of<actor_states_t, actor_chain_type>)
  DETRAY_HOST_DEVICE bool propagate(
      state_base<is_owning> &propagation,
      actor_states_t actor_state_refs = dtuple<>{}) const {
    auto &navigation = propagation.navigation();
    auto &stepping = propagation.stepping();
    auto &context = propagation.context();
    const auto &track = stepping();
    assert(!track.is_invalid());

    DETRAY_VERBOSE_HOST("Starting propagation for track:\n" << track);

    bool is_init = false;
    if (this->is_paused(propagation)) {
      DETRAY_VERBOSE_HOST("Resuming propagation...");
    } else {
      // Initialize the navigation
      DETRAY_VERBOSE_HOST("Initialize navigation...");
      m_navigator.init(track, navigation, m_cfg.navigation, context);
      propagation.heartbeat(navigation.is_alive());

      is_init = true;
    }

    // Run while there is a heartbeat. In order to help the compiler
    // optimize this, and in order to make the code more GPU-friendly,
    // this code is run as a flat loop, but this loop has a defined
    // structure. Indeed, the structure is always to run either the
    // actors or the stepper (in alternating order) followed by the
    // navigation update.
    //
    // A = actors
    // N = navigation update
    // S = propagation step
    //
    // ANSNANSNANSNANSNANSNANS...
    scalar_type path_length{0.f};
    unsigned int stall_counter{0u};
    for (unsigned int i = 0; i % 2 == 0 || propagation.is_alive(); ++i) {
      if (i % 2 == 0) {
        DETRAY_VERBOSE_HOST_DEVICE("Propagation step: %d", i / 2);
        DETRAY_VERBOSE_HOST_DEVICE("-> Path length: %f mm",
                                   stepping.path_length());

        // Run all registered actors/aborters
        run_actors(actor_state_refs, propagation);

        // Don't run another navigation update, if already exited
        if (!propagation.is_alive()) {
          continue;
        }

        path_length = stepping.path_length();

        assert(!track.is_invalid());
      } else {
        assert(!track.is_invalid());

        // Set access to the volume material for the stepper
        auto vol = navigation.current_volume();
        const material<scalar_type> *vol_mat_ptr =
            vol.has_material() ? vol.material_parameters(track.pos()) : nullptr;

        // Break automatic step size scaling by the stepper when a
        // surface was reached and whenever the navigation is
        // (re-)initialized
        const bool reset_stepsize{navigation.is_on_surface() || is_init};

        // Take the step
        DETRAY_VERBOSE_HOST("Calling stepper...");
        // GCC reports 'propagator' as maybe-uninitialized here once this call
        // is inlined at -O3. It is a false positive: m_cfg is bound in the
        // constructor's init list and m_stepper/m_navigator have default member
        // initializers, so the object is always fully initialized. The pragma
        // has to sit here, at the diagnostic's own location, rather than in the
        // including header -- that way it applies no matter which header first
        // pulls this file into a translation unit.
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
        propagation.heartbeat(propagation.heartbeat() &&
                              m_stepper.step(navigation(), stepping,
                                             m_cfg.stepping, reset_stepsize,
                                             vol_mat_ptr));
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif

        // Reduce navigation trust level according to stepper update
        DETRAY_VERBOSE_HOST("-> Evaluate stepper navigation policy:");
        typename stepper_t::policy_type{}(stepping.policy_state(), propagation);

        if (i > 0) {
          is_init = false;
        }

        // Check if the propagation makes progress
        if (math::fabs(stepping.path_length()) <=
            math::fabs(path_length) +
                m_cfg.navigation.intersection.path_tolerance) {
          if (stall_counter >= 10u) {
            propagation.heartbeat(false);
            navigation.abort("Propagation stalled");
            DETRAY_ERROR_HOST("Propagation stalled");
          } else if (stall_counter > 2u) {
            // Print a warning if the propagation starts stalling
            // (no overlap)
            DETRAY_WARN_HOST("Propagation is stalling (counter "
                             << stall_counter << ")");
            DETRAY_VERBOSE_DEVICE("Propagation is stalling (counter %d)",
                                  stall_counter);
            DETRAY_WARN_HOST(print(propagation));
            DETRAY_WARN_HOST("-> Track: " << stepping());
          }
          DETRAY_DEBUG_HOST_DEVICE("-> Step stalled. Counter %d",
                                   stall_counter);
          stall_counter++;
        } else {
          stall_counter = 0u;
        }
      }

      // Find next candidate
      DETRAY_VERBOSE_HOST("Calling navigator...");
      const bool nav_is_init =
          m_navigator.update(track, navigation, m_cfg.navigation, context);
      is_init = is_init || nav_is_init;

      propagation.heartbeat(propagation.heartbeat() && navigation.is_alive());

      if (i % 2 == 0 && i > 0 && propagation.debug()) {
        DETRAY_VERBOSE_HOST(print(propagation));
      }
    }

    // Pass on the whether the propagation was successful
    DETRAY_VERBOSE_HOST("Finished propagation for track:\n" << track);
    if (finished(propagation)) {
      DETRAY_VERBOSE_HOST_DEVICE("Status: SUCCESS");
    } else if (is_paused(propagation)) {
      DETRAY_VERBOSE_HOST_DEVICE("Status: PAUSED");
    } else {
      DETRAY_VERBOSE_HOST_DEVICE("Status: ABORT");
    }

    return finished(propagation) || is_paused(propagation);
  }

  /// Propagate using an empty actor chain.
  /// @param propagation State of the propagation flow.
  /// @return True if propagation completed or was paused; false if aborted.
  template <bool is_owning>
  DETRAY_HOST_DEVICE bool propagate(state_base<is_owning> &propagation) {
    // Will not be used
    actor_chain<>::state empty_state{};
    // Run propagation
    return propagate(propagation, empty_state);
  }

  /// Format propagation state for diagnostic output.
  /// @param propagation State to describe.
  /// @return Human-readable status, navigation, and track information.
  template <bool is_owning>
  DETRAY_HOST std::string print(state_base<is_owning> &propagation) const {
    const auto &navigation = propagation.navigation();
    const auto &stepping = propagation.stepping();

    std::stringstream debug_stream{};
    debug_stream << std::left << std::setw(10);
    debug_stream << "\nstatus: " << navigation.status() << std::endl;

    debug_stream << "volume: " << std::setw(10);
    if (detail::is_invalid_value(navigation.volume())) {
      debug_stream << "invalid";
    } else {
      debug_stream << navigation.volume();
    }
    debug_stream << std::endl;

    debug_stream << "navigation:" << std::endl;
    if (navigation.is_on_surface()) {
      debug_stream << std::setw(10)
                   << " -> on surface: " << navigation.geometry_identifier();
    } else {
      debug_stream << std::setw(10) << " -> target: "
                   << navigation.target().surface().identifier();
    }
    debug_stream << std::endl;
    debug_stream << std::setw(10) << " -> path: " << navigation() << " mm"
                 << std::endl;

    debug_stream << "stepping:" << std::endl;
    debug_stream << std::setw(10) << " -> step size: " << stepping.step_size()
                 << " mm" << std::endl;
    debug_stream << " -> " << detail::ray<algebra_type>(stepping());

    return debug_stream.str();
  }
};
}  // namespace detray
