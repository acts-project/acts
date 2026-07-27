/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// TODO: Unexpectedly appeared in device builds for the propagator in l.256
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic warning "-Wmaybe-uninitialized"
#endif

// Project include(s).
#include "kalman_actor.hpp"
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/track_collection.hpp"
#include "traccc/edm/track_parameters.hpp"
#include "traccc/edm/track_state_collection.hpp"
#include "traccc/finding/measurement_selector.hpp"
#include "traccc/fitting/fitting_config.hpp"
#include "traccc/fitting/kalman_filter/kalman_actor.hpp"
#include "traccc/fitting/kalman_filter/kalman_step_aborter.hpp"
#include "traccc/fitting/kalman_filter/statistics_updater.hpp"
#include "traccc/fitting/kalman_filter/two_filters_smoother.hpp"
#include "traccc/fitting/status_codes.hpp"
#include "traccc/utils/logging.hpp"
#include "traccc/utils/particle.hpp"
#include "traccc/utils/prob.hpp"
#include "traccc/utils/propagation.hpp"

// detray include(s)
#include <detray/definitions/navigation.hpp>

// vecmem include(s)
#include <type_traits>

#include <vecmem/containers/device_vector.hpp>

// System include(s).
#include <limits>

namespace traccc {

/// Kalman fitter algorithm to fit a single track
template <typename stepper_t, typename navigator_t>
class kalman_fitter {
 public:
  // Detector type
  using detector_type = typename navigator_t::detector_type;
  using surface_type = typename detector_type::surface_type;
  using algebra_type = typename detector_type::algebra_type;
  using scalar_type = detray::dscalar<algebra_type>;

  /// Configuration type
  using config_type = fitting_config;

  // Field type
  using bfield_type = typename stepper_t::magnetic_field_type;

  // Actor types
  using pathlimit_aborter = detray::actor::pathlimit_aborter<scalar_type>;
  using momentum_aborter = detray::actor::momentum_aborter<scalar_type>;
  using interactor = detray::actor::pointwise_material_interactor<algebra_type>;
  using forward_fit_actor =
      traccc::kalman_actor<algebra_type, surface_type,
                           kalman_actor_direction::FORWARD_ONLY>;
  using backward_fit_actor =
      traccc::kalman_actor<algebra_type, surface_type,
                           kalman_actor_direction::BACKWARD_ONLY>;
  using surface_sequencer = detray::actor::surface_sequencer<surface_type>;

  using forward_updater =
      detray::actor::parameter_updater<algebra_type, interactor,
                                       forward_fit_actor>;
  using backward_updater =
      detray::actor::parameter_updater<algebra_type, backward_fit_actor,
                                       interactor>;

  static_assert(std::is_same_v<typename forward_fit_actor::state,
                               typename backward_fit_actor::state>);

  using forward_actor_chain_type =
      detray::actor_chain<momentum_aborter, pathlimit_aborter, forward_updater,
                          kalman_step_aborter>;

  using backward_actor_chain_type =
      detray::actor_chain<momentum_aborter, pathlimit_aborter, backward_updater,
                          kalman_step_aborter>;

  // Navigator type for backward propagator
  using direct_navigator_type = detray::direct_navigator<detector_type>;

  // Propagator type
  using forward_propagator_type =
      detray::propagator<stepper_t, navigator_t, forward_actor_chain_type>;

  using backward_propagator_type =
      detray::propagator<stepper_t, direct_navigator_type,
                         backward_actor_chain_type>;

  /// Constructor with a detector
  ///
  /// @param det the detector object
  TRACCC_HOST_DEVICE
  kalman_fitter(const detector_type& det, const bfield_type& field,
                const config_type& cfg)
      : m_detector(det), m_field(field), m_cfg(cfg) {}

  /// Kalman fitter state
  struct state {
    /// State constructor
    ///
    /// @param track_states the vector of track states
    TRACCC_HOST_DEVICE
    explicit state(
        const typename edm::track_collection<algebra_type>::device::proxy_type&
            track,
        const typename edm::track_state_collection<algebra_type>::device&
            track_states,
        const edm::measurement_collection::const_device& measurements,
        vecmem::data::vector_view<surface_type> sequence_buffer,
        const detray::propagation::config& prop_cfg,
        const measurement_selector::config& calib_cfg)
        : m_updater_state{prop_cfg},
          m_fit_actor_state{
              track, track_states, measurements,
              vecmem::device_vector<surface_type>(sequence_buffer), calib_cfg},
          m_fit_res{track},
          m_sequence_buffer{sequence_buffer} {}

    /// @return the actor chain state
    TRACCC_HOST_DEVICE
    typename forward_actor_chain_type::state_ref_tuple operator()() {
      return detray::tie(m_step_aborter_state, m_momentum_aborter_state,
                         m_pathlimit_aborter_state, m_interactor_state,
                         m_fit_actor_state, m_updater_state);
    }

    /// @return the actor chain state
    TRACCC_HOST_DEVICE
    typename backward_actor_chain_type::state_ref_tuple backward_actor_state() {
      return detray::tie(m_step_aborter_state, m_momentum_aborter_state,
                         m_pathlimit_aborter_state, m_fit_actor_state,
                         m_interactor_state, m_updater_state);
    }

    /// Individual actor states
    typename pathlimit_aborter::state m_pathlimit_aborter_state{};
    typename momentum_aborter::state m_momentum_aborter_state{};
    typename detray::actor::parameter_updater_state<algebra_type>
        m_updater_state{};
    typename interactor::state m_interactor_state{};
    typename forward_fit_actor::state m_fit_actor_state;
    kalman_step_aborter::state m_step_aborter_state{};

    /// Fitting result per track
    typename edm::track_collection<algebra_type>::device::proxy_type m_fit_res;

    /// View object for identifier sequence
    vecmem::data::vector_view<surface_type> m_sequence_buffer;
  };

  /// Run the kalman fitter for a given number of iterations
  ///
  /// @tparam seed_parameters_t the type of seed track parameter
  ///
  /// @param seed_params seed track parameter
  /// @param fitter_state the state of kalman fitter
  template <typename seed_parameters_t>
  [[nodiscard]] TRACCC_HOST_DEVICE kalman_fitter_status
  fit(const seed_parameters_t& seed_params, state& fitter_state) const {
    seed_parameters_t params = seed_params;
    fitter_state.m_fit_actor_state.reset();
    fitter_state.m_fit_actor_state.do_precise_hole_count =
        m_cfg.do_precise_hole_count;

    // Run the kalman filtering for a given number of iterations
    for (std::size_t i = 0; i < m_cfg.n_iterations; i++) {
      if (kalman_fitter_status res = fit_iteration(params, fitter_state);
          res != kalman_fitter_status::SUCCESS) {
        return res;
      }

      // TODO: For multiple iterations, seed parameter should be set to
      // the first track state which has either filtered or smoothed
      // state. If the first track state is a hole, we need to back
      // extrapolate from the filtered or smoothed state of next valid
      // track state.
      params =
          fitter_state.m_fit_actor_state.m_track_states
              .at(fitter_state.m_fit_actor_state.m_track.constituent_links()
                      .at(0)
                      .index)
              .filtered_params();
      // Reset the iterator of kalman actor
      fitter_state.m_updater_state.init(params);
      fitter_state.m_fit_actor_state.reset();
    }

    return kalman_fitter_status::SUCCESS;
  }

  template <typename seed_parameters_t>
  [[nodiscard]] TRACCC_HOST_DEVICE kalman_fitter_status
  fit_iteration(seed_parameters_t params, state& fitter_state) const {
    inflate_covariance(params, m_cfg.covariance_inflation_factor);

    const kalman_fitter_status res_fw = filter(params, fitter_state);

    // Reset hole count
    const unsigned int n_holes_fw{fitter_state.m_fit_actor_state.n_holes};
    fitter_state.m_fit_actor_state.n_holes = 0u;

    // Run smoothing
    kalman_fitter_status res_bw{kalman_fitter_status::ERROR_OTHER};
    if (res_fw == kalman_fitter_status::SUCCESS) {
      res_bw = smooth(fitter_state);

      // In case the smoother does not apply hole counting
      fitter_state.m_fit_actor_state.n_holes =
          math::max(fitter_state.m_fit_actor_state.n_holes, n_holes_fw);
    }

    TRACCC_VERBOSE_HOST_DEVICE("Number of holes after forward fit: %d",
                               n_holes_fw);
    TRACCC_VERBOSE_HOST_DEVICE("Updated number of holes after smoothing: %d",
                               fitter_state.m_fit_actor_state.n_holes);

    // Update track fitting qualities
    update_statistics(fitter_state);

    // Check the fit result
    check_fitting_result(fitter_state, res_fw, res_bw);

    return res_bw;
  }

  /// Run the kalman fitter for an iteration
  ///
  /// @tparam seed_parameters_t the type of seed track parameter
  ///
  /// @param seed_params seed track parameter
  /// @param fitter_state the state of kalman fitter
  template <typename seed_parameters_t>
  [[nodiscard]] TRACCC_HOST_DEVICE kalman_fitter_status
  filter(const seed_parameters_t& seed_params, state& fitter_state) const {
    TRACCC_VERBOSE_HOST_DEVICE("Run forward fit...");

    // Create propagator
    forward_propagator_type propagator(m_cfg.propagation);

    // Set initial track parameters for parameter transport
    fitter_state.m_updater_state.init(seed_params);

    // Set minimum momentum
    fitter_state.m_momentum_aborter_state.min_pT(
        static_cast<scalar_type>(m_cfg.min_pT));
    fitter_state.m_momentum_aborter_state.min_p(
        static_cast<scalar_type>(m_cfg.min_p));

    // Set path limit
    fitter_state.m_pathlimit_aborter_state.set_path_limit(
        m_cfg.propagation.stepping.path_limit);

    // Create propagator state
    typename forward_propagator_type::state propagation(
        seed_params, m_field, m_detector, m_cfg.propagation.context);

    propagation.set_particle(
        detail::correct_particle_hypothesis(m_cfg.ptc_hypothesis, seed_params));

    // Set overstep tolerance, stepper constraint and mask tolerance
    propagation.stepping()
        .template set_constraint<detray::step::constraint::e_accuracy>(
            m_cfg.propagation.stepping.step_constraint);

    // Reset fitter statistics
    fitter_state.m_fit_res.reset_quality();

    // Run forward filtering
    propagator.propagate(propagation, fitter_state());

    // Encountered error during fitting?
    if (fitter_state.m_fit_actor_state.fit_result !=
        kalman_fitter_status::SUCCESS) {
      return fitter_state.m_fit_actor_state.fit_result;
    }
    // Encountered error during propagation?
    if (!propagator.finished(propagation)) {
      return kalman_fitter_status::ERROR_PROPAGATION_FAILURE;
    }

    return kalman_fitter_status::SUCCESS;
  }

  /// Run smoothing after kalman filtering
  ///
  /// @brief The smoother is based on "Application of Kalman filtering to
  /// track and vertex fitting", R.Frühwirth, NIM A.
  ///
  /// @param fitter_state the state of kalman fitter
  [[nodiscard]] TRACCC_HOST_DEVICE kalman_fitter_status
  smooth(state& fitter_state) const {
    TRACCC_VERBOSE_HOST_DEVICE("Run smoothing...");

    if (fitter_state.m_fit_actor_state.sequencer().overflow()) {
      TRACCC_ERROR_HOST_DEVICE("Surface sequence overlow");
      return kalman_fitter_status::ERROR_GEOID_SEQUENCE_OVERFLOW;
    }
    if (fitter_state.m_fit_actor_state.sequencer().sequence().empty()) {
      TRACCC_ERROR_HOST_DEVICE("Surface sequence empty");
      return kalman_fitter_status::ERROR_UPDATER_SKIPPED_STATE;
    }

    // Since the smoothed track parameter of the last surface can be
    // considered to be the filtered one, we can reversly iterate the
    // algorithm to obtain the smoothed parameter of other surfaces
    fitter_state.m_fit_actor_state.reset();

    constexpr bool backward_mode{true};
    while (!fitter_state.m_fit_actor_state.finished() &&
           fitter_state.m_fit_actor_state(backward_mode).is_hole()) {
      fitter_state.m_fit_actor_state.next();
    }

    if (fitter_state.m_fit_actor_state.finished()) {
      return kalman_fitter_status::ERROR_UPDATER_SKIPPED_STATE;
    }

    auto last = fitter_state.m_fit_actor_state(backward_mode);

    const scalar theta = last.filtered_params().theta();
    if (!std::isfinite(theta)) {
      TRACCC_ERROR_HOST_DEVICE(
          "Theta is infinite after forward fit (Matrix inversion)");
      return kalman_fitter_status::ERROR_INVERSION;
    }

    if (!std::isfinite(last.filtered_params().phi())) {
      TRACCC_ERROR_HOST_DEVICE(
          "Phi is infinite after forward fit (Matrix inversion)");
      return kalman_fitter_status::ERROR_INVERSION;
    }

    if (theta <= 0.f || theta >= 2.f * constant<traccc::scalar>::pi) {
      TRACCC_ERROR_HOST_DEVICE("Hit theta pole after forward fit : %f", theta);
      TRACCC_ERROR_HOST("Params: " << last.filtered_params());
      return kalman_fitter_status::ERROR_THETA_POLE;
    }

    last.smoothed_params() = last.filtered_params();

    TRACCC_DEBUG_HOST("Start smoothing at: " << last.smoothed_params());

    // Configure actors
    fitter_state.m_updater_state.init(last.smoothed_params());
    fitter_state.m_updater_state.noise_estimation_cfg()
        .estimate_scattering_noise = false;

    backward_propagator_type propagator(m_cfg.propagation);

    typename backward_propagator_type::state propagation(
        last.smoothed_params(), m_field, m_detector,
        fitter_state.m_sequence_buffer, m_cfg.propagation.context);

    // Configure backward propagation state
    propagation.navigation().set_direction(
        detray::navigation::direction::e_backward);
    propagation.navigation().reset();

    // Synchronize the current geo ID with the input track parameter
    TRACCC_DEBUG_HOST("Expecting: " << last.smoothed_params().surface_link());
    while (propagation.navigation().has_next_external() &&
           propagation.navigation().next_external().identifier() !=
               last.smoothed_params().surface_link()) {
      TRACCC_DEBUG_HOST(
          "Advancing to next external surface from: "
          << propagation.navigation().next_external().identifier());
      propagation.navigation().advance();
    }

    // No valid states produced by forward pass
    if (propagation.navigation().finished()) {
      TRACCC_ERROR_HOST_DEVICE("No matching surfaces!");
      return kalman_fitter_status::ERROR_UPDATER_SKIPPED_STATE;
    }

    // Prepare the parameter transport
    propagation.set_particle(detail::correct_particle_hypothesis(
        m_cfg.ptc_hypothesis, last.smoothed_params()));

    assert(
        std::signbit(propagation.stepping().particle_hypothesis().charge()) ==
        std::signbit(fitter_state.m_updater_state.bound_params().qop()));

    inflate_covariance(fitter_state.m_updater_state.bound_params(),
                       m_cfg.covariance_inflation_factor);

    // Run the smoothing
    propagator.propagate(propagation, fitter_state.backward_actor_state());

    // Encountered error in the smoother?
    if (fitter_state.m_fit_actor_state.fit_result !=
        kalman_fitter_status::SUCCESS) {
      TRACCC_ERROR_HOST_DEVICE("Fit failed!");
      return fitter_state.m_fit_actor_state.fit_result;
    }
    // Encountered error during propagation?
    if (!propagator.finished(propagation)) {
      TRACCC_ERROR_HOST_DEVICE("-> Propagation failed!");
      return kalman_fitter_status::ERROR_PROPAGATION_FAILURE;
    }

    return kalman_fitter_status::SUCCESS;
  }

  TRACCC_HOST_DEVICE
  void update_statistics(state& fitter_state) const {
    typename edm::track_collection<algebra_type>::device::proxy_type& fit_res =
        fitter_state.m_fit_res;
    typename edm::track_state_collection<algebra_type>::device& track_states =
        fitter_state.m_fit_actor_state.m_track_states;

    // Only count statistics for fully fitted tracks
    if (fitter_state.m_fit_actor_state.fit_result !=
        kalman_fitter_status::SUCCESS) {
      return;
    }

    TRACCC_DEBUG_HOST_DEVICE("Updating fit statistics");

    // Fit parameter = smoothed track parameter of the first smoothed track
    // state
    for (const edm::track_constituent_link& link :
         fit_res.constituent_links()) {
      assert(link.type == edm::track_constituent_link::track_state);
      if (track_states.at(link.index).is_smoothed()) {
        fit_res.params() = track_states.at(link.index).smoothed_params();
        break;
      }
    }

    for (const edm::track_constituent_link& link :
         fit_res.constituent_links()) {
      assert(link.type == edm::track_constituent_link::track_state);
      edm::track_state trk_state = track_states.at(link.index);
      const detray::tracking_surface sf{
          m_detector, fitter_state.m_fit_actor_state.m_measurements
                          .at(trk_state.measurement_index())
                          .surface_link()};
      statistics_updater<algebra_type>{}(
          fit_res, trk_state, fitter_state.m_fit_actor_state.m_measurements);
    }

    // Subtract the NDoF with the degree of freedom of the bound track (=5)
    fit_res.ndf() -= 5.f;
    fit_res.pval() = prob(fit_res.chi2(), fit_res.ndf());

    // The number of holes
    fit_res.nholes() = fitter_state.m_fit_actor_state.n_holes;
  }

  TRACCC_HOST_DEVICE
  void check_fitting_result(state& fitter_state,
                            const kalman_fitter_status res_fw,
                            const kalman_fitter_status res_bw) const {
    typename edm::track_collection<algebra_type>::device::proxy_type& fit_res =
        fitter_state.m_fit_res;
    const typename edm::track_state_collection<algebra_type>::device&
        track_states = fitter_state.m_fit_actor_state.m_track_states;

    TRACCC_DEBUG_HOST_DEVICE("Checking fit results...");

    // Check the fitter result code
    if (res_fw != kalman_fitter_status::SUCCESS) {
      // Check if track states were skipped
      if (res_fw == kalman_fitter_status::ERROR_UPDATER_SKIPPED_STATE) {
        TRACCC_ERROR_HOST_DEVICE(
            "Skipped sensitive surfaces during forward fit: %d",
            fitter_state.m_fit_actor_state.count_missed_fit());
        fit_res.fit_outcome() = track_fit_outcome::FAILURE_NOT_ALL_FITTED;
        return;
      }

      // Check for propgation failure, otherwise assume fitter error
      TRACCC_ERROR_HOST("Error during fitting: " << fitter_debug_msg{res_fw}());
      fit_res.fit_outcome() =
          res_fw == kalman_fitter_status::ERROR_PROPAGATION_FAILURE
              ? track_fit_outcome::FAILURE_FORWARD_PROPAGATION
              : track_fit_outcome::FAILURE_FITTER;
      return;
    }

    if (res_bw != kalman_fitter_status::SUCCESS) {
      // Check if track states were skipped
      if (res_bw == kalman_fitter_status::ERROR_SMOOTHER_SKIPPED_STATE) {
        TRACCC_ERROR_HOST_DEVICE(
            "Skipped sensitive surfaces during smoothing: %d",
            fitter_state.m_fit_actor_state.count_missed_smoother());
        fit_res.fit_outcome() = track_fit_outcome::FAILURE_NOT_ALL_SMOOTHED;
        return;
      }

      // Check for propgation failure, otherwise assume smoother error
      TRACCC_ERROR_HOST(
          "Error during smoothing: " << fitter_debug_msg{res_bw}());
      fit_res.fit_outcome() =
          res_bw == kalman_fitter_status::ERROR_PROPAGATION_FAILURE
              ? track_fit_outcome::FAILURE_BACKWARD_PROPAGATION
              : track_fit_outcome::FAILURE_SMOOTHER;
      return;
    }

    if (const unsigned int n_holes{fitter_state.m_fit_actor_state.n_holes};
        n_holes > 2u) {
      TRACCC_WARNING_HOST_DEVICE("Exceeded max. hole count: %d", n_holes);
    }

    // NDF should always be positive for fitting
    if (fit_res.ndf() > 0) {
      for (const edm::track_constituent_link& link :
           fit_res.constituent_links()) {
        assert(link.type == edm::track_constituent_link::track_state);
        edm::track_state trk_state = track_states.at(link.index);

        // Fitting fails if any of non-hole track states is not smoothed
        if (!trk_state.is_hole() && !trk_state.is_smoothed()) {
          TRACCC_ERROR_HOST_DEVICE("Not all smoothed");
          fit_res.fit_outcome() = track_fit_outcome::FAILURE_NOT_ALL_SMOOTHED;
          return;
        }
      }

      // Fitting succeeds if any of non-hole track states is not smoothed
      TRACCC_DEBUG_HOST_DEVICE("Fit status: SUCCESS");
      fit_res.fit_outcome() = track_fit_outcome::SUCCESS;
      return;
    }

    TRACCC_ERROR_HOST_DEVICE("Negative NDF");
    fit_res.fit_outcome() = track_fit_outcome::FAILURE_NON_POSITIVE_NDF;
    return;
  }

  TRACCC_HOST_DEVICE
  const config_type& config() const { return m_cfg; }

 private:
  // Detector object
  const detector_type& m_detector;
  // Field object
  const bfield_type m_field;

  // Configuration object
  config_type m_cfg;
};

}  // namespace traccc
