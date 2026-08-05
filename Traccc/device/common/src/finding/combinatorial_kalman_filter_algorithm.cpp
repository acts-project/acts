/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/finding/device/combinatorial_kalman_filter_algorithm.hpp"

// Project include(s).
#include "traccc/finding/candidate_link.hpp"

// System include(s).
#include <stdexcept>

namespace traccc::device {

struct combinatorial_kalman_filter_algorithm::data {
  /// The (finding) algorithm configuration
  finding_config m_config;
};

combinatorial_kalman_filter_algorithm::combinatorial_kalman_filter_algorithm(
    const finding_config& config, const traccc::memory_resource& mr,
    const vecmem::copy& copy, std::unique_ptr<const Logger> logger,
    std::unique_ptr<kalman_fitting_algorithm> kf_fitter)
    : messaging(std::move(logger)),
      algorithm_base{mr, copy},
      m_data{std::make_unique<data>(config)},
      m_kf_fitter{std::move(kf_fitter)} {
  if (config.run_smoother == smoother_type::e_kalman &&
      m_kf_fitter == nullptr) {
    throw std::invalid_argument(
        "A Kalman fitting algorithm must be provided when the Kalman "
        "smoother is enabled");
  }

  if (config.min_step_length_for_next_surface <=
      math::fabs(
          config.propagation.navigation.intersection.overstep_tolerance)) {
    throw std::invalid_argument(
        "Min step length for the next surface should be higher than the "
        "overstep tolerance");
  }
}

combinatorial_kalman_filter_algorithm::
    ~combinatorial_kalman_filter_algorithm() = default;

auto combinatorial_kalman_filter_algorithm::operator()(
    const detector_buffer& det, const magnetic_field& bfield,
    const edm::measurement_collection::const_view& measurements_view,
    const bound_track_parameters_collection_types::const_view& seeds) const
    -> output_type {
  using algebra_t = default_algebra;
  using scalar_t = detray::dscalar<algebra_t>;

  assert(m_data);
  assert(input_is_valid(measurements_view));

  const finding_config& cfg = m_data->m_config;

  /*****************************************************************
   * Measurement Operations
   *****************************************************************/

  // Get the number of measurements. In an asynchronous way if possible.
  edm::measurement_collection::const_view::size_type n_measurements = 0u;
  if (mr().host) {
    vecmem::async_size size = copy().get_size(measurements_view, *(mr().host));
    // Here we could give control back to the caller, once our code allows
    // for it. (coroutines...)
    n_measurements = size.get();
  } else {
    n_measurements = copy().get_size(measurements_view);
  }

  // Get upper bounds of measurement ranges per surface
  vecmem::data::vector_buffer<unsigned int> meas_ranges_buffer =
      build_measurement_ranges_buffer(det, n_measurements, measurements_view);

  // Get the number of track seeds. In an asynchronous way if possible.
  bound_track_parameters_collection_types::const_view::size_type n_seeds = 0u;
  if (mr().host) {
    vecmem::async_size size = copy().get_size(seeds, *(mr().host));
    // Here we could give control back to the caller, once our code allows
    // for it. (coroutines...)
    n_seeds = size.get();
  } else {
    n_seeds = copy().get_size(seeds);
  }

  // Prepare input parameters with seeds
  bound_track_parameters_collection_types::buffer in_params_buffer(n_seeds,
                                                                   mr().main);
  copy().setup(in_params_buffer)->ignore();
  copy()(seeds, in_params_buffer, vecmem::copy::type::device_to_device)
      ->ignore();

  // Create track buffer
  vecmem::vector<unsigned int> empty_vec{};
  typename edm::track_container<default_algebra>::buffer
      track_candidates_buffer{
          {empty_vec, mr().main, mr().host},
          {0u, mr().main, vecmem::data::buffer_type::resizable},
          measurements_view};

  // Prepare the payload for the fitting kernel (only for valid fitting alg.)
  kalman_fitting_algorithm::fit_payload smoothing_payload{det, bfield};

  if (cfg.run_smoother == smoother_type::e_kalman) {
    assert(m_kf_fitter != nullptr);

    // Setup the surface sequence buffer
    const unsigned int n_surfaces_per_track{
        std::max(cfg.max_track_candidates_per_track *
                     cfg.kalman_smoother.surface_sequence_size_factor,
                 cfg.kalman_smoother.min_surface_sequence_capacity)};
    std::vector<unsigned int> seqs_sizes(n_seeds, n_surfaces_per_track);

    // Not needed for PKF
    vecmem::data::vector_view<unsigned int> param_ids_view{};
    vecmem::data::vector_view<unsigned int> param_liveness_view{};

    kalman_fitting_algorithm::fit_payload tmp =
        m_kf_fitter->prepare_fit_payload(
            det, bfield, seqs_sizes,
            {param_ids_view, param_liveness_view, track_candidates_buffer});

    // Save the (non-templated) host payload.
    smoothing_payload.payload = tmp.payload;

    // Save all the type erased payloads into it.
    smoothing_payload.surfaces = std::move(tmp.surfaces);
    smoothing_payload.host_tpayload = tmp.host_tpayload;
    smoothing_payload.device_tpayload = std::move(tmp.device_tpayload);
  }

  /*****************************************************************
   * Progressive Kalman Filter
   *****************************************************************/

  if (cfg.run_pkf && (cfg.max_num_branches_per_seed == 1 ||
                      cfg.max_num_branches_per_surface == 1)) {
    // Create track candidate buffer
    vecmem::vector<unsigned int> n_constituent_links(mr().host);
    n_constituent_links.resize(n_seeds, cfg.max_track_candidates_per_track);

    track_candidates_buffer = typename edm::track_container<algebra_t>::buffer(
        {n_constituent_links, mr().main, mr().host,
         vecmem::data::buffer_type::resizable},
        {cfg.run_smoother == smoother_type::e_none
             ? 0u
             : n_seeds * cfg.max_track_candidates_per_track,
         mr().main, vecmem::data::buffer_type::resizable},
        measurements_view);
    copy().setup(track_candidates_buffer.tracks)->ignore();
    copy().setup(track_candidates_buffer.states)->ignore();

    smoothing_payload.payload.tracks = track_candidates_buffer;

    // Get output track statistics
    vecmem::data::vector_buffer<track_stats<scalar_t>> track_stats_buffer{
        n_seeds, mr().main};
    copy().setup(track_stats_buffer)->ignore();

    // Get the output track state data, depending on which data is
    // required
    vecmem::data::vector_buffer<track_state_candidate> track_cand_buffer{
        0, mr().main};
    vecmem::data::vector_buffer<filtered_track_state_candidate<algebra_t>>
        filtered_track_cand_buffer{0, mr().main};
    vecmem::data::vector_buffer<full_track_state_candidate<algebra_t>>
        full_track_cand_buffer{0, mr().main};

    // Allocate memory for the required data collection mode
    const unsigned int max_cands{seeds.size() *
                                 cfg.max_track_candidates_per_track};
    if (cfg.run_smoother == smoother_type::e_none) {
      track_cand_buffer = vecmem::data::vector_buffer<track_state_candidate>{
          max_cands, mr().main};
    } else if (cfg.run_smoother == smoother_type::e_kalman) {
      filtered_track_cand_buffer = vecmem::data::vector_buffer<
          filtered_track_state_candidate<algebra_t>>{max_cands, mr().main};
    } else if (cfg.run_smoother == smoother_type::e_mbf) {
      full_track_cand_buffer =
          vecmem::data::vector_buffer<full_track_state_candidate<algebra_t>>{
              max_cands, mr().main};
    }

    progressive_kalman_filter_kernel(
        n_seeds, cfg, det, bfield,
        {
            .seeds_view = in_params_buffer,
            .measurements_view = measurements_view,
            .measurement_ranges_view = meas_ranges_buffer,
            .track_stats_view = track_stats_buffer,
            .track_cand_view = track_cand_buffer,
            .filtered_track_cand_view = filtered_track_cand_buffer,
            .full_track_cand_view = full_track_cand_buffer,
            .tracks_view = {track_candidates_buffer},
        },
        smoothing_payload);
  }

  /*****************************************************************
   * Combinatorial Kalman Filter
   *****************************************************************/

  else {
    vecmem::data::vector_buffer<unsigned int> param_liveness_buffer =
        vecmem::data::vector_buffer<unsigned int>(n_seeds, mr().main);
    copy().setup(param_liveness_buffer)->ignore();
    copy().memset(param_liveness_buffer, 1)->ignore();

    // Number of tracks per seed
    vecmem::data::vector_buffer<unsigned int> n_tracks_per_seed_buffer(
        n_seeds, mr().main);
    copy().setup(n_tracks_per_seed_buffer)->ignore();

    // Compute the effective number of initial links per seed. If the
    // branching factor (`max_num_branches_per_surface`) is arbitrary there
    // is no useful upper bound on the number of links, but if the branching
    // factor is exactly one, we can never have more links per seed than the
    // number of CKF steps, which is a useful upper bound.
    const unsigned int effective_initial_links_per_seed =
        cfg.max_num_branches_per_surface == 1
            ? std::min(cfg.initial_links_per_seed,
                       cfg.max_track_candidates_per_track)
            : cfg.initial_links_per_seed;

    // Create a buffer of candidate links
    unsigned int link_buffer_capacity =
        effective_initial_links_per_seed * n_seeds;
    vecmem::data::vector_buffer<candidate_link> links_buffer(
        link_buffer_capacity, mr().main, vecmem::data::buffer_type::resizable);
    copy().setup(links_buffer)->ignore();

    // Buffers needed for MBF smoother (if enabled).
    vecmem::data::vector_buffer<bound_matrix<algebra_t>> jacobian_buffer,
        tmp_jacobian_buffer;
    bound_track_parameters_collection_types::buffer
        link_predicted_parameter_buffer,
        link_filtered_parameter_buffer;

    /*
     * If we are aiming to run the MBF smoother at the end of the track
     * finding, we need some space to store the intermediate Jacobians
     * and parameters. Allocate that space here.
     */
    if (cfg.run_smoother == smoother_type::e_mbf) {
      jacobian_buffer = {link_buffer_capacity, mr().main};
      link_predicted_parameter_buffer = {link_buffer_capacity, mr().main};
      link_filtered_parameter_buffer = {link_buffer_capacity, mr().main};
    }

    // Create a buffer of tip links
    const unsigned int tips_buffer_capacity =
        cfg.max_num_branches_per_seed * n_seeds;
    vecmem::data::vector_buffer<unsigned int> tips_buffer{
        tips_buffer_capacity, mr().main, vecmem::data::buffer_type::resizable};
    copy().setup(tips_buffer)->ignore();
    vecmem::data::vector_buffer<unsigned int> tip_length_buffer{
        tips_buffer_capacity, mr().main};
    copy().setup(tip_length_buffer)->ignore();

    std::map<unsigned int, unsigned int> step_to_link_idx_map;
    step_to_link_idx_map[0u] = 0u;

    unsigned int n_in_params = n_seeds;
    for (unsigned int step = 0;
         step < cfg.max_track_candidates_per_track && n_in_params > 0; ++step) {
      /*****************************************************************
       * Find valid tracks
       *****************************************************************/

      unsigned int n_candidates = 0;

      // Buffer for kalman-updated parameters spawned by the
      // measurement candidates
      const unsigned int n_max_candidates =
          n_in_params * cfg.max_num_branches_per_surface;

      bound_track_parameters_collection_types::buffer updated_params_buffer(
          n_max_candidates, mr().main);
      copy().setup(updated_params_buffer)->ignore();
      vecmem::data::vector_buffer<unsigned int> updated_liveness_buffer(
          n_max_candidates, mr().main);
      copy().setup(updated_liveness_buffer)->ignore();

      // Reset the number of tracks per seed
      copy().memset(n_tracks_per_seed_buffer, 0)->ignore();

      const unsigned int links_size = step_to_link_idx_map[step];

      // Ensure that the links buffer is large enough to hold all new
      // links.
      if (links_size + n_max_candidates > link_buffer_capacity) {
        const unsigned int new_link_buffer_capacity =
            std::max(2 * link_buffer_capacity, links_size + n_max_candidates);
        TRACCC_INFO("Link buffer (capacity "
                    << link_buffer_capacity << ") is too small to hold "
                    << links_size << " current and " << n_max_candidates
                    << " new links; increasing capacity to "
                    << new_link_buffer_capacity);

        link_buffer_capacity = new_link_buffer_capacity;

        vecmem::data::vector_buffer<candidate_link> new_links_buffer(
            link_buffer_capacity, mr().main,
            vecmem::data::buffer_type::resizable);
        copy().setup(new_links_buffer)->ignore();
        copy()(links_buffer, new_links_buffer)->wait();

        links_buffer = std::move(new_links_buffer);

        if (cfg.run_smoother == smoother_type::e_mbf) {
          // Create new, larger buffers for the MBF smoother data.
          vecmem::data::vector_buffer<bound_matrix<algebra_t>>
              new_jacobian_buffer{link_buffer_capacity, mr().main};
          bound_track_parameters_collection_types::buffer
              new_link_predicted_parameter_buffer{link_buffer_capacity,
                                                  mr().main},
              new_link_filtered_parameter_buffer{link_buffer_capacity,
                                                 mr().main};

          // Copy old data to new buffers.
          copy()(jacobian_buffer, new_jacobian_buffer)->ignore();
          copy()(link_predicted_parameter_buffer,
                 new_link_predicted_parameter_buffer)
              ->ignore();
          copy()(link_filtered_parameter_buffer,
                 new_link_filtered_parameter_buffer)
              ->wait();

          // Replace old buffers with the new ones.
          jacobian_buffer = std::move(new_jacobian_buffer);
          link_predicted_parameter_buffer =
              std::move(new_link_predicted_parameter_buffer);
          link_filtered_parameter_buffer =
              std::move(new_link_filtered_parameter_buffer);
        }
      }

      {
        vecmem::data::vector_buffer<unsigned int>
            out_params_per_in_param_buffer(n_in_params, mr().main);
        copy().setup(out_params_per_in_param_buffer)->ignore();
        vecmem::data::vector_buffer<unsigned int> out_params_index_buffer(
            n_in_params, mr().main);
        copy().setup(out_params_index_buffer)->ignore();
        vecmem::data::vector_buffer<candidate_link> tmp_links_buffer(
            n_max_candidates, mr().main);
        copy().setup(tmp_links_buffer)->ignore();
        bound_track_parameters_collection_types::buffer tmp_params_buffer(
            n_max_candidates, mr().main);
        copy().setup(tmp_params_buffer)->ignore();

        // Launch the track finding kernel.
        find_tracks_kernel(
            n_in_params, cfg, det,
            {.measurements_view = measurements_view,
             .in_params_view = in_params_buffer,
             .in_params_liveness_view = param_liveness_buffer,
             .n_in_params = n_in_params,
             .measurement_ranges_view = meas_ranges_buffer,
             .links_view = links_buffer,
             .prev_links_idx = (step == 0 ? 0 : step_to_link_idx_map[step - 1]),
             .step = step,
             .out_params_per_in_param_view = out_params_per_in_param_buffer,
             .tips_view = tips_buffer,
             .tip_lengths_view = tip_length_buffer,
             .n_tracks_per_seed_view = n_tracks_per_seed_buffer,
             .tmp_params_view = tmp_params_buffer,
             .tmp_links_view = tmp_links_buffer});

        // Launch the track condensing kernel.
        condense_tracks_kernel(
            n_in_params, out_params_per_in_param_buffer,
            out_params_index_buffer,
            {.n_in_params = n_in_params,
             .step = step,
             .curr_links_idx = step_to_link_idx_map[step],
             .max_num_branches_per_surface = cfg.max_num_branches_per_surface,
             .min_track_candidates_per_track =
                 cfg.min_track_candidates_per_track,
             .max_track_candidates_per_track =
                 cfg.max_track_candidates_per_track,
             .links_view = links_buffer,
             .in_params_view = in_params_buffer,
             .in_tmp_params_view = tmp_params_buffer,
             .in_tmp_links_view = tmp_links_buffer,
             .in_params_index_view = out_params_index_buffer,
             .out_params_view = updated_params_buffer,
             .out_params_liveness_view = updated_liveness_buffer,
             .tips_view = tips_buffer,
             .tip_lengths_view = tip_length_buffer,
             .jacobian_view = jacobian_buffer,
             .tmp_jacobian_view = tmp_jacobian_buffer,
             .link_predicted_parameter_view = link_predicted_parameter_buffer,
             .link_filtered_parameter_view = link_filtered_parameter_buffer});

        std::swap(in_params_buffer, updated_params_buffer);
        std::swap(param_liveness_buffer, updated_liveness_buffer);

        if (mr().host) {
          vecmem::async_size size = copy().get_size(links_buffer, *(mr().host));
          // Here we could give control back to the caller, once our
          // code allows for it. (coroutines...)
          step_to_link_idx_map[step + 1] = size.get();
        } else {
          step_to_link_idx_map[step + 1] = copy().get_size(links_buffer);
        }
        n_candidates =
            step_to_link_idx_map[step + 1] - step_to_link_idx_map[step];
      }

      // Set up the buffer used in the next few steps
      vecmem::data::vector_buffer<unsigned int> param_ids_buffer =
          vecmem::data::vector_buffer<unsigned int>(n_candidates, mr().main);
      copy().setup(param_ids_buffer)->ignore();

      /*
       * On later steps, we can duplicate removal which will attempt to
       * find tracks that are propagated multiple times and deduplicate
       * them.
       */
      if ((n_candidates > 0) &&
          (step >= cfg.duplicate_removal_minimum_length)) {
        vecmem::data::vector_buffer<unsigned int> link_last_measurement_buffer(
            n_candidates, mr().main);
        copy().setup(link_last_measurement_buffer)->ignore();

        /*
         * First, we sort the tracks by the index of their final
         * measurement which is critical to ensure good performance.
         */
        fill_finding_duplicate_removal_sort_keys_kernel(
            n_candidates,
            {.links_view = links_buffer,
             .param_liveness_view = param_liveness_buffer,
             .link_last_measurement_view = link_last_measurement_buffer,
             .param_ids_view = param_ids_buffer,
             .n_links = n_candidates,
             .curr_links_idx = step_to_link_idx_map[step],
             .n_measurements = n_measurements});

        sort_param_ids_by_last_measurement(link_last_measurement_buffer,
                                           param_ids_buffer);

        /*
         * Then, we run the actual duplicate removal kernel.
         */
        remove_duplicates_kernel(
            n_candidates, cfg,
            {.links_view = links_buffer,
             .link_last_measurement_view = link_last_measurement_buffer,
             .param_ids_view = param_ids_buffer,
             .param_liveness_view = param_liveness_buffer,
             .n_links = n_candidates,
             .curr_links_idx = step_to_link_idx_map[step],
             .n_measurements = n_measurements,
             .step = step});
      }

      // If no more CKF step is expected, the tips and links are
      // populated, and any further time-consuming action is avoided
      if (step == (cfg.max_track_candidates_per_track - 1)) {
        break;
      }

      if (n_candidates > 0) {
        /*****************************************************************
         * Get key and value for parameter sorting
         *****************************************************************/
        {
          vecmem::data::vector_buffer<device::sort_key> keys_buffer(
              n_candidates, mr().main);
          copy().setup(keys_buffer)->ignore();

          fill_finding_propagation_sort_keys_kernel(
              n_candidates, {.params_view = in_params_buffer,
                             .param_liveness_view = param_liveness_buffer,
                             .keys_view = keys_buffer,
                             .ids_view = param_ids_buffer});

          // Sort the key and values
          sort_param_ids_by_keys(keys_buffer, param_ids_buffer);
        }

        /*****************************************************************
         * Propagate to the next surface
         *****************************************************************/
        if (cfg.run_smoother == smoother_type::e_mbf) {
          tmp_jacobian_buffer = {n_candidates, mr().main};
        }

        propagate_to_next_surface_kernel(
            n_candidates, cfg, det, bfield,
            {.params_view = in_params_buffer,
             .params_liveness_view = param_liveness_buffer,
             .param_ids_view = param_ids_buffer,
             .links_view = links_buffer,
             .prev_links_idx = step_to_link_idx_map[step],
             .step = step,
             .n_in_params = n_candidates,
             .tips_view = tips_buffer,
             .tip_lengths_view = tip_length_buffer,
             .tmp_jacobian_view = tmp_jacobian_buffer});
      }

      n_in_params = n_candidates;
    }

    tmp_jacobian_buffer = {};

    TRACCC_DEBUG(
        "Final link buffer usage was "
        << copy().get_size(links_buffer) << " out of " << link_buffer_capacity
        << " ("
        << ((100.f * static_cast<float>(copy().get_size(links_buffer))) /
            static_cast<float>(link_buffer_capacity))
        << "%)");

    /*****************************************************************
     * Kernel6: Build tracks
     *****************************************************************/

    // Get the number of tips
    unsigned int n_tips_total = 0u;
    if (mr().host) {
      vecmem::async_size size = copy().get_size(tips_buffer, *(mr().host));
      // Here we could give control back to the caller, once our code
      // allows for it. (coroutines...)
      n_tips_total = size.get();
    } else {
      n_tips_total = copy().get_size(tips_buffer);
    }

    vecmem::data::vector_buffer<unsigned int> tip_to_output_map;

    if (n_tips_total > 0 && cfg.max_num_tracks_per_measurement > 0) {
      // TODO: DOCS

      vecmem::data::vector_buffer<unsigned int>
          best_tips_per_measurement_index_buffer(
              cfg.max_num_tracks_per_measurement * n_measurements, mr().main);
      copy().setup(best_tips_per_measurement_index_buffer)->ignore();

      vecmem::data::vector_buffer<unsigned long long int>
          best_tips_per_measurement_insertion_mutex_buffer(n_measurements,
                                                           mr().main);
      copy().setup(best_tips_per_measurement_insertion_mutex_buffer)->ignore();
      copy()
          .memset(best_tips_per_measurement_insertion_mutex_buffer, 0)
          ->ignore();

      {
        vecmem::data::vector_buffer<scalar>
            best_tips_per_measurement_pval_buffer(
                cfg.max_num_tracks_per_measurement * n_measurements, mr().main);
        copy().setup(best_tips_per_measurement_pval_buffer)->ignore();

        gather_best_tips_per_measurement_kernel(
            n_tips_total, {tips_buffer, links_buffer, measurements_view,
                           best_tips_per_measurement_insertion_mutex_buffer,
                           best_tips_per_measurement_index_buffer,
                           best_tips_per_measurement_pval_buffer,
                           cfg.max_num_tracks_per_measurement});
      }

      vecmem::data::vector_buffer<unsigned int> votes_per_tip_buffer(
          n_tips_total, mr().main);
      copy().setup(votes_per_tip_buffer)->ignore();
      copy().memset(votes_per_tip_buffer, 0)->ignore();

      gather_measurement_votes_kernel(
          cfg.max_num_tracks_per_measurement * n_measurements,
          {.insertion_mutex = best_tips_per_measurement_insertion_mutex_buffer,
           .tip_index = best_tips_per_measurement_index_buffer,
           .votes_per_tip = votes_per_tip_buffer,
           .max_num_tracks_per_measurement =
               cfg.max_num_tracks_per_measurement});

      tip_to_output_map =
          vecmem::data::vector_buffer<unsigned int>(n_tips_total, mr().main);
      copy().setup(tip_to_output_map)->ignore();

      vecmem::data::vector_buffer<unsigned int> new_tip_length_buffer{
          n_tips_total, mr().main, vecmem::data::buffer_type::resizable};
      copy().setup(new_tip_length_buffer)->ignore();

      update_tip_length_buffer_kernel(
          n_tips_total, {.old_tip_length = tip_length_buffer,
                         .new_tip_length = new_tip_length_buffer,
                         .measurement_votes = votes_per_tip_buffer,
                         .tip_to_output_map = tip_to_output_map,
                         .min_measurement_voting_fraction =
                             cfg.min_measurement_voting_fraction});

      tip_length_buffer = std::move(new_tip_length_buffer);
    }

    vecmem::vector<unsigned int> tips_length_host(mr().host);
    vecmem::copy::event_type ev = copy()(tip_length_buffer, tips_length_host);
    // Here we could give control back to the caller, once our code allows
    // for it. (coroutines...)
    ev->wait();
    // The following is only necessary if filtering was not turned on. Since
    // with filtering on, the buffer is resizable, so the host vector would
    // already have the correct size.
    if (cfg.max_num_tracks_per_measurement == 0) {
      tips_length_host.resize(n_tips_total);
    }

    unsigned int n_states = 0u;
    if (cfg.run_smoother == smoother_type::e_mbf) {
      n_states =
          std::accumulate(tips_length_host.begin(), tips_length_host.end(), 0u);
    }

    // Create track buffer
    track_candidates_buffer =
        typename edm::track_container<default_algebra>::buffer(
            {tips_length_host, mr().main, mr().host},
            {n_states, mr().main, vecmem::data::buffer_type::resizable},
            measurements_view);
    copy().setup(track_candidates_buffer.tracks)->ignore();
    copy().setup(track_candidates_buffer.states)->ignore();

    smoothing_payload.payload.tracks = track_candidates_buffer;

    // @Note: nBlocks can be zero in case there is no tip. This happens when
    // chi2_max cfg is set tightly and no tips are found
    if (n_tips_total > 0) {
      build_tracks_kernel(
          n_tips_total, cfg.run_smoother == smoother_type::e_mbf,
          cfg.meas_calibration,
          {.seeds_view = seeds,
           .links_view = links_buffer,
           .tips_view = tips_buffer,
           .tracks_view = {track_candidates_buffer},
           .tip_to_output_map = tip_to_output_map,
           .jacobian_ptr = jacobian_buffer.ptr(),
           .link_predicted_parameter_view = link_predicted_parameter_buffer,
           .link_filtered_parameter_view = link_filtered_parameter_buffer});
    }
  }

  /*****************************************************************
   * Run Fitting
   *****************************************************************/

  /// Run a Kalman filter in back propagation mode to smoothe the tracks
  if (cfg.run_smoother == smoother_type::e_kalman) {
    assert(m_kf_fitter != nullptr);

    m_kf_fitter->fit_backward_kernel(cfg.kalman_smoother, smoothing_payload);

  } else if (cfg.run_pkf && cfg.run_smoother == smoother_type::e_mbf) {
    std::cout << "FATAL: MBF smoother not implemented for "
                 "progressive filter!"
              << std::endl;
  }

  return track_candidates_buffer;
}

}  // namespace traccc::device
