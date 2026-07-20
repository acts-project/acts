/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// HACK: Fix for intel/llvm#15544
// As of Intel LLVM 2025.0, enabling an AMD SYCL target inadvertently sets the
// `__CUDA_ARCH__` preprocessor definition which breaks all sorts of internal
// logic in Thrust. Thus, we very selectively undefine the `__CUDA_ARCH__`
// definition when we are are compiling SYCL code using the Intel LLVM
// compiler. This can be removed when intel/llvm#15443 makes it into a OneAPI
// release.
#include <limits>

#include <vecmem/memory/device_atomic_ref.hpp>
#if defined(__INTEL_LLVM_COMPILER) && defined(SYCL_LANGUAGE_VERSION)
#undef __CUDA_ARCH__
#endif

// Project include(s).
#include "traccc/device/array_insertion_mutex.hpp"
#include "traccc/edm/track_state_helpers.hpp"
#include "traccc/finding/measurement_selector.hpp"
#include "traccc/fitting/kalman_filter/gain_matrix_updater.hpp"
#include "traccc/fitting/kalman_filter/is_line_visitor.hpp"
#include "traccc/fitting/status_codes.hpp"
#include "traccc/utils/logging.hpp"

// Detray include(s)
#include <detray/geometry/tracking_surface.hpp>

// Thrust include(s).
#include <thrust/binary_search.h>
#include <thrust/execution_policy.h>

namespace traccc::device {

template <typename detector_t, concepts::thread_id1 thread_id_t,
          concepts::barrier barrier_t>
TRACCC_HOST_DEVICE inline void find_tracks(
    const thread_id_t& thread_id, const barrier_t& barrier,
    const finding_config& cfg, typename detector_t::const_view_type det_data,
    const find_tracks_payload& payload,
    const find_tracks_shared_payload& shared_payload) {
  using algebra_t = typename detector_t::algebra_type;

  const unsigned int in_param_id = thread_id.getGlobalThreadIdX();
  const bool last_step = payload.step == cfg.max_track_candidates_per_track - 1;

  /*
   * Initialize all of the device vectors from their vecmem views.
   */
  detector_t det(det_data);
  edm::measurement_collection::const_device measurements(
      payload.measurements_view);
  bound_track_parameters_collection_types::const_device in_params(
      payload.in_params_view);
  vecmem::device_vector<const unsigned int> in_params_liveness(
      payload.in_params_liveness_view);
  vecmem::device_vector<candidate_link> links(payload.links_view);
  vecmem::device_vector<candidate_link> tmp_links(payload.tmp_links_view);
  bound_track_parameters_collection_types::device tmp_params(
      payload.tmp_params_view);
  vecmem::device_vector<unsigned int> out_params_per_in_param(
      payload.out_params_per_in_param_view);
  vecmem::device_vector<const unsigned int> meas_ranges(
      payload.measurement_ranges_view);
  vecmem::device_vector<unsigned int> tips(payload.tips_view);
  vecmem::device_vector<unsigned int> tip_lengths(payload.tip_lengths_view);
  vecmem::device_vector<unsigned int> n_tracks_per_seed(
      payload.n_tracks_per_seed_view);

  /*
   * Initialize the block-shared data; in particular, set the total size of
   * the candidate buffer to zero, and then set the number of candidates for
   * each parameter to zero.
   */
  if (thread_id.getLocalThreadIdX() == 0) {
    shared_payload.shared_candidates_size = 0;
    shared_payload.shared_num_out_params = 0;
  }

  shared_payload.shared_insertion_mutex[thread_id.getLocalThreadIdX()] =
      encode_insertion_mutex(false, 0, 0.f);

  barrier.blockBarrier();

  /*
   * Step 1 of this kernel is to determine which indices belong to which
   * parameter. Because the measurements are guaranteed to be grouped, we can
   * simply find the first measurement's index and the total number of
   * indices.
   *
   * This entire step is executed on a one-thread-one-parameter model.
   */
  unsigned int init_meas = 0;
  unsigned int num_meas = 0;

  if (in_param_id < payload.n_in_params &&
      in_params_liveness.at(in_param_id) > 0u) {
    /*
     * Get the geo ID of this thread's parameters, then find the first
     * measurement that matches it.
     */
    const unsigned int sf_idx{in_params.at(in_param_id).surface_link().index()};
    init_meas = sf_idx == 0u ? 0u : meas_ranges[sf_idx - 1];
    num_meas = meas_ranges[sf_idx] - init_meas;
    TRACCC_VERBOSE_DEVICE("No. measurements: %d", num_meas);
    TRACCC_VERBOSE_DEVICE("Testing measurement: %d", init_meas);

    /*
     * If we cannot find any corresponding measurements, set the number of
     * measurements to zero.
     */
    init_meas = (num_meas == 0u) ? 0u : init_meas;
  }

  /*
   * Step 2 of this kernel involves processing the candidate measurements and
   * updating them on their corresponding surface.
   *
   * Because the number of measurements per parameter can vary wildly
   * (between 0 and 20), a naive one-thread-one-parameter model would incur a
   * lot of thread divergence here. Instead, we use a load-balanced model in
   * which threads process each others' measurements.
   *
   * The core idea is that each thread places its measurements into a shared
   * pool. We keep track of how many measurements each thread has placed into
   * the pool.
   */
  unsigned int curr_meas = 0;

  /*
   * This loop keeps running until all threads have processed all of their
   * measurements.
   */
  while (barrier.blockOr(curr_meas < num_meas ||
                         shared_payload.shared_candidates_size > 0)) {
    /*
     * The outer loop consists of three general components. The first
     * components is that each thread starts to fill a shared buffer of
     * measurements. The buffer is twice the size of the block to
     * accomodate any overflow.
     *
     * Threads insert their measurements into the shared buffer until they
     * either run out of measurements, or until the shared buffer is full.
     */
    for (; curr_meas < num_meas &&
           shared_payload.shared_candidates_size < thread_id.getBlockDimX();
         curr_meas++) {
      unsigned int idx =
          vecmem::device_atomic_ref<unsigned int,
                                    vecmem::device_address_space::local>(
              shared_payload.shared_candidates_size)
              .fetch_add(1u);

      /*
       * The buffer elemements are tuples of the measurement index and
       * the index of the thread that originally inserted that
       * measurement.
       */
      shared_payload.shared_candidates[idx] = {init_meas + curr_meas,
                                               thread_id.getLocalThreadIdX()};
    }

    barrier.blockBarrier();

    std::optional<std::tuple<
        typename edm::track_state_collection<algebra_t>::device::object_type,
        unsigned int, unsigned int>>
        result = std::nullopt;

    /*
     * The shared buffer is now full; each thread picks out zero or one of
     * the measurements and processes it.
     */
    if (thread_id.getLocalThreadIdX() < shared_payload.shared_candidates_size) {
      const unsigned int owner_local_thread_id =
          shared_payload.shared_candidates[thread_id.getLocalThreadIdX()]
              .second;
      const unsigned int owner_global_thread_id =
          owner_local_thread_id +
          thread_id.getBlockDimX() * thread_id.getBlockIdX();
      assert(in_params_liveness.at(owner_global_thread_id) != 0u);
      const bound_track_parameters<>& in_par =
          in_params.at(owner_global_thread_id);
      const unsigned int meas_idx =
          shared_payload.shared_candidates[thread_id.getLocalThreadIdX()].first;
      const unsigned int prev_link_idx =
          payload.prev_links_idx + owner_global_thread_id;
      const unsigned int seed_idx = payload.step > 0
                                        ? links.at(prev_link_idx).seed_idx
                                        : owner_global_thread_id;

      if (n_tracks_per_seed.at(seed_idx) < cfg.max_num_branches_per_seed) {
        const detray::tracking_surface sf{det, in_par.surface_link()};
        const bool is_line = detail::is_line(sf);

        const edm::measurement meas = measurements.at(meas_idx);

        const traccc::scalar chi2 = measurement_selector::predicted_chi2(
            meas, in_par, cfg.meas_calibration, is_line);

        if (chi2 <= cfg.chi2_max && chi2 >= 0.f) {
          edm::track_state trk_state =
              edm::make_track_state<algebra_t>(measurements, meas_idx);
          trk_state.filtered_chi2() = chi2;

          // Kalman filter status code
          kalman_fitter_status res{kalman_fitter_status::ERROR_OTHER};

          if (payload.step == 0 && !sf.has_material()) {
            // Only do this for the actual seed measurement
            res = kalman_fitter_status::SUCCESS;

            trk_state.filtered_params() = in_par;

            // Update measurement covariance
            const auto V =
                measurement_selector::calibrated_measurement_covariance<
                    algebra_t, 2>(meas, cfg.meas_calibration);

            auto& filtered_cov = trk_state.filtered_params().covariance();
            getter::element(filtered_cov, e_bound_loc0, e_bound_loc0) =
                getter::element(V, 0, 0);
            getter::element(filtered_cov, e_bound_loc1, e_bound_loc1) =
                getter::element(V, 1, 1);
          } else {
            // Run the Kalman update on a copy of the track
            // parameters
            res = gain_matrix_updater<algebra_t>{}(
                trk_state, meas, in_par, cfg.meas_calibration, is_line);
          }

          TRACCC_DEBUG_DEVICE("KF status: %d", res);

          /*
           * The $\chi^2$ value from the Kalman update should be less
           * than `chi2_max`, and the fit should have succeeded. If
           * both conditions are true, we emplace the state, the
           * measurement index, and the thread ID into an optional
           * value.
           *
           * NOTE: Using the optional value here allows us to remove
           * the depth of if-statements which is important for code
           * quality but, more importantly, allows us to more easily
           * use block-wide synchronization primitives.
           */
          if (res == kalman_fitter_status::SUCCESS) {
            TRACCC_VERBOSE_DEVICE("Found measurement: %d", meas_idx);
            result.emplace(trk_state, meas_idx, owner_local_thread_id);
          }
        }
      }
    }

    /*
     * Now comes the stage in which we add the parameters to the temporary
     * array, in such a way that we keep the best ones. This loop has a
     * barrier to ensure both thread safety and forward progress.
     *
     * NOTE: This has to be a loop because the software is set up such
     * that only one thread can write to the array of one input
     * parameter per loop cycle. Thus, the loop is here to resolve any
     * contention.
     */
    while (barrier.blockOr(result.has_value())) {
      /*
       * Threads which have no parameter stored (either because they
       * never had one or because they already deposited) do not have to
       * do anything.
       */
      if (result.has_value()) {
        /*
         * First, we reconstruct some necessary information from the
         * data that we stored previously.
         */
        const unsigned int meas_idx = std::get<1>(*result);
        const unsigned int owner_local_thread_id = std::get<2>(*result);
        const unsigned int owner_global_thread_id =
            owner_local_thread_id +
            thread_id.getBlockDimX() * thread_id.getBlockIdX();
        const traccc::scalar chi2 = std::get<0>(*result).filtered_chi2();
        assert(chi2 >= 0.f);
        unsigned long long int* mutex_ptr =
            &shared_payload.shared_insertion_mutex[owner_local_thread_id];
        const unsigned int prev_link_idx =
            payload.prev_links_idx + owner_global_thread_id;

        /*
         * The current thread will attempt to get a lock on the
         * output array for the input parameter ID which it is now
         * holding. If it manages to do so, the `index` variable will
         * be set to a value smaller than or equal to the maximum
         * number of elements; otherwise, it will be set to
         * `UINT_MAX`.
         */
        unsigned int index = std::numeric_limits<unsigned int>::max();
        unsigned long long int desired = 0;

        /*
         * We fetch and decode whatever the mutex state is at the
         * current time. The mutex is a 64-bit integer containing the
         * following:
         *
         * [00:31] A 32-bit IEEE 754 floating point number that equals
         *         the highest $\chi^2$ value among parameters
         *         currently stored.
         * [32:62] A 31-bit unsigned integer representing the number
         *         of parameters currently stored.
         * [63:63] A boolean that, if true, indicates that a thread is
         *         currently operating on the array guarded.
         */
        unsigned long long int assumed = *mutex_ptr;
        auto [locked, size, max] = decode_insertion_mutex(assumed);

        /*
         * If the array is already full _and_ our parameter has a
         * higher $\chi^2$ value than any of the elements in the
         * array, we can discard the current track state.
         */
        if (size >= cfg.max_num_branches_per_surface && chi2 > max) {
          result.reset();
        }

        /*
         * If we still have a track after the previous check, we will
         * try to add this. We can only do this if the mutex is not
         * locked.
         */
        if (result.has_value() && !locked) {
          desired = encode_insertion_mutex(true, size, max);

          /*
           * Attempt to CAS the mutex with the same value as before
           * but with the lock bit switched. If this succeeds (e.g.
           * the return value is as we assumed) then we have succes
           * fully locked and we set the `index` variable, which
           * indicates that we have the lock.
           */
          if (vecmem::device_atomic_ref<unsigned long long,
                                        vecmem::device_address_space::local>(
                  *mutex_ptr)
                  .compare_exchange_strong(assumed, desired)) {
            index = size;
          }
        }

        /*
         * If `index` is not `UINT32_MAX`, we are in the green to
         * write to the parameter array!
         */
        if (index != std::numeric_limits<unsigned int>::max()) {
          assert(result.has_value());
          assert(index <= cfg.max_num_branches_per_surface);

          /*
           * We will now proceed to find the index in the temporary
           * array that we will write to. There are two distinct
           * cases:
           *
           * 1. If `index` is the maximum branching value, then the
           *    array is already full, and we need to replace the
           *    worst existing parameter.
           * 2. If `index` is less than the maximum branching value,
           *    we can trivially insert the value at index.
           */
          unsigned int l_pos = std::numeric_limits<unsigned int>::max();
          const unsigned int p_offset =
              owner_global_thread_id * cfg.max_num_branches_per_surface;
          // No need to initialize this next variable. It always gets
          // a valid value in the proceeding expressions.
          float new_max;

          if (index == cfg.max_num_branches_per_surface) {
            /*
             * Handle the case in which we need to replace a
             * value; find the worst existing parameter and then
             * replace it. Also keep track of what the new maximum
             * $\chi^2$ value will be.
             */
            traccc::scalar highest = static_cast<traccc::scalar>(-1.f);

            for (unsigned int i = 0; i < cfg.max_num_branches_per_surface;
                 ++i) {
              const traccc::scalar old_chi2 = tmp_links.at(p_offset + i).chi2;

              if (old_chi2 > highest) {
                highest = old_chi2;
                l_pos = i;
              }
            }

            assert(l_pos != std::numeric_limits<unsigned int>::max());

            new_max = static_cast<float>(chi2);

            for (unsigned int i = 0; i < cfg.max_num_branches_per_surface;
                 ++i) {
              const traccc::scalar old_chi2 = tmp_links.at(p_offset + i).chi2;

              if (i != l_pos && old_chi2 > new_max) {
                new_max = static_cast<float>(old_chi2);
              }

              assert(old_chi2 <= tmp_links.at(p_offset + l_pos).chi2);
            }

            assert(chi2 <= new_max);
          } else {
            l_pos = index;
            new_max = std::max(static_cast<float>(chi2), max);
          }

          assert(l_pos < cfg.max_num_branches_per_surface);

          /*
           * Now, simply insert the temporary link at the found
           * position. Different cases for step 0 and other steps.
           */
          const unsigned int n_skipped =
              payload.step == 0 ? 0 : links.at(prev_link_idx).n_skipped;
          const unsigned int seed_idx = payload.step > 0
                                            ? links.at(prev_link_idx).seed_idx
                                            : owner_global_thread_id;
          const scalar prev_chi2_sum =
              payload.step > 0 ? links.at(prev_link_idx).chi2_sum : 0.f;
          const unsigned int prev_ndf_sum =
              payload.step > 0 ? links.at(prev_link_idx).ndf_sum : 0;

          /*
           * We add the new link under one of two conditions. First,
           * if the best-of array is not yet full (i.e. if the index
           * variable is not the maximum number of branches), we add
           * it unconditionally. If the array _is_ full, we add the
           * link if we have a link with a lower $\chi^2$ value than
           * the previous highest value. If the value is exactly
           * equal, we use a tie-breaking mechanism, namely a
           * comparison of the measurement indices. This tie breaker
           * should be extremely rare and should not bias the
           * physics results, but helps ensure that the output of
           * this algorithm is deterministic.
           */
          if (index != cfg.max_num_branches_per_surface ||
              chi2 < tmp_links.at(p_offset + l_pos).chi2 ||
              (chi2 == tmp_links.at(p_offset + l_pos).chi2 &&
               meas_idx < tmp_links.at(p_offset + l_pos).meas_idx)) {
            tmp_links.at(p_offset + l_pos) = {
                .step = payload.step,
                .previous_candidate_idx = prev_link_idx,
                .meas_idx = meas_idx,
                .seed_idx = seed_idx,
                .n_skipped = n_skipped,
                .n_consecutive_skipped = 0,
                .chi2 = chi2,
                .chi2_sum = prev_chi2_sum + chi2,
                .ndf_sum =
                    prev_ndf_sum +
                    measurements.at(std::get<0>(*result).measurement_index())
                        .dimensions()};

            tmp_params.at(p_offset + l_pos) =
                std::get<0>(*result).filtered_params();
          }

          /*
           * Reset the temporary state storage, as this is no longer
           * necessary; this implies that this thread will not try
           * to insert anything in the next loop iteration.
           */
          result.reset();

          unsigned int new_size =
              size < cfg.max_num_branches_per_surface ? size + 1 : size;

          /*
           * Release the lock using another atomic CAS operation.
           * Because nobody should be writing to this value, it
           * should always succeed!
           */
          [[maybe_unused]] bool cas_result =
              vecmem::device_atomic_ref<unsigned long long,
                                        vecmem::device_address_space::local>(
                  *mutex_ptr)
                  .compare_exchange_strong(
                      desired,
                      encode_insertion_mutex(false, new_size, new_max));

          assert(cas_result);
        }
      }
    }

    barrier.blockBarrier();

    /*
     * The reason the buffer is twice the size of the block is that we
     * might end up having some spill-over; this spill-over should be moved
     * to the front of the buffer.
     */
    shared_payload.shared_candidates[thread_id.getLocalThreadIdX()] =
        shared_payload.shared_candidates[thread_id.getLocalThreadIdX() +
                                         thread_id.getBlockDimX()];

    if (thread_id.getLocalThreadIdX() == 0) {
      if (shared_payload.shared_candidates_size >= thread_id.getBlockDimX()) {
        shared_payload.shared_candidates_size -= thread_id.getBlockDimX();
      } else {
        shared_payload.shared_candidates_size = 0;
      }
    }
  }

  /*
   * NOTE: A synchronization point here is unnecessary, as it is implicit in
   * the condition of the while-loop above.
   */

  unsigned int prev_link_idx = std::numeric_limits<unsigned int>::max();
  unsigned int seed_idx = std::numeric_limits<unsigned int>::max();
  unsigned int n_skipped = std::numeric_limits<unsigned int>::max();
  unsigned int n_consecutive_skipped = std::numeric_limits<unsigned int>::max();
  unsigned int prev_ndf_sum = 0u;
  scalar prev_chi2_sum = 0.f;

  unsigned int local_num_params = 0;

  bool in_param_can_create_hole = false;

  const bool in_param_is_live = in_param_id < payload.n_in_params &&
                                in_params_liveness.at(in_param_id) > 0u;

  if (in_param_is_live) {
    prev_link_idx = payload.prev_links_idx + in_param_id;
    seed_idx =
        payload.step > 0 ? links.at(prev_link_idx).seed_idx : in_param_id;
    n_skipped = payload.step == 0 ? 0 : links.at(prev_link_idx).n_skipped;
    n_consecutive_skipped =
        payload.step == 0 ? 0 : links.at(prev_link_idx).n_consecutive_skipped;
    in_param_can_create_hole =
        (n_skipped < cfg.max_num_skipping_per_cand) &&
        (n_consecutive_skipped < cfg.max_num_consecutive_skipped) &&
        (!last_step);
    prev_ndf_sum = payload.step == 0 ? 0 : links.at(prev_link_idx).ndf_sum;
    prev_chi2_sum = payload.step == 0 ? 0.f : links.at(prev_link_idx).chi2_sum;
  }

  /*
   * Compute the offset at which this block will write, as well as the
   * index at which this block will write.
   */
  if (in_param_is_live) {
    assert(prev_link_idx != std::numeric_limits<unsigned int>::max());
    assert(n_skipped != std::numeric_limits<unsigned int>::max());

    local_num_params = std::get<1>(decode_insertion_mutex(
        shared_payload.shared_insertion_mutex[thread_id.getLocalThreadIdX()]));

    /*
     * If we found zero parameters and we can create a hole, add that hole
     * to the temporary link and parameter lists.
     */
    if (local_num_params == 0 && in_param_can_create_hole) {
      const unsigned int in_offset =
          thread_id.getGlobalThreadIdX() * cfg.max_num_branches_per_surface;

      tmp_links.at(in_offset) = {
          .step = payload.step,
          .previous_candidate_idx = prev_link_idx,
          .meas_idx = std::numeric_limits<unsigned int>::max(),
          .seed_idx = seed_idx,
          .n_skipped = n_skipped + 1,
          .n_consecutive_skipped = n_consecutive_skipped + 1,
          .chi2 = std::numeric_limits<traccc::scalar>::max(),
          .chi2_sum = prev_chi2_sum,
          .ndf_sum = prev_ndf_sum};

      tmp_params.at(in_offset) = in_params.at(in_param_id);

      /*
       * If we created a hole, we now have a single output parameter!
       */
      local_num_params = 1;
    }

    /*
     * NOTE: We always create at least one state, because we also create
     * hole states for nodes which don't find any good compatible
     * measurements.
     */
    if (local_num_params > 0) {
      vecmem::device_atomic_ref<unsigned int> num_tracks_per_seed(
          n_tracks_per_seed.at(seed_idx));
      const unsigned int params_to_add = std::min(
          local_num_params,
          cfg.max_num_branches_per_seed -
              std::min(cfg.max_num_branches_per_seed,
                       num_tracks_per_seed.fetch_add(local_num_params)));

      out_params_per_in_param.at(in_param_id) = params_to_add;

      vecmem::device_atomic_ref<unsigned int,
                                vecmem::device_address_space::local>(
          shared_payload.shared_num_out_params)
          .fetch_add(params_to_add);
    } else {
      out_params_per_in_param.at(in_param_id) = 0;

      assert(!in_param_can_create_hole);
      const unsigned int n_cands = payload.step - n_skipped;

      if (n_cands >= cfg.min_track_candidates_per_track) {
        TRACCC_WARNING_DEVICE("Create tip: Max no. holes");
        auto tip_pos = tips.push_back(prev_link_idx);
        tip_lengths.at(tip_pos) = n_cands;
      }
    }
  } else if (in_param_id < payload.n_in_params) {
    out_params_per_in_param.at(in_param_id) = 0;
  }

  barrier.blockBarrier();

  if (thread_id.getLocalThreadIdX() == 0) {
    links.bulk_append_implicit(shared_payload.shared_num_out_params);
  }
}

}  // namespace traccc::device
