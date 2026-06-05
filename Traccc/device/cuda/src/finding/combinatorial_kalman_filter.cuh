/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "../sanity/contiguous_on.cuh"
#include "../utils/barrier.hpp"
#include "../utils/cuda_error_handling.hpp"
#include "../utils/thread_id.hpp"
#include "../utils/utils.hpp"
#include "./kernels/build_tracks.cuh"
#include "./kernels/fill_finding_duplicate_removal_sort_keys.cuh"
#include "./kernels/fill_finding_propagation_sort_keys.cuh"
#include "./kernels/find_tracks.cuh"
#include "./kernels/gather_best_tips_per_measurement.cuh"
#include "./kernels/gather_measurement_votes.cuh"
#include "./kernels/propagate_to_next_surface.hpp"
#include "./kernels/remove_duplicates.cuh"
#include "./kernels/update_tip_length_buffer.cuh"

// Project include(s).
#include "traccc/edm/device/identity_projector.hpp"
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/track_container.hpp"
#include "traccc/edm/track_parameters.hpp"
#include "traccc/finding/candidate_link.hpp"
#include "traccc/finding/details/combinatorial_kalman_filter_types.hpp"
#include "traccc/finding/device/geo_id_surface_comparator.hpp"
#include "traccc/finding/finding_config.hpp"
#include "traccc/utils/logging.hpp"
#include "traccc/utils/memory_resource.hpp"
#include "traccc/utils/projections.hpp"
#include "traccc/utils/propagation.hpp"

// VecMem include(s).
#include <vecmem/memory/unique_ptr.hpp>
#include <vecmem/utils/copy.hpp>

// Thrust include(s).
#include <thrust/copy.h>
#include <thrust/execution_policy.h>
#include <thrust/fill.h>
#include <thrust/scan.h>
#include <thrust/sort.h>
#include <thrust/unique.h>

namespace traccc::cuda::details {

/// Templated implementation of the track finding algorithm.
///
/// Concrete track finding algorithms can use this function with the appropriate
/// specializations, to find tracks on top of a specific detector type, magnetic
/// field type, and track finding configuration.
///
/// @tparam detector_t The (device) detector type to use
/// @tparam bfield_t   The magnetic field type to use
///
/// @param det               A view of the detector object
/// @param field             The magnetic field object
/// @param measurements_view All measurements in an event
/// @param seeds_view        All seeds in an event to start the track finding
///                          with
/// @param config            The track finding configuration
/// @param mr                The memory resource(s) to use
/// @param copy              The copy object to use
/// @param log               The logger to use for message logging
/// @param stream            The CUDA stream to use for the operations
/// @param warp_size         The warp size of the used CUDA device
///
/// @return A buffer of the found track candidates
///
template <typename detector_t, typename bfield_t>
edm::track_container<typename detector_t::algebra_type>::buffer
combinatorial_kalman_filter(
    const typename detector_t::const_view_type& det, const bfield_t& field,
    const edm::measurement_collection::const_view& measurements_view,
    const bound_track_parameters_collection_types::const_view& seeds,
    const finding_config& config, const memory_resource& mr, vecmem::copy& copy,
    const Logger& log, stream& str, unsigned int warp_size) {

    const edm::measurement_collection::const_device measurements{
        measurements_view};

    assert(config.min_step_length_for_next_surface >
               math::fabs(config.propagation.navigation.intersection
                              .overstep_tolerance) &&
           "Min step length for the next surface should be higher than the "
           "overstep tolerance");
    assert(is_contiguous_on<
           vecmem::device_vector<const detray::geometry::identifier>>(
        device::identity_projector{}, mr.main, copy, str,
        measurements_view.template get<6>()));

    // Create a logger.
    auto logger = [&log]() -> const Logger& { return log; };

    /// Access the underlying CUDA stream.
    cudaStream_t stream = get_stream(str);

    vecmem::unique_alloc_ptr<unsigned int> size_staging_ptr =
        vecmem::make_unique_alloc<unsigned int>(*(mr.host));

    /// Thrust policy to use.
    auto thrust_policy =
        thrust::cuda::par_nosync(std::pmr::polymorphic_allocator(&(mr.main)))
            .on(stream);

    /*****************************************************************
     * Measurement Operations
     *****************************************************************/

    unsigned int n_measurements;

    if (mr.host) {
        const vecmem::async_size size =
            copy.get_size(measurements_view, *(mr.host));
        n_measurements = size.get();
    } else {
        n_measurements = copy.get_size(measurements_view);
    }

    // Access the detector view as a detector object
    detector_t device_det(det);
    const unsigned int n_surfaces{device_det.surfaces().size()};

    // Get upper bounds of measurement ranges per surface
    vecmem::data::vector_buffer<unsigned int> meas_ranges_buffer{n_surfaces,
                                                                 mr.main};
    copy.setup(meas_ranges_buffer)->ignore();
    vecmem::device_vector<unsigned int> measurement_ranges(meas_ranges_buffer);

    thrust::upper_bound(thrust_policy, measurements.surface_link().begin(),
                        // We have to use this ugly form here, because if the
                        // measurement collection is resizable (which it often
                        // is), the end() function cannot be used in host code.
                        measurements.surface_link().begin() + n_measurements,
                        device_det.surfaces().begin(),
                        device_det.surfaces().end(), measurement_ranges.begin(),
                        device::geo_id_surface_comparator{});

    unsigned int n_seeds;

    if (mr.host) {
        const vecmem::async_size size = copy.get_size(seeds, *(mr.host));
        n_seeds = size.get();
    } else {
        n_seeds = copy.get_size(seeds);
    }

    // Prepare input parameters with seeds
    bound_track_parameters_collection_types::buffer in_params_buffer(n_seeds,
                                                                     mr.main);
    copy.setup(in_params_buffer)->ignore();
    copy(seeds, in_params_buffer, vecmem::copy::type::device_to_device)
        ->ignore();
    vecmem::data::vector_buffer<unsigned int> param_liveness_buffer(n_seeds,
                                                                    mr.main);
    copy.setup(param_liveness_buffer)->ignore();
    copy.memset(param_liveness_buffer, 1)->ignore();

    // Number of tracks per seed
    vecmem::data::vector_buffer<unsigned int> n_tracks_per_seed_buffer(n_seeds,
                                                                       mr.main);
    copy.setup(n_tracks_per_seed_buffer)->ignore();

    // Compute the effective number of initial links per seed. If the
    // branching factor (`max_num_branches_per_surface`) is arbitrary there
    // is no useful upper bound on the number of links, but if the branching
    // factor is exactly one, we can never have more links per seed than the
    // number of CKF steps, which is a useful upper bound.
    const unsigned int effective_initial_links_per_seed =
        config.max_num_branches_per_surface == 1
            ? std::min(config.initial_links_per_seed,
                       config.max_track_candidates_per_track)
            : config.initial_links_per_seed;

    // Create a buffer for links
    unsigned int link_buffer_capacity =
        effective_initial_links_per_seed * n_seeds;
    vecmem::data::vector_buffer<candidate_link> links_buffer(
        link_buffer_capacity, mr.main, vecmem::data::buffer_type::resizable);
    copy.setup(links_buffer)->ignore();

    vecmem::unique_alloc_ptr<bound_matrix<typename detector_t::algebra_type>[]>
        jacobian_ptr = nullptr;
    bound_track_parameters_collection_types::buffer
        link_predicted_parameter_buffer(0, mr.main);
    bound_track_parameters_collection_types::buffer
        link_filtered_parameter_buffer(0, mr.main);

    /*
     * If we are aiming to run the MBF smoother at the end of the track
     * finding, we need some space to store the intermediate Jacobians
     * and parameters. Allocate that space here.
     */
    if (config.run_mbf_smoother) {
        jacobian_ptr = vecmem::make_unique_alloc<
            bound_matrix<typename detector_t::algebra_type>[]>(
            mr.main, link_buffer_capacity);
        link_predicted_parameter_buffer =
            bound_track_parameters_collection_types::buffer(
                link_buffer_capacity, mr.main);
        link_filtered_parameter_buffer =
            bound_track_parameters_collection_types::buffer(
                link_buffer_capacity, mr.main);
    }

    // Create a buffer of tip links
    vecmem::data::vector_buffer<unsigned int> tips_buffer{
        config.max_num_branches_per_seed * n_seeds, mr.main,
        vecmem::data::buffer_type::resizable};
    copy.setup(tips_buffer)->ignore();
    vecmem::data::vector_buffer<unsigned int> tip_length_buffer{
        config.max_num_branches_per_seed * n_seeds, mr.main};
    copy.setup(tip_length_buffer)->ignore();

    std::map<unsigned int, unsigned int> step_to_link_idx_map;
    step_to_link_idx_map[0] = 0;

    vecmem::unique_alloc_ptr<bound_matrix<typename detector_t::algebra_type>[]>
        tmp_jacobian_ptr = nullptr;

    unsigned int n_in_params = n_seeds;
    for (unsigned int step = 0;
         step < config.max_track_candidates_per_track && n_in_params > 0;
         step++) {

        /*****************************************************************
         * Kernel2: Find valid tracks
         *****************************************************************/

        unsigned int n_candidates = 0;

        // Buffer for kalman-updated parameters spawned by the
        // measurement candidates
        const unsigned int n_max_candidates =
            n_in_params * config.max_num_branches_per_surface;

        bound_track_parameters_collection_types::buffer updated_params_buffer(
            n_max_candidates, mr.main);
        copy.setup(updated_params_buffer)->ignore();

        vecmem::data::vector_buffer<unsigned int> updated_liveness_buffer(
            n_max_candidates, mr.main);
        copy.setup(updated_liveness_buffer)->ignore();

        // Reset the number of tracks per seed
        copy.memset(n_tracks_per_seed_buffer, 0)->ignore();

        const unsigned int links_size = step_to_link_idx_map[step];

        if (links_size + n_max_candidates > link_buffer_capacity) {
            const unsigned int new_link_buffer_capacity = std::max(
                2 * link_buffer_capacity, links_size + n_max_candidates);

            TRACCC_INFO("Link buffer (capacity "
                        << link_buffer_capacity << ") is too small to hold "
                        << links_size << " current and " << n_max_candidates
                        << " new links; increasing capacity to "
                        << new_link_buffer_capacity);

            link_buffer_capacity = new_link_buffer_capacity;

            vecmem::data::vector_buffer<candidate_link> new_links_buffer(
                link_buffer_capacity, mr.main,
                vecmem::data::buffer_type::resizable);

            copy.setup(new_links_buffer)->ignore();
            copy(links_buffer, new_links_buffer)->wait();

            links_buffer = std::move(new_links_buffer);

            if (config.run_mbf_smoother) {
                vecmem::unique_alloc_ptr<
                    bound_matrix<typename detector_t::algebra_type>[]>
                    new_jacobian_ptr = vecmem::make_unique_alloc<
                        bound_matrix<typename detector_t::algebra_type>[]>(
                        mr.main, link_buffer_capacity);
                bound_track_parameters_collection_types::buffer
                    new_link_predicted_parameter_buffer{link_buffer_capacity,
                                                        mr.main};
                bound_track_parameters_collection_types::buffer
                    new_link_filtered_parameter_buffer{link_buffer_capacity,
                                                       mr.main};

                TRACCC_CUDA_ERROR_CHECK(cudaMemcpyAsync(
                    new_jacobian_ptr.get(), jacobian_ptr.get(),
                    links_size *
                        sizeof(bound_matrix<typename detector_t::algebra_type>),
                    cudaMemcpyDeviceToDevice, stream));

                copy(link_predicted_parameter_buffer,
                     new_link_predicted_parameter_buffer)
                    ->wait();
                copy(link_filtered_parameter_buffer,
                     new_link_filtered_parameter_buffer)
                    ->wait();

                TRACCC_CUDA_ERROR_CHECK(cudaStreamSynchronize(stream));

                jacobian_ptr = std::move(new_jacobian_ptr);
                link_predicted_parameter_buffer =
                    std::move(new_link_predicted_parameter_buffer);
                link_filtered_parameter_buffer =
                    std::move(new_link_filtered_parameter_buffer);
            }
        }

        {
            vecmem::data::vector_buffer<candidate_link> tmp_links_buffer(
                n_max_candidates, mr.main);
            copy.setup(tmp_links_buffer)->ignore();
            bound_track_parameters_collection_types::buffer tmp_params_buffer(
                n_max_candidates, mr.main);
            copy.setup(tmp_params_buffer)->ignore();

            // Allocate the kernel's payload in host memory.
            using payload_t = device::find_tracks_payload<detector_t>;
            const payload_t host_payload{
                .det_data = det,
                .measurements_view = measurements_view,
                .in_params_view = in_params_buffer,
                .in_params_liveness_view = param_liveness_buffer,
                .n_in_params = n_in_params,
                .measurement_ranges_view = meas_ranges_buffer,
                .links_view = links_buffer,
                .prev_links_idx =
                    (step == 0 ? 0 : step_to_link_idx_map[step - 1]),
                .curr_links_idx = step_to_link_idx_map[step],
                .step = step,
                .out_params_view = updated_params_buffer,
                .out_params_liveness_view = updated_liveness_buffer,
                .tips_view = tips_buffer,
                .tip_lengths_view = tip_length_buffer,
                .n_tracks_per_seed_view = n_tracks_per_seed_buffer,
                .tmp_params_view = tmp_params_buffer,
                .tmp_links_view = tmp_links_buffer,
                .jacobian_ptr = jacobian_ptr.get(),
                .tmp_jacobian_ptr = tmp_jacobian_ptr.get(),
                .link_predicted_parameter_view =
                    link_predicted_parameter_buffer,
                .link_filtered_parameter_view = link_filtered_parameter_buffer};

            // The number of threads, blocks and shared memory to use.
            const unsigned int nThreads = warp_size * 2;
            const unsigned int nBlocks =
                (n_in_params + nThreads - 1) / nThreads;
            const std::size_t shared_size =
                nThreads * sizeof(unsigned long long int) +
                2 * nThreads * sizeof(std::pair<unsigned int, unsigned int>);

            // Run the kernel.
            find_tracks<detector_t>(nBlocks, nThreads, shared_size, stream,
                                    config, host_payload);
            TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());

            std::swap(in_params_buffer, updated_params_buffer);
            std::swap(param_liveness_buffer, updated_liveness_buffer);

            TRACCC_CUDA_ERROR_CHECK(cudaMemcpyAsync(
                size_staging_ptr.get(), links_buffer.size_ptr(),
                sizeof(unsigned int), cudaMemcpyDeviceToHost, stream));

            str.synchronize();

            step_to_link_idx_map[step + 1] = *size_staging_ptr;
            n_candidates =
                step_to_link_idx_map[step + 1] - step_to_link_idx_map[step];
        }

        /*
         * On later steps, we can duplicate removal which will attempt to find
         * tracks that are propagated multiple times and deduplicate them.
         */
        if (n_candidates > 0 &&
            step >= config.duplicate_removal_minimum_length) {
            vecmem::data::vector_buffer<unsigned int>
                link_last_measurement_buffer(n_candidates, mr.main);
            vecmem::data::vector_buffer<unsigned int> param_ids_buffer(
                n_candidates, mr.main);

            /*
             * First, we sort the tracks by the index of their final
             * measurement which is critical to ensure good performance.
             */
            {
                const unsigned int nThreads = 256;
                const unsigned int nBlocks =
                    (n_candidates + nThreads - 1) / nThreads;

                kernels::fill_finding_duplicate_removal_sort_keys<<<
                    nBlocks, nThreads, 0, stream>>>(
                    {.links_view = links_buffer,
                     .param_liveness_view = param_liveness_buffer,
                     .link_last_measurement_view = link_last_measurement_buffer,
                     .param_ids_view = param_ids_buffer,
                     .n_links = n_candidates,
                     .curr_links_idx = step_to_link_idx_map[step],
                     .n_measurements = n_measurements});

                TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
            }

            vecmem::device_vector<unsigned int> keys_device(
                link_last_measurement_buffer);
            vecmem::device_vector<unsigned int> param_ids_device(
                param_ids_buffer);
            thrust::sort_by_key(thrust_policy, keys_device.begin(),
                                keys_device.end(), param_ids_device.begin());

            /*
             * Then, we run the actual duplicate removal kernel.
             */
            {
                const unsigned int nThreads = 256;
                const unsigned int nBlocks =
                    (n_candidates + nThreads - 1) / nThreads;

                kernels::remove_duplicates<<<nBlocks, nThreads, 0, stream>>>(
                    config,
                    {.links_view = links_buffer,
                     .link_last_measurement_view = link_last_measurement_buffer,
                     .param_ids_view = param_ids_buffer,
                     .param_liveness_view = param_liveness_buffer,
                     .n_links = n_candidates,
                     .curr_links_idx = step_to_link_idx_map[step],
                     .n_measurements = n_measurements,
                     .step = step});
            }

            TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
        }

        // If no more CKF step is expected, the tips and links are populated,
        // and any further time-consuming action is avoided
        if (step == config.max_track_candidates_per_track - 1) {
            break;
        }

        if (n_candidates > 0) {
            /*****************************************************************
             * Kernel3: Get key and value for parameter sorting
             *****************************************************************/

            vecmem::data::vector_buffer<unsigned int> param_ids_buffer(
                n_candidates, mr.main);
            copy.setup(param_ids_buffer)->ignore();

            {
                vecmem::data::vector_buffer<device::sort_key> keys_buffer(
                    n_candidates, mr.main);
                copy.setup(keys_buffer)->wait();

                const unsigned int nThreads = warp_size * 2;
                const unsigned int nBlocks =
                    (n_candidates + nThreads - 1) / nThreads;
                kernels::fill_finding_propagation_sort_keys<<<nBlocks, nThreads,
                                                              0, stream>>>(
                    {.params_view = in_params_buffer,
                     .param_liveness_view = param_liveness_buffer,
                     .keys_view = keys_buffer,
                     .ids_view = param_ids_buffer});
                TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());

                // Sort the key and values
                vecmem::device_vector<device::sort_key> keys_device(
                    keys_buffer);
                vecmem::device_vector<unsigned int> param_ids_device(
                    param_ids_buffer);
                thrust::sort_by_key(thrust_policy, keys_device.begin(),
                                    keys_device.end(),
                                    param_ids_device.begin());
                str.synchronize();
            }

            /*****************************************************************
             * Kernel4: Propagate to the next surface
             *****************************************************************/

            {
                if (config.run_mbf_smoother) {
                    tmp_jacobian_ptr = vecmem::make_unique_alloc<
                        bound_matrix<typename detector_t::algebra_type>[]>(
                        mr.main, n_candidates);
                }

                // Allocate the kernel's payload in host memory.
                using payload_t = device::propagate_to_next_surface_payload<
                    traccc::details::ckf_propagator_t<detector_t, bfield_t>,
                    bfield_t>;
                const payload_t host_payload{
                    .det_data = det,
                    .field_data = field,
                    .params_view = in_params_buffer,
                    .params_liveness_view = param_liveness_buffer,
                    .param_ids_view = param_ids_buffer,
                    .links_view = links_buffer,
                    .prev_links_idx = step_to_link_idx_map[step],
                    .step = step,
                    .n_in_params = n_candidates,
                    .tips_view = tips_buffer,
                    .tip_lengths_view = tip_length_buffer,
                    .tmp_jacobian_ptr = tmp_jacobian_ptr.get()};

                const unsigned int nThreads = warp_size * 4;
                const unsigned int nBlocks =
                    (n_candidates + nThreads - 1) / nThreads;
                propagate_to_next_surface<
                    traccc::details::ckf_propagator_t<detector_t, bfield_t>,
                    bfield_t>(nBlocks, nThreads, 0, stream, config,
                              host_payload);
                TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());

                str.synchronize();
            }
        }

        n_in_params = n_candidates;
    }

    tmp_jacobian_ptr.reset();

    TRACCC_DEBUG("Final link buffer usage was "
                 << copy.get_size(links_buffer) << " out of "
                 << link_buffer_capacity << " ("
                 << ((100.f * static_cast<float>(copy.get_size(links_buffer))) /
                     static_cast<float>(link_buffer_capacity))
                 << "%)");

    /*****************************************************************
     * Kernel5: Build tracks
     *****************************************************************/

    // Get the number of tips
    unsigned int n_tips_total;

    if (mr.host) {
        const vecmem::async_size size = copy.get_size(tips_buffer, *(mr.host));
        n_tips_total = size.get();
    } else {
        n_tips_total = copy.get_size(tips_buffer);
    }

    vecmem::vector<unsigned int> tips_length_host(mr.host);
    vecmem::data::vector_buffer<unsigned int> tip_to_output_map;

    unsigned int n_tips_total_filtered = n_tips_total;

    if (n_tips_total > 0 && config.max_num_tracks_per_measurement > 0) {
        // TODO: DOCS

        vecmem::data::vector_buffer<unsigned int>
            best_tips_per_measurement_index_buffer(
                config.max_num_tracks_per_measurement * n_measurements,
                mr.main);
        copy.setup(best_tips_per_measurement_index_buffer)->wait();

        vecmem::data::vector_buffer<unsigned long long int>
            best_tips_per_measurement_insertion_mutex_buffer(n_measurements,
                                                             mr.main);
        copy.setup(best_tips_per_measurement_insertion_mutex_buffer)->wait();

        // NOTE: This memset assumes that an all-zero bit vector interpreted
        // as a floating point value has value zero, which is true for IEEE
        // 754 but might not be true for arbitrary float formats.
        copy.memset(best_tips_per_measurement_insertion_mutex_buffer, 0)
            ->wait();

        {
            vecmem::data::vector_buffer<scalar>
                best_tips_per_measurement_pval_buffer(
                    config.max_num_tracks_per_measurement * n_measurements,
                    mr.main);
            copy.setup(best_tips_per_measurement_pval_buffer)->wait();

            // NOTE: Normally, launching small blocks is a performance
            // antipattern, but there is little use to having larger blocks
            // here.
            const unsigned int num_threads = 32;
            const unsigned int num_blocks =
                (n_tips_total + num_threads - 1) / num_threads;

            kernels::gather_best_tips_per_measurement<<<num_blocks, num_threads,
                                                        0, stream>>>(
                {tips_buffer, links_buffer, measurements_view,
                 best_tips_per_measurement_insertion_mutex_buffer,
                 best_tips_per_measurement_index_buffer,
                 best_tips_per_measurement_pval_buffer,
                 config.max_num_tracks_per_measurement});

            TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
        }

        vecmem::data::vector_buffer<unsigned int> votes_per_tip_buffer(
            n_tips_total, mr.main);
        copy.setup(votes_per_tip_buffer)->wait();
        copy.memset(votes_per_tip_buffer, 0)->wait();

        {
            const unsigned int num_threads = 512;
            const unsigned int num_blocks =
                (config.max_num_tracks_per_measurement * n_measurements +
                 num_threads - 1) /
                num_threads;

            kernels::gather_measurement_votes<<<num_blocks, num_threads, 0,
                                                stream>>>(
                {best_tips_per_measurement_insertion_mutex_buffer,
                 best_tips_per_measurement_index_buffer, votes_per_tip_buffer,
                 config.max_num_tracks_per_measurement});

            TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
        }

        tip_to_output_map =
            vecmem::data::vector_buffer<unsigned int>(n_tips_total, mr.main);
        copy.setup(tip_to_output_map)->wait();

        {
            const unsigned int num_threads = 512;
            const unsigned int num_blocks =
                (n_tips_total + num_threads - 1) / num_threads;

            vecmem::data::vector_buffer<unsigned int> new_tip_length_buffer{
                n_tips_total, mr.main, vecmem::data::buffer_type::resizable};
            copy.setup(new_tip_length_buffer)->wait();

            kernels::update_tip_length_buffer<<<num_blocks, num_threads, 0,
                                                stream>>>(
                {tip_length_buffer, new_tip_length_buffer, votes_per_tip_buffer,
                 tip_to_output_map, config.min_measurement_voting_fraction});

            TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());

            if (mr.host) {
                vecmem::async_size size =
                    copy.get_size(new_tip_length_buffer, *(mr.host));
                // Here we could give control back to the caller, once our code
                // allows for it. (coroutines...)
                n_tips_total_filtered = size.get();
            } else {
                n_tips_total_filtered = copy.get_size(new_tip_length_buffer);
            }

            tip_length_buffer = std::move(new_tip_length_buffer);
        }
    }

    copy(tip_length_buffer, tips_length_host)->wait();
    tips_length_host.resize(n_tips_total_filtered);

    unsigned int n_states;

    if (config.run_mbf_smoother) {
        n_states = std::accumulate(tips_length_host.begin(),
                                   tips_length_host.end(), 0u);
    } else {
        n_states = 0;
    }

    // Create track candidate buffer
    typename edm::track_container<typename detector_t::algebra_type>::buffer
        track_candidates_buffer{
            {tips_length_host, mr.main, mr.host},
            {n_states, mr.main, vecmem::data::buffer_type::resizable},
            measurements_view};
    copy.setup(track_candidates_buffer.tracks)->ignore();
    copy.setup(track_candidates_buffer.states)->ignore();

    // @Note: nBlocks can be zero in case there is no tip. This happens when
    // chi2_max config is set tightly and no tips are found
    if (n_tips_total > 0) {
        const unsigned int nThreads = warp_size * 2;
        const unsigned int nBlocks = (n_tips_total + nThreads - 1) / nThreads;

        const device::build_tracks_payload payload{
            .seeds_view = seeds,
            .links_view = links_buffer,
            .tips_view = tips_buffer,
            .tracks_view = {track_candidates_buffer},
            .tip_to_output_map = tip_to_output_map,
            .jacobian_ptr = jacobian_ptr.get(),
            .link_predicted_parameter_view = link_predicted_parameter_buffer,
            .link_filtered_parameter_view = link_filtered_parameter_buffer,
        };
        kernels::build_tracks<<<nBlocks, nThreads, 0, stream>>>(
            config.run_mbf_smoother, config.meas_calibration, payload);
        TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());

        str.synchronize();
    }

    return track_candidates_buffer;
}

}  // namespace traccc::cuda::details
