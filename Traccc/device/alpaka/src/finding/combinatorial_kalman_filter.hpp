/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "../utils/barrier.hpp"
#include "../utils/parallel_algorithms.hpp"
#include "../utils/thread_id.hpp"
#include "../utils/utils.hpp"

// Project include(s).
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/track_container.hpp"
#include "traccc/finding/actors/ckf_aborter.hpp"
#include "traccc/finding/candidate_link.hpp"
#include "traccc/finding/details/combinatorial_kalman_filter_types.hpp"
#include "traccc/finding/device/build_tracks.hpp"
#include "traccc/finding/device/fill_finding_duplicate_removal_sort_keys.hpp"
#include "traccc/finding/device/fill_finding_propagation_sort_keys.hpp"
#include "traccc/finding/device/find_tracks.hpp"
#include "traccc/finding/device/geo_id_surface_comparator.hpp"
#include "traccc/finding/device/propagate_to_next_surface.hpp"
#include "traccc/finding/device/remove_duplicates.hpp"
#include "traccc/finding/finding_config.hpp"
#include "traccc/utils/memory_resource.hpp"
#include "traccc/utils/projections.hpp"
#include "traccc/utils/propagation.hpp"

// VecMem include(s).
#include <vecmem/utils/copy.hpp>

namespace traccc::alpaka::details {
namespace kernels {

/// Alpaka kernel functor for @c traccc::device::find_tracks
template <typename detector_t>
struct find_tracks {
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc, const finding_config cfg,
        const device::find_tracks_payload<detector_t>* payload) const {

        auto& shared_num_out_params =
            ::alpaka::declareSharedVar<unsigned int, __COUNTER__>(acc);
        auto& shared_out_offset =
            ::alpaka::declareSharedVar<unsigned int, __COUNTER__>(acc);
        auto& shared_candidates_size =
            ::alpaka::declareSharedVar<unsigned int, __COUNTER__>(acc);
        unsigned long long int* const s =
            ::alpaka::getDynSharedMem<unsigned long long int>(acc);
        unsigned long long int* shared_insertion_mutex = s;

        alpaka::barrier<TAcc> barrier(&acc);
        details::thread_id1 thread_id(acc);

        const unsigned int blockDimX = thread_id.getBlockDimX();
        std::pair<unsigned int, unsigned int>* shared_candidates =
            reinterpret_cast<std::pair<unsigned int, unsigned int>*>(
                &shared_insertion_mutex[blockDimX]);

        device::find_tracks<detector_t>(
            thread_id, barrier, cfg, *payload,
            {shared_num_out_params, shared_out_offset, shared_insertion_mutex,
             shared_candidates, shared_candidates_size});
    }
};

/// Alpaka kernel functor for
/// @c traccc::device::fill_finding_duplicate_removal_sort_keys
struct fill_finding_duplicate_removal_sort_keys {
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc,
        const device::fill_finding_duplicate_removal_sort_keys_payload& payload)
        const {

        const device::global_index_t globalThreadIdx =
            ::alpaka::getIdx<::alpaka::Grid, ::alpaka::Threads>(acc)[0];
        device::fill_finding_duplicate_removal_sort_keys(globalThreadIdx,
                                                         payload);
    }
};

/// Alpaka kernel functor for @c traccc::device::remove_duplicates
struct remove_duplicates {
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc, const finding_config& cfg,
        const device::remove_duplicates_payload& payload) const {

        const device::global_index_t globalThreadIdx =
            ::alpaka::getIdx<::alpaka::Grid, ::alpaka::Threads>(acc)[0];
        device::remove_duplicates(globalThreadIdx, cfg, payload);
    }
};

/// Alpaka kernel functor for
/// @c traccc::device::fill_finding_propagation_sort_keys
struct fill_finding_propagation_sort_keys {
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc,
        const device::fill_finding_propagation_sort_keys_payload payload)
        const {

        const device::global_index_t globalThreadIdx =
            ::alpaka::getIdx<::alpaka::Grid, ::alpaka::Threads>(acc)[0];
        device::fill_finding_propagation_sort_keys(globalThreadIdx, payload);
    }
};

/// Alpaka kernel functor for @c traccc::device::propagate_to_next_surface
template <typename propagator_t, typename bfield_t>
struct propagate_to_next_surface {
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc, const finding_config cfg,
        const device::propagate_to_next_surface_payload<propagator_t, bfield_t>*
            payload) const {

        device::global_index_t globalThreadIdx =
            ::alpaka::getIdx<::alpaka::Grid, ::alpaka::Threads>(acc)[0];
        device::propagate_to_next_surface<propagator_t, bfield_t>(
            globalThreadIdx, cfg, *payload);
    }
};

/// Alpaka kernel functor for @c traccc::device::build_tracks
struct build_tracks {
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc, bool run_mbf,
        const measurement_selector::config& calib_cfg,
        const device::build_tracks_payload payload) const {

        device::global_index_t globalThreadIdx =
            ::alpaka::getIdx<::alpaka::Grid, ::alpaka::Threads>(acc)[0];
        device::build_tracks(globalThreadIdx, run_mbf, calib_cfg, payload);
    }
};

}  // namespace kernels

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
/// @param queue             The Alpaka queue to use
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
    const Logger& log, Queue& queue) {

    const edm::measurement_collection::const_device measurements{
        measurements_view};

    assert(config.min_step_length_for_next_surface >
               math::fabs(config.propagation.navigation.intersection
                              .overstep_tolerance) &&
           "Min step length for the next surface should be higher than the "
           "overstep tolerance");

    // Create a logger.
    auto logger = [&log]() -> const Logger& { return log; };

    // Number of threads per block to use.
    const Idx threadsPerBlock = getWarpSize<Acc>() * 2;

    /*****************************************************************
     * Measurement Operations
     *****************************************************************/

    const unsigned int n_measurements = copy.get_size(measurements_view);

    // Access the detector view as a detector object
    detector_t device_det(det);
    const unsigned int n_surfaces{device_det.surfaces().size()};

    // Get upper bounds of measurement ranges per surface
    vecmem::data::vector_buffer<unsigned int> meas_ranges_buffer{n_surfaces,
                                                                 mr.main};
    copy.setup(meas_ranges_buffer)->ignore();
    vecmem::device_vector<unsigned int> measurement_ranges(meas_ranges_buffer);

    // Get upper bounds of measurement ranges
    details::upper_bound(
        queue, mr, measurements.surface_link().begin(),
        // We have to use this ugly form here, because if the
        // measurement collection is resizable (which it often
        // is), the end() function cannot be used in host code.
        measurements.surface_link().begin() + n_measurements,
        device_det.surfaces().begin(), device_det.surfaces().end(),
        measurement_ranges.begin(), device::geo_id_surface_comparator{});

    const unsigned int n_seeds = copy.get_size(seeds);

    // Prepare input parameters with seeds
    bound_track_parameters_collection_types::buffer in_params_buffer(n_seeds,
                                                                     mr.main);
    copy.setup(in_params_buffer)->wait();
    copy(seeds, in_params_buffer, vecmem::copy::type::device_to_device)->wait();
    vecmem::data::vector_buffer<unsigned int> param_liveness_buffer(n_seeds,
                                                                    mr.main);
    copy.setup(param_liveness_buffer)->wait();
    copy.memset(param_liveness_buffer, 1)->wait();

    // Number of tracks per seed
    vecmem::data::vector_buffer<unsigned int> n_tracks_per_seed_buffer(n_seeds,
                                                                       mr.main);
    copy.setup(n_tracks_per_seed_buffer)->wait();

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
    copy.setup(links_buffer)->wait();

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
    if (false && config.run_mbf_smoother) {
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
    copy.setup(tips_buffer)->wait();
    vecmem::data::vector_buffer<unsigned int> tip_length_buffer{
        config.max_num_branches_per_seed * n_seeds, mr.main};
    copy.setup(tip_length_buffer)->wait();

    vecmem::unique_alloc_ptr<bound_matrix<typename detector_t::algebra_type>[]>
        tmp_jacobian_ptr = nullptr;

    std::map<unsigned int, unsigned int> step_to_link_idx_map;
    step_to_link_idx_map[0] = 0;

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
        copy.setup(updated_params_buffer)->wait();

        vecmem::data::vector_buffer<unsigned int> updated_liveness_buffer(
            n_max_candidates, mr.main);
        copy.setup(updated_liveness_buffer)->wait();

        // Reset the number of tracks per seed
        copy.memset(n_tracks_per_seed_buffer, 0)->wait();

        const unsigned int links_size = copy.get_size(links_buffer);

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

            copy.setup(new_links_buffer)->wait();
            copy(links_buffer, new_links_buffer)->wait();

            links_buffer = std::move(new_links_buffer);
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
            // Now copy it to device memory.
            vecmem::data::vector_buffer<payload_t> device_payload(1u, mr.main);
            copy.setup(device_payload)->wait();
            copy(vecmem::data::vector_view<const payload_t>(1u, &host_payload),
                 device_payload)
                ->wait();

            // The number of threads to use per block in the track finding.
            const Idx blocksPerGrid =
                (n_in_params + threadsPerBlock - 1) / threadsPerBlock;
            const auto workDiv =
                makeWorkDiv<Acc>(blocksPerGrid, threadsPerBlock);

            // Submit the kernel to the queue.
            ::alpaka::exec<Acc>(queue, workDiv,
                                kernels::find_tracks<detector_t>{}, config,
                                device_payload.ptr());
            ::alpaka::wait(queue);

            std::swap(in_params_buffer, updated_params_buffer);
            std::swap(param_liveness_buffer, updated_liveness_buffer);

            step_to_link_idx_map[step + 1] = copy.get_size(links_buffer);
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
                const auto workDiv = makeWorkDiv<Acc>(nBlocks, nThreads);

                ::alpaka::exec<Acc>(
                    queue, workDiv,
                    kernels::fill_finding_duplicate_removal_sort_keys{},
                    device::fill_finding_duplicate_removal_sort_keys_payload{
                        .links_view = links_buffer,
                        .param_liveness_view = param_liveness_buffer,
                        .link_last_measurement_view =
                            link_last_measurement_buffer,
                        .param_ids_view = param_ids_buffer,
                        .n_links = n_candidates,
                        .curr_links_idx = step_to_link_idx_map[step],
                        .n_measurements = n_measurements});
                ::alpaka::wait(queue);
            }

            vecmem::device_vector<unsigned int> keys_device(
                link_last_measurement_buffer);
            vecmem::device_vector<unsigned int> param_ids_device(
                param_ids_buffer);
            details::sort_by_key(queue, mr, keys_device.begin(),
                                 keys_device.end(), param_ids_device.begin());

            /*
             * Then, we run the actual duplicate removal kernel.
             */
            {
                const unsigned int nThreads = 256;
                const unsigned int nBlocks =
                    (n_candidates + nThreads - 1) / nThreads;
                const auto workDiv = makeWorkDiv<Acc>(nBlocks, nThreads);

                ::alpaka::exec<Acc>(
                    queue, workDiv, kernels::remove_duplicates{}, config,
                    device::remove_duplicates_payload{
                        .links_view = links_buffer,
                        .link_last_measurement_view =
                            link_last_measurement_buffer,
                        .param_ids_view = param_ids_buffer,
                        .param_liveness_view = param_liveness_buffer,
                        .n_links = n_candidates,
                        .curr_links_idx = step_to_link_idx_map[step],
                        .n_measurements = n_measurements,
                        .step = step});
                ::alpaka::wait(queue);
            }
        }

        if (step == config.max_track_candidates_per_track - 1) {
            break;
        }

        if (n_candidates > 0) {
            /*****************************************************************
             * Kernel3: Get key and value for parameter sorting
             *****************************************************************/

            vecmem::data::vector_buffer<unsigned int> param_ids_buffer(
                n_candidates, mr.main);
            copy.setup(param_ids_buffer)->wait();

            {
                vecmem::data::vector_buffer<device::sort_key> keys_buffer(
                    n_candidates, mr.main);
                copy.setup(keys_buffer)->wait();

                Idx blocksPerGrid =
                    (n_candidates + threadsPerBlock - 1) / threadsPerBlock;
                auto workDiv = makeWorkDiv<Acc>(blocksPerGrid, threadsPerBlock);

                ::alpaka::exec<Acc>(
                    queue, workDiv,
                    kernels::fill_finding_propagation_sort_keys{},
                    device::fill_finding_propagation_sort_keys_payload{
                        in_params_buffer, param_liveness_buffer, keys_buffer,
                        param_ids_buffer});
                ::alpaka::wait(queue);

                // Sort the key and values
                vecmem::device_vector<device::sort_key> keys_device(
                    keys_buffer);
                vecmem::device_vector<unsigned int> param_ids_device(
                    param_ids_buffer);
                details::sort_by_key(queue, mr, keys_device.begin(),
                                     keys_device.end(),
                                     param_ids_device.begin());
            }

            /*****************************************************************
             * Kernel4: Propagate to the next surface
             *****************************************************************/

            {
                if (false && config.run_mbf_smoother) {
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
                // Now copy it to device memory.
                vecmem::data::vector_buffer<payload_t> device_payload(1u,
                                                                      mr.main);
                copy.setup(device_payload)->wait();
                copy(vecmem::data::vector_view<const payload_t>(1u,
                                                                &host_payload),
                     device_payload)
                    ->wait();

                // The number of threads to use per block in the propagation.
                const Idx blocksPerGrid =
                    (n_candidates + threadsPerBlock - 1) / threadsPerBlock;
                const auto workDiv =
                    makeWorkDiv<Acc>(blocksPerGrid, threadsPerBlock);

                // Launch the kernel to propagate all active tracks to the next
                // surface.
                ::alpaka::exec<Acc>(
                    queue, workDiv,
                    kernels::propagate_to_next_surface<
                        traccc::details::ckf_propagator_t<detector_t, bfield_t>,
                        bfield_t>{},
                    config, device_payload.ptr());
                ::alpaka::wait(queue);
            }
        }

        n_in_params = n_candidates;
    }

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
    const unsigned int n_tips_total = copy.get_size(tips_buffer);

    std::vector<unsigned int> tips_length_host;

    if (n_tips_total > 0) {
        copy(tip_length_buffer, tips_length_host)->wait();
        tips_length_host.resize(n_tips_total);
    }

    // Create track candidate buffer
    typename edm::track_container<typename detector_t::algebra_type>::buffer
        track_candidates_buffer{
            {tips_length_host, mr.main, mr.host}, {}, measurements_view};
    copy.setup(track_candidates_buffer.tracks)->wait();

    if (n_tips_total > 0) {
        const Idx blocksPerGrid =
            (n_tips_total + threadsPerBlock - 1) / threadsPerBlock;
        const auto workDiv = makeWorkDiv<Acc>(blocksPerGrid, threadsPerBlock);

        ::alpaka::exec<Acc>(
            queue, workDiv, kernels::build_tracks{},
            false && config.run_mbf_smoother, config.meas_calibration,
            device::build_tracks_payload{
                .seeds_view = seeds,
                .links_view = links_buffer,
                .tips_view = tips_buffer,
                .tracks_view = track_candidates_buffer,
                .tip_to_output_map = {},
                .jacobian_ptr = jacobian_ptr.get(),
                .link_predicted_parameter_view =
                    link_predicted_parameter_buffer,
                .link_filtered_parameter_view = link_filtered_parameter_buffer,
            });
        ::alpaka::wait(queue);
    }

    return track_candidates_buffer;
}

}  // namespace traccc::alpaka::details

namespace alpaka::trait {

/// Specify how much dynamic shared memory is needed for the
/// @c traccc::alpaka::details::kernels::find_tracks kernel.
template <typename TAcc, typename detector_t>
struct BlockSharedMemDynSizeBytes<
    traccc::alpaka::details::kernels::find_tracks<detector_t>, TAcc> {
    template <typename TVec, typename... TArgs>
    ALPAKA_FN_HOST_ACC static auto getBlockSharedMemDynSizeBytes(
        traccc::alpaka::details::kernels::find_tracks<
            detector_t> const& /* kernel */,
        TVec const& blockThreadExtent, TVec const& /* threadElemExtent */,
        TArgs const&... /* args */
        ) -> std::size_t {
        return static_cast<std::size_t>(blockThreadExtent.prod()) *
                   sizeof(unsigned long long int) +
               2 * static_cast<std::size_t>(blockThreadExtent.prod()) *
                   sizeof(std::pair<unsigned int, unsigned int>);
    }
};

}  // namespace alpaka::trait

namespace alpaka {

/// Convince Alpaka that
/// @c traccc::device::fill_finding_duplicate_removal_sort_keys_payload
/// is trivially copyable
template <>
struct IsKernelArgumentTriviallyCopyable<
    traccc::device::fill_finding_duplicate_removal_sort_keys_payload, void>
    : std::true_type {};

/// Convince Alpaka that
/// @c traccc::device::remove_duplicates_payload
/// is trivially copyable
template <>
struct IsKernelArgumentTriviallyCopyable<
    traccc::device::remove_duplicates_payload, void> : std::true_type {};

}  // namespace alpaka
