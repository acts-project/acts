/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// SYCL include(s).
#include <sycl/sycl.hpp>

// Local include(s).
#include "../sanity/contiguous_on.hpp"
#include "../utils/barrier.hpp"
#include "../utils/calculate1DimNdRange.hpp"
#include "../utils/global_index.hpp"
#include "../utils/oneDPL.hpp"
#include "../utils/thread_id.hpp"

// Project include(s).
#include "traccc/edm/device/identity_projector.hpp"
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/track_container.hpp"
#include "traccc/finding/actors/ckf_aborter.hpp"
#include "traccc/finding/candidate_link.hpp"
#include "traccc/finding/details/combinatorial_kalman_filter_types.hpp"
#include "traccc/finding/device/build_tracks.hpp"
#include "traccc/finding/device/fill_finding_duplicate_removal_sort_keys.hpp"
#include "traccc/finding/device/fill_finding_propagation_sort_keys.hpp"
#include "traccc/finding/device/find_tracks.hpp"
#include "traccc/finding/device/gather_best_tips_per_measurement.hpp"
#include "traccc/finding/device/gather_measurement_votes.hpp"
#include "traccc/finding/device/geo_id_surface_comparator.hpp"
#include "traccc/finding/device/propagate_to_next_surface.hpp"
#include "traccc/finding/device/remove_duplicates.hpp"
#include "traccc/finding/device/update_tip_length_buffer.hpp"
#include "traccc/finding/finding_config.hpp"
#include "traccc/utils/memory_resource.hpp"
#include "traccc/utils/projections.hpp"
#include "traccc/utils/propagation.hpp"

// VecMem include(s).
#include <vecmem/utils/copy.hpp>
#include <vecmem/utils/sycl/local_accessor.hpp>

namespace traccc::sycl::details {
namespace kernels {
template <typename T>
struct find_tracks {};
template <typename T>
struct fill_finding_duplicate_removal_sort_keys {};
template <typename T>
struct remove_duplicates {};
template <typename T>
struct fill_finding_propagation_sort_keys {};
template <typename T>
struct propagate_to_next_surface {};
template <typename T>
struct gather_best_tips_per_measurement {};

template <typename T>
struct gather_measurement_votes {};

template <typename T>
struct update_tip_length_buffer {};
template <typename T>
struct build_tracks {};

template <typename T>
struct upper_bound {};
template <typename T>
struct sort_by_key_1 {};
template <typename T>
struct sort_by_key_2 {};

}  // namespace kernels

/// Templated implementation of the track finding algorithm.
///
/// Concrete track finding algorithms can use this function with the appropriate
/// specializations, to find tracks on top of a specific detector type, magnetic
/// field type, and track finding configuration.
///
/// @tparam kernel_t   Structure to generate unique kernel names with
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
/// @param queue             The SYCL queue to use
///
/// @return A buffer of the found track candidates
///
template <typename kernel_t, typename detector_t, typename bfield_t>
edm::track_container<typename detector_t::algebra_type>::buffer
combinatorial_kalman_filter(
    const typename detector_t::const_view_type& det, const bfield_t& field,
    const edm::measurement_collection::const_view& measurements_view,
    const bound_track_parameters_collection_types::const_view& seeds,
    const finding_config& config, const memory_resource& mr, vecmem::copy& copy,
    ::sycl::queue& queue) {

    const edm::measurement_collection::const_device measurements{
        measurements_view};

    assert(config.min_step_length_for_next_surface >
               math::fabs(config.propagation.navigation.intersection
                              .overstep_tolerance) &&
           "Min step length for the next surface should be higher than the "
           "overstep tolerance");
    assert(is_contiguous_on<
           vecmem::device_vector<const detray::geometry::identifier>>(
        device::identity_projector{}, mr.main, copy, queue,
        measurements_view.template get<6>()));

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

    oneapi::dpl::upper_bound(
        oneapi::dpl::execution::device_policy<kernels::upper_bound<kernel_t>>{
            queue},
        measurements.surface_link().begin(),
        // We have to use this ugly form here, because if the
        // measurement collection is resizable (which it often
        // is), the end() function cannot be used in host code.
        measurements.surface_link().begin() + n_measurements,
        device_det.surfaces().begin(), device_det.surfaces().end(),
        measurement_ranges.begin(), device::geo_id_surface_comparator{});
    queue.wait_and_throw();

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
    copy.setup(tips_buffer)->wait();
    vecmem::data::vector_buffer<unsigned int> tip_length_buffer{
        config.max_num_branches_per_seed * n_seeds, mr.main};
    copy.setup(tip_length_buffer)->wait();

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

            link_buffer_capacity = new_link_buffer_capacity;

            vecmem::data::vector_buffer<candidate_link> new_links_buffer(
                link_buffer_capacity, mr.main,
                vecmem::data::buffer_type::resizable);

            copy.setup(new_links_buffer)->wait();
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

                copy(
                    vecmem::data::vector_view<
                        bound_matrix<typename detector_t::algebra_type>>{
                        links_size, jacobian_ptr.get()},
                    vecmem::data::vector_view<
                        bound_matrix<typename detector_t::algebra_type>>{
                        link_buffer_capacity, new_jacobian_ptr.get()})
                    ->wait();
                copy(link_predicted_parameter_buffer,
                     new_link_predicted_parameter_buffer)
                    ->wait();
                copy(link_filtered_parameter_buffer,
                     new_link_filtered_parameter_buffer)
                    ->wait();

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
            copy.setup(tmp_links_buffer)->wait();
            bound_track_parameters_collection_types::buffer tmp_params_buffer(
                n_max_candidates, mr.main);
            copy.setup(tmp_params_buffer)->wait();

            // The number of threads to use per block in the track finding.
            static const unsigned int nFindTracksThreads = 64;

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

            // Submit the kernel to the queue.
            queue
                .submit([&](::sycl::handler& h) {
                    // Allocate shared memory for the kernel.
                    vecmem::sycl::local_accessor<unsigned long long int>
                        shared_insertion_mutex(nFindTracksThreads, h);
                    vecmem::sycl::local_accessor<
                        std::pair<unsigned int, unsigned int>>
                        shared_candidates(2 * nFindTracksThreads, h);
                    vecmem::sycl::local_accessor<unsigned int>
                        shared_candidates_size(1, h);
                    vecmem::sycl::local_accessor<unsigned int>
                        shared_num_out_params(1, h);
                    vecmem::sycl::local_accessor<unsigned int>
                        shared_out_offset(1, h);
                    // Launch the kernel.
                    h.parallel_for<kernels::find_tracks<kernel_t>>(
                        calculate1DimNdRange(n_in_params, nFindTracksThreads),
                        [config, payload = device_payload.ptr(),
                         shared_insertion_mutex, shared_candidates,
                         shared_candidates_size, shared_num_out_params,
                         shared_out_offset](::sycl::nd_item<1> item) {
                            // SYCL wrappers used in the algorithm.
                            const details::barrier barrier{item};
                            const details::thread_id thread_id{item};

                            // Call the device function to find tracks.
                            device::find_tracks<detector_t>(
                                thread_id, barrier, config, *payload,
                                {shared_num_out_params[0], shared_out_offset[0],
                                 &(shared_insertion_mutex[0]),
                                 &(shared_candidates[0]),
                                 shared_candidates_size[0]});
                        });
                })
                .wait_and_throw();

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
            queue
                .submit([&](::sycl::handler& h) {
                    h.parallel_for<
                        kernels::fill_finding_duplicate_removal_sort_keys<
                            kernel_t>>(
                        calculate1DimNdRange(n_candidates, 256),
                        [links_view = vecmem::get_data(links_buffer),
                         param_liveness_view =
                             vecmem::get_data(param_liveness_buffer),
                         link_last_measurement_view =
                             vecmem::get_data(link_last_measurement_buffer),
                         param_ids_view = vecmem::get_data(param_ids_buffer),
                         n_candidates,
                         curr_links_idx = step_to_link_idx_map[step],
                         n_measurements](::sycl::nd_item<1> item) {
                            device::fill_finding_duplicate_removal_sort_keys(
                                details::global_index(item),
                                {.links_view = links_view,
                                 .param_liveness_view = param_liveness_view,
                                 .link_last_measurement_view =
                                     link_last_measurement_view,
                                 .param_ids_view = param_ids_view,
                                 .n_links = n_candidates,
                                 .curr_links_idx = curr_links_idx,
                                 .n_measurements = n_measurements});
                        });
                })
                .wait_and_throw();

            vecmem::device_vector<unsigned int> keys_device(
                link_last_measurement_buffer);
            vecmem::device_vector<unsigned int> param_ids_device(
                param_ids_buffer);
            oneapi::dpl::sort_by_key(
                oneapi::dpl::execution::device_policy<
                    kernels::sort_by_key_1<kernel_t>>{queue},
                keys_device.begin(), keys_device.end(),
                param_ids_device.begin());
            queue.wait_and_throw();

            /*
             * Then, we run the actual duplicate removal kernel.
             */
            queue
                .submit([&](::sycl::handler& h) {
                    h.parallel_for<kernels::remove_duplicates<kernel_t>>(
                        calculate1DimNdRange(n_candidates, 256),
                        [config, links_view = vecmem::get_data(links_buffer),
                         link_last_measurement_view =
                             vecmem::get_data(link_last_measurement_buffer),
                         param_ids_view = vecmem::get_data(param_ids_buffer),
                         param_liveness_view =
                             vecmem::get_data(param_liveness_buffer),
                         n_candidates,
                         curr_links_idx = step_to_link_idx_map[step],
                         n_measurements, step](::sycl::nd_item<1> item) {
                            device::remove_duplicates(
                                details::global_index(item), config,
                                {.links_view = links_view,
                                 .link_last_measurement_view =
                                     link_last_measurement_view,
                                 .param_ids_view = param_ids_view,
                                 .param_liveness_view = param_liveness_view,
                                 .n_links = n_candidates,
                                 .curr_links_idx = curr_links_idx,
                                 .n_measurements = n_measurements,
                                 .step = step});
                        });
                })
                .wait_and_throw();
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

                queue
                    .submit([&](::sycl::handler& h) {
                        h.parallel_for<
                            kernels::fill_finding_propagation_sort_keys<
                                kernel_t>>(
                            calculate1DimNdRange(n_candidates, 256),
                            [in_params = vecmem::get_data(in_params_buffer),
                             param_liveness =
                                 vecmem::get_data(param_liveness_buffer),
                             keys = vecmem::get_data(keys_buffer),
                             param_ids = vecmem::get_data(param_ids_buffer)](
                                ::sycl::nd_item<1> item) {
                                device::fill_finding_propagation_sort_keys(
                                    details::global_index(item),
                                    {in_params, param_liveness, keys,
                                     param_ids});
                            });
                    })
                    .wait_and_throw();

                // Sort the keys and values.
                vecmem::device_vector<device::sort_key> keys_device(
                    keys_buffer);
                vecmem::device_vector<unsigned int> param_ids_device(
                    param_ids_buffer);
                oneapi::dpl::sort_by_key(
                    oneapi::dpl::execution::device_policy<
                        kernels::sort_by_key_2<kernel_t>>{queue},
                    keys_device.begin(), keys_device.end(),
                    param_ids_device.begin());
                queue.wait_and_throw();
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
                // Now copy it to device memory.
                vecmem::data::vector_buffer<payload_t> device_payload(1u,
                                                                      mr.main);
                copy.setup(device_payload)->wait();
                copy(vecmem::data::vector_view<const payload_t>(1u,
                                                                &host_payload),
                     device_payload)
                    ->wait();

                // Launch the kernel to propagate all active tracks to the next
                // surface.
                queue
                    .submit([&](::sycl::handler& h) {
                        h.parallel_for<
                            kernels::propagate_to_next_surface<kernel_t>>(
                            calculate1DimNdRange(n_candidates, 64),
                            [config, payload = device_payload.ptr()](
                                ::sycl::nd_item<1> item) {
                                device::propagate_to_next_surface<
                                    traccc::details::ckf_propagator_t<
                                        detector_t, bfield_t>,
                                    bfield_t>(details::global_index(item),
                                              config, *payload);
                            });
                    })
                    .wait_and_throw();
            }
        }

        n_in_params = n_candidates;
    }

    /*****************************************************************
     * Kernel5: Build tracks
     *****************************************************************/

    // Get the number of tips
    const unsigned int n_tips_total = copy.get_size(tips_buffer);

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

            const device::gather_best_tips_per_measurement_payload<
                typename detector_t::algebra_type>
                payload{tips_buffer,
                        links_buffer,
                        measurements_view,
                        best_tips_per_measurement_insertion_mutex_buffer,
                        best_tips_per_measurement_index_buffer,
                        best_tips_per_measurement_pval_buffer,
                        config.max_num_tracks_per_measurement};
            queue
                .submit([&](::sycl::handler& h) {
                    h.parallel_for<
                        kernels::gather_best_tips_per_measurement<kernel_t>>(
                        calculate1DimNdRange(n_tips_total, 32),
                        [payload](::sycl::nd_item<1> item) {
                            device::gather_best_tips_per_measurement(
                                details::global_index(item),
                                details::barrier{item}, payload);
                        });
                })
                .wait_and_throw();
        }

        vecmem::data::vector_buffer<unsigned int> votes_per_tip_buffer(
            n_tips_total, mr.main);
        copy.setup(votes_per_tip_buffer)->wait();
        copy.memset(votes_per_tip_buffer, 0)->wait();

        {
            const device::gather_measurement_votes_payload payload{
                best_tips_per_measurement_insertion_mutex_buffer,
                best_tips_per_measurement_index_buffer, votes_per_tip_buffer,
                config.max_num_tracks_per_measurement};

            queue
                .submit([&](::sycl::handler& h) {
                    h.parallel_for<kernels::gather_measurement_votes<kernel_t>>(
                        calculate1DimNdRange(
                            config.max_num_tracks_per_measurement *
                                n_measurements,
                            512),
                        [payload](::sycl::nd_item<1> item) {
                            device::gather_measurement_votes(
                                details::global_index(item), payload);
                        });
                })
                .wait_and_throw();
        }

        tip_to_output_map =
            vecmem::data::vector_buffer<unsigned int>(n_tips_total, mr.main);
        copy.setup(tip_to_output_map)->wait();

        {
            vecmem::data::vector_buffer<unsigned int> new_tip_length_buffer{
                n_tips_total, mr.main, vecmem::data::buffer_type::resizable};
            copy.setup(new_tip_length_buffer)->wait();

            const device::update_tip_length_buffer_payload payload{
                tip_length_buffer, new_tip_length_buffer, votes_per_tip_buffer,
                tip_to_output_map, config.min_measurement_voting_fraction};

            queue
                .submit([&](::sycl::handler& h) {
                    h.parallel_for<kernels::update_tip_length_buffer<kernel_t>>(
                        calculate1DimNdRange(n_tips_total, 512),
                        [payload](::sycl::nd_item<1> item) {
                            device::update_tip_length_buffer(
                                details::global_index(item), payload);
                        });
                })
                .wait_and_throw();

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

    if (n_tips_total > 0) {
        copy(tip_length_buffer, tips_length_host)->wait();
        tips_length_host.resize(n_tips_total);
    }

    // Create track candidate buffer
    typename edm::track_container<typename detector_t::algebra_type>::buffer
        track_candidates_buffer{
            {tips_length_host, mr.main, mr.host},
            {n_states, mr.main, vecmem::data::buffer_type::resizable},
            measurements_view};
    copy.setup(track_candidates_buffer.tracks)->wait();
    copy.setup(track_candidates_buffer.states)->wait();

    if (n_tips_total > 0) {
        queue
            .submit([&](::sycl::handler& h) {
                h.parallel_for<kernels::build_tracks<kernel_t>>(
                    calculate1DimNdRange(n_tips_total, 64),
                    [config, seeds, links = vecmem::get_data(links_buffer),
                     tips = vecmem::get_data(tips_buffer),
                     tracks = typename edm::track_container<
                         typename detector_t::algebra_type>::
                         view(track_candidates_buffer),
                     tip_to_output_map = vecmem::get_data(tip_to_output_map),
                     jacobian_ptr = jacobian_ptr.get(),
                     link_predicted_parameters =
                         vecmem::get_data(link_predicted_parameter_buffer),
                     link_filtered_parameters =
                         vecmem::get_data(link_filtered_parameter_buffer)](
                        ::sycl::nd_item<1> item) {
                        device::build_tracks(
                            details::global_index(item),
                            config.run_mbf_smoother, config.meas_calibration,
                            {.seeds_view = seeds,
                             .links_view = links,
                             .tips_view = tips,
                             .tracks_view = tracks,
                             .tip_to_output_map = tip_to_output_map,
                             .jacobian_ptr = jacobian_ptr,
                             .link_predicted_parameter_view =
                                 link_predicted_parameters,
                             .link_filtered_parameter_view =
                                 link_filtered_parameters});
                    });
            })
            .wait_and_throw();
    }

    return track_candidates_buffer;
}

}  // namespace traccc::sycl::details
