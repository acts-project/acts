/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "../utils/cuda_error_handling.hpp"
#include "../utils/global_index.hpp"
#include "../utils/utils.hpp"
#include "./kernels/fill_fitting_sort_keys.hpp"
#include "./kernels/fit_backward.hpp"
#include "./kernels/fit_forward.hpp"
#include "./kernels/fit_prelude.hpp"

// Project include(s).
#include "traccc/edm/device/sort_key.hpp"
#include "traccc/edm/track_container.hpp"
#include "traccc/fitting/details/kalman_fitting_types.hpp"
#include "traccc/fitting/device/fill_fitting_sort_keys.hpp"
#include "traccc/fitting/fitting_config.hpp"
#include "traccc/utils/memory_resource.hpp"

// VecMem include(s).
#include <vecmem/utils/copy.hpp>

// Thrust include(s).
#include <thrust/execution_policy.h>
#include <thrust/sort.h>

// System include(s).
#include <numeric>

namespace traccc::cuda::details {

/// Templated implementation of the CUDA track fitting algorithm.
///
/// @tparam detector_t The (device) detector type to use
/// @tparam bfield_t   The magnetic field type to use
///
/// @param[in] det_view     A view of the detector geometry
/// @param[in] field_view   A view of the magnetic field
/// @param[in] track_candidates_view All track candidates to fit
/// @param[in] config       The fitting configuration
/// @param[in] mr           Memory resource(s) to use
/// @param[in] copy         The copy object to use for memory transfers
/// @param[in] queue        The Alpaka queue to use for execution
///
/// @return A container of the fitted track states
///
template <typename detector_t, typename bfield_t>
typename edm::track_container<typename detector_t::algebra_type>::buffer
kalman_fitting(
    const typename detector_t::const_view_type& det_view,
    const bfield_t& field_view,
    const typename edm::track_container<
        typename detector_t::algebra_type>::const_view& track_candidates_view,
    const fitting_config& config, const memory_resource& mr, vecmem::copy& copy,
    stream& str, unsigned int warp_size) {

    // Get a convenience variable for the stream that we'll be using.
    cudaStream_t stream = details::get_stream(str);

    // Get the number of tracks.
    const unsigned int n_tracks = copy.get_size(track_candidates_view.tracks);

    // Get the sizes of the track candidates in each track.
    const std::vector<unsigned int> candidate_sizes =
        copy.get_sizes(track_candidates_view.tracks);
    const unsigned int n_states =
        std::accumulate(candidate_sizes.begin(), candidate_sizes.end(), 0u);

    // Create the result buffer.
    typename edm::track_container<typename detector_t::algebra_type>::buffer
        track_states_buffer{
            {candidate_sizes, mr.main, mr.host,
             vecmem::data::buffer_type::resizable},
            {n_states, mr.main, vecmem::data::buffer_type::resizable},
            track_candidates_view.measurements};
    copy.setup(track_states_buffer.tracks)->ignore();
    copy.setup(track_states_buffer.states)->ignore();

    // Return early, if there are no tracks.
    if (n_tracks == 0) {
        return track_states_buffer;
    }

    std::vector<unsigned int> seqs_sizes(candidate_sizes.size());
    std::transform(candidate_sizes.begin(), candidate_sizes.end(),
                   seqs_sizes.begin(), [&config](const unsigned int sz) {
                       return std::max(sz * config.surface_sequence_size_factor,
                                       config.min_surface_sequence_capacity);
                   });
    vecmem::data::jagged_vector_buffer<typename detector_t::surface_type>
        seqs_buffer{seqs_sizes, mr.main, mr.host,
                    vecmem::data::buffer_type::resizable};
    copy.setup(seqs_buffer)->ignore();

    // Create the buffers for sorting the parameter IDs.
    vecmem::data::vector_buffer<device::sort_key> keys_buffer(n_tracks,
                                                              mr.main);
    vecmem::data::vector_buffer<unsigned int> param_ids_buffer(n_tracks,
                                                               mr.main);
    vecmem::data::vector_buffer<unsigned int> param_liveness_buffer(n_tracks,
                                                                    mr.main);
    vecmem::copy::event_type keys_setup_event = copy.setup(keys_buffer);
    vecmem::copy::event_type param_ids_setup_event =
        copy.setup(param_ids_buffer);
    vecmem::copy::event_type param_liveness_setup_event =
        copy.setup(param_liveness_buffer);
    keys_setup_event->ignore();
    param_ids_setup_event->ignore();
    param_liveness_setup_event->ignore();

    // Launch parameters for all the kernels.
    const unsigned int nThreads = warp_size * 4;
    const unsigned int nBlocks = (n_tracks + nThreads - 1) / nThreads;

    // Fill the keys and param_ids buffers.
    fill_fitting_sort_keys(nBlocks, nThreads, stream,
                           track_candidates_view.tracks, keys_buffer,
                           param_ids_buffer);

    // Sort the key to get the sorted parameter ids
    vecmem::device_vector<device::sort_key> keys_device(keys_buffer);
    vecmem::device_vector<unsigned int> param_ids_device(param_ids_buffer);
    thrust::sort_by_key(
        thrust::cuda::par_nosync(std::pmr::polymorphic_allocator(&mr.main))
            .on(stream),
        keys_device.begin(), keys_device.end(), param_ids_device.begin());

    // Run the fitting, using the sorted parameter IDs.
    fit_prelude(nBlocks, nThreads, 0, stream, param_ids_buffer,
                track_candidates_view, track_states_buffer,
                param_liveness_buffer);
    TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
    str.synchronize();

    // Allocate the fitting kernels's payload in host memory.
    using fitter_t = traccc::details::kalman_fitter_t<detector_t, bfield_t>;
    device::fit_payload<fitter_t> host_payload{
        .det_data = det_view,
        .field_data = field_view,
        .param_ids_view = param_ids_buffer,
        .param_liveness_view = param_liveness_buffer,
        .tracks_view = track_states_buffer,
        .surfaces_view = seqs_buffer};

    for (std::size_t i = 0; i < config.n_iterations; ++i) {
        // Run the track fitting
        fit_forward<fitter_t>(nBlocks, nThreads, 0, stream, config,
                              host_payload);
        TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
        fit_backward<fitter_t>(nBlocks, nThreads, 0, stream, config,
                               host_payload);
        TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
    }

    // Return the fitted tracks.
    return track_states_buffer;
}

}  // namespace traccc::cuda::details
