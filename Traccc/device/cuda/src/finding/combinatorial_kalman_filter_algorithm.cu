/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "../sanity/contiguous_on.cuh"
#include "../utils/magnetic_field_types.hpp"
#include "./kernels/build_tracks.cuh"
#include "./kernels/condense_tracks.cuh"
#include "./kernels/fill_finding_duplicate_removal_sort_keys.cuh"
#include "./kernels/fill_finding_propagation_sort_keys.cuh"
#include "./kernels/find_tracks.cuh"
#include "./kernels/gather_best_tips_per_measurement.cuh"
#include "./kernels/gather_measurement_votes.cuh"
#include "./kernels/progressive_kalman_filter.hpp"
#include "./kernels/propagate_to_next_surface.hpp"
#include "./kernels/remove_duplicates.cuh"
#include "./kernels/update_tip_length_buffer.cuh"
#include "traccc/cuda/finding/combinatorial_kalman_filter_algorithm.hpp"

// Project include(s).
#include "traccc/bfield/magnetic_field_types.hpp"
#include "traccc/edm/device/identity_projector.hpp"
#include "traccc/finding/details/combinatorial_kalman_filter_types.hpp"
#include "traccc/finding/device/geo_id_surface_comparator.hpp"
#include "traccc/geometry/detector_buffer.hpp"
#include "traccc/utils/detector_buffer_bfield_visitor.hpp"

// Thrust include(s).
#include <thrust/execution_policy.h>
#include <thrust/scan.h>
#include <thrust/sort.h>

namespace traccc::cuda {

bool combinatorial_kalman_filter_algorithm::input_is_valid(
    const edm::measurement_collection::const_view& measurements) const {

    static constexpr std::size_t GEOMID_INDEX = 6u;
    return is_contiguous_on<
        vecmem::device_vector<const detray::geometry::identifier>>(
        device::identity_projector{}, mr().main, copy(), stream(),
        measurements.template get<GEOMID_INDEX>());
}

vecmem::data::vector_buffer<edm::measurement_collection::const_view::size_type>
combinatorial_kalman_filter_algorithm::build_measurement_ranges_buffer(
    const detector_buffer& det,
    const edm::measurement_collection::const_view::size_type n_measurements,
    const edm::measurement_collection::const_view& measurements) const {

    return detector_buffer_visitor<detector_type_list>(
        det, [&]<typename detector_traits_t>(
                 const typename detector_traits_t::view& det) {
            // Construct an appropriate device detector object.
            typename detector_traits_t::device device_det{det};

            // Create the result buffer.
            vecmem::data::vector_buffer<
                edm::measurement_collection::const_view::size_type>
                result{device_det.surfaces().size(), mr().main};
            copy().setup(result)->ignore();

            // Create a measurement device object for convenience.
            const edm::measurement_collection::const_device measurements_device{
                measurements};

            // Fill it with Thrust's help.
            thrust::upper_bound(
                thrust::cuda::par_nosync(
                    std::pmr::polymorphic_allocator(&(mr().main)))
                    .on(details::get_stream(stream())),
                measurements_device.surface_link().begin(),
                // We have to use this ugly form here, because if the
                // measurement collection is resizable (which it often
                // is), the end() function cannot be used in host code.
                measurements_device.surface_link().begin() + n_measurements,
                device_det.surfaces().begin(), device_det.surfaces().end(),
                result.ptr(), device::geo_id_surface_comparator{});

            // Return the filled buffer.
            return result;
        });
}

void combinatorial_kalman_filter_algorithm::progressive_kalman_filter_kernel(
    unsigned int n_seeds, const finding_config& config,
    const detector_buffer& detector, const magnetic_field& field,
    const device::progressive_kalman_filter_payload& payload,
    const device::kalman_fitting_algorithm::fit_payload& smoothing_payload)
    const {

    // Establish the kernel launch parameters.
    const unsigned int deviceThreads = warp_size() * 4;
    const unsigned int deviceBlocks =
        (n_seeds + deviceThreads - 1) / deviceThreads;

    detector_buffer_magnetic_field_visitor<detector_type_list,
                                           cuda::bfield_type_list<scalar>>(
        detector, field,
        [&]<typename detector_traits_t, typename bfield_view_t>(
            const typename detector_traits_t::view& det,
            const bfield_view_t& bfield) {
            using detector_t = typename detector_traits_t::device;
            using surface_t = typename detector_t::surface_type;

            // If the Kalman smoother should be run, obtain the real allocation
            vecmem::data::jagged_vector_view<surface_t> sf_sequences;
            if (config.run_smoother == smoother_type::e_kalman) {
                sf_sequences =
                    smoothing_payload.surfaces
                        .as<vecmem::data::jagged_vector_buffer<surface_t>>();
            }

            progressive_kalman_filter<
                traccc::details::pkf_propagator_t<detector_t, bfield_view_t>>(
                deviceBlocks, deviceThreads, 0u, details::get_stream(stream()),
                config, det, bfield, sf_sequences, payload);
        });

    TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void combinatorial_kalman_filter_algorithm::find_tracks_kernel(
    unsigned int n_threads, const finding_config& config,
    const detector_buffer& detector,
    const device::find_tracks_payload& payload) const {

    // Establish the kernel launch parameters.
    const unsigned int deviceThreads = warp_size() * 2;
    const unsigned int deviceBlocks =
        (n_threads + deviceThreads - 1) / deviceThreads;
    const std::size_t deviceSharedMem =
        deviceThreads * sizeof(unsigned long long int) +
        2 * deviceThreads * sizeof(std::pair<unsigned int, unsigned int>);

    // Launch the kernel for the appropriate detector type.
    detector_buffer_visitor<detector_type_list>(
        detector, [&]<typename detector_traits_t>(
                      const typename detector_traits_t::view& det) {
            find_tracks<typename detector_traits_t::device>(
                deviceBlocks, deviceThreads, deviceSharedMem,
                details::get_stream(stream()), config, det, payload);
        });
    TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void combinatorial_kalman_filter_algorithm::condense_tracks_kernel(
    unsigned int n_threads,
    const vecmem::data::vector_view<const unsigned int>&
        out_params_per_in_param,
    vecmem::data::vector_view<unsigned int>& params_index,
    const device::condense_tracks_payload& payload) const {

    // Prepare the "in_params_index_view" input parameter of the track
    // condensing kernel, with Thrust's help.
    const vecmem::device_vector<const unsigned int>
        out_params_per_in_param_vector(out_params_per_in_param);
    vecmem::device_vector<unsigned int> params_index_vector(params_index);
    thrust::inclusive_scan(
        thrust::cuda::par_nosync(std::pmr::polymorphic_allocator(&(mr().main)))
            .on(details::get_stream(stream())),
        out_params_per_in_param_vector.begin(),
        out_params_per_in_param_vector.end(), params_index_vector.begin());

    // Establish the kernel launch parameters.
    const unsigned int deviceThreads = warp_size() * 8;
    const unsigned int deviceBlocks =
        (n_threads + deviceThreads - 1) / deviceThreads;

    // Launch the kernel.
    assert(params_index.ptr() == payload.in_params_index_view.ptr());
    condense_tracks(deviceBlocks, deviceThreads, 0,
                    details::get_stream(stream()), payload);
    TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void combinatorial_kalman_filter_algorithm::
    fill_finding_duplicate_removal_sort_keys_kernel(
        unsigned int n_threads,
        const device::fill_finding_duplicate_removal_sort_keys_payload& payload)
        const {

    // Establish the kernel launch parameters.
    const unsigned int deviceThreads = warp_size() * 8;
    const unsigned int deviceBlocks =
        (n_threads + deviceThreads - 1) / deviceThreads;

    // Launch the kernel.
    kernels::fill_finding_duplicate_removal_sort_keys<<<
        deviceBlocks, deviceThreads, 0, details::get_stream(stream())>>>(
        payload);
    TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void combinatorial_kalman_filter_algorithm::sort_param_ids_by_last_measurement(
    vecmem::data::vector_view<unsigned int>& link_last_measurement,
    vecmem::data::vector_view<unsigned int>& param_ids) const {

    assert(link_last_measurement.capacity() == param_ids.capacity());
    assert(link_last_measurement.size_ptr() == nullptr);
    assert(param_ids.size_ptr() == nullptr);
    thrust::sort_by_key(
        thrust::cuda::par_nosync(std::pmr::polymorphic_allocator(&(mr().main)))
            .on(details::get_stream(stream())),
        link_last_measurement.ptr(),
        link_last_measurement.ptr() + link_last_measurement.capacity(),
        param_ids.ptr());
}

void combinatorial_kalman_filter_algorithm::remove_duplicates_kernel(
    unsigned int n_threads, const finding_config& config,
    const device::remove_duplicates_payload& payload) const {

    // Establish the kernel launch parameters.
    const unsigned int deviceThreads = warp_size() * 8;
    const unsigned int deviceBlocks =
        (n_threads + deviceThreads - 1) / deviceThreads;

    // Launch the kernel.
    kernels::remove_duplicates<<<deviceBlocks, deviceThreads, 0,
                                 details::get_stream(stream())>>>(config,
                                                                  payload);
    TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void combinatorial_kalman_filter_algorithm::
    fill_finding_propagation_sort_keys_kernel(
        unsigned int n_threads,
        const device::fill_finding_propagation_sort_keys_payload& payload)
        const {

    // Establish the kernel launch parameters.
    const unsigned int deviceThreads = warp_size() * 2;
    const unsigned int deviceBlocks =
        (n_threads + deviceThreads - 1) / deviceThreads;

    // Launch the kernel.
    kernels::fill_finding_propagation_sort_keys<<<
        deviceBlocks, deviceThreads, 0, details::get_stream(stream())>>>(
        payload);
    TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void combinatorial_kalman_filter_algorithm::sort_param_ids_by_keys(
    vecmem::data::vector_view<device::sort_key>& keys,
    vecmem::data::vector_view<unsigned int>& param_ids) const {

    assert(keys.capacity() == param_ids.capacity());
    assert(keys.size_ptr() == nullptr);
    assert(param_ids.size_ptr() == nullptr);
    thrust::sort_by_key(
        thrust::cuda::par_nosync(std::pmr::polymorphic_allocator(&(mr().main)))
            .on(details::get_stream(stream())),
        keys.ptr(), keys.ptr() + keys.capacity(), param_ids.ptr());
}

void combinatorial_kalman_filter_algorithm::propagate_to_next_surface_kernel(
    unsigned int n_threads, const finding_config& config,
    const detector_buffer& detector, const magnetic_field& field,
    const device::propagate_to_next_surface_payload& payload) const {

    // Establish the kernel launch parameters.
    const unsigned int deviceThreads = warp_size() * 4;
    const unsigned int deviceBlocks =
        (n_threads + deviceThreads - 1) / deviceThreads;

    // Launch the kernel for the appropriate detector and magnetic field type.
    detector_buffer_magnetic_field_visitor<detector_type_list,
                                           cuda::bfield_type_list<scalar>>(
        detector, field,
        [&]<typename detector_traits_t, typename bfield_view_t>(
            const typename detector_traits_t::view& det,
            const bfield_view_t& bfield) {
            propagate_to_next_surface<
                traccc::details::ckf_propagator_t<
                    typename detector_traits_t::device, bfield_view_t>,
                bfield_view_t>(deviceBlocks, deviceThreads, 0u,
                               details::get_stream(stream()), config, det,
                               bfield, payload);
        });
    TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void combinatorial_kalman_filter_algorithm::
    gather_best_tips_per_measurement_kernel(
        unsigned int n_threads,
        const device::gather_best_tips_per_measurement_payload<default_algebra>&
            payload) const {

    // Establish the kernel launch parameters.
    const unsigned int deviceThreads = warp_size();
    const unsigned int deviceBlocks =
        (n_threads + deviceThreads - 1) / deviceThreads;

    // Launch the kernel.
    kernels::gather_best_tips_per_measurement<<<
        deviceBlocks, deviceThreads, 0, details::get_stream(stream())>>>(
        payload);
    TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void combinatorial_kalman_filter_algorithm::gather_measurement_votes_kernel(
    unsigned int n_threads,
    const device::gather_measurement_votes_payload& payload) const {

    // Establish the kernel launch parameters.
    const unsigned int deviceThreads = warp_size() * 16;
    const unsigned int deviceBlocks =
        (n_threads + deviceThreads - 1) / deviceThreads;

    // Launch the kernel.
    kernels::gather_measurement_votes<<<deviceBlocks, deviceThreads, 0,
                                        details::get_stream(stream())>>>(
        payload);
    TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void combinatorial_kalman_filter_algorithm::update_tip_length_buffer_kernel(
    unsigned int n_threads,
    const device::update_tip_length_buffer_payload& payload) const {

    // Establish the kernel launch parameters.
    const unsigned int deviceThreads = warp_size() * 16;
    const unsigned int deviceBlocks =
        (n_threads + deviceThreads - 1) / deviceThreads;

    // Launch the kernel.
    kernels::update_tip_length_buffer<<<deviceBlocks, deviceThreads, 0,
                                        details::get_stream(stream())>>>(
        payload);
    TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void combinatorial_kalman_filter_algorithm::build_tracks_kernel(
    unsigned int n_threads, bool run_mbf_smoother,
    const measurement_selector::config& calib_cfg,
    const device::build_tracks_payload& payload) const {

    // Establish the kernel launch parameters.
    const unsigned int deviceThreads = warp_size() * 2;
    const unsigned int deviceBlocks =
        (n_threads + deviceThreads - 1) / deviceThreads;

    // Launch the kernel.
    kernels::build_tracks<<<deviceBlocks, deviceThreads, 0,
                            details::get_stream(stream())>>>(
        run_mbf_smoother, calib_cfg, payload);
    TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

}  // namespace traccc::cuda
