/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic warning "-Wmaybe-uninitialized"
#endif

// Local include(s).
#include "traccc/alpaka/finding/combinatorial_kalman_filter_algorithm.hpp"

#include "../utils/barrier.hpp"
#include "../utils/get_queue.hpp"
#include "../utils/magnetic_field_types.hpp"
#include "../utils/parallel_algorithms.hpp"
#include "../utils/thread_id.hpp"
#include "../utils/utils.hpp"

// Project include(s).
#include "traccc/bfield/magnetic_field_types.hpp"
#include "traccc/finding/details/combinatorial_kalman_filter_types.hpp"
#include "traccc/finding/device/build_tracks.hpp"
#include "traccc/finding/device/fill_finding_duplicate_removal_sort_keys.hpp"
#include "traccc/finding/device/fill_finding_propagation_sort_keys.hpp"
#include "traccc/finding/device/find_tracks.hpp"
#include "traccc/finding/device/geo_id_surface_comparator.hpp"
#include "traccc/finding/device/progressive_kalman_filter.hpp"
#include "traccc/finding/device/propagate_to_next_surface.hpp"
#include "traccc/finding/device/remove_duplicates.hpp"
#include "traccc/utils/detector_buffer_bfield_visitor.hpp"

namespace traccc::alpaka {
namespace kernels {

/// Alpaka kernel functor for @c traccc::device::find_tracks
template <typename detector_t>
struct find_tracks {
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc, const finding_config& cfg,
        const typename detector_t::const_view_type* det_data,
        const device::find_tracks_payload& payload) const {

        auto& shared_num_out_params =
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
            thread_id, barrier, cfg, *det_data, payload,
            {.shared_num_out_params = shared_num_out_params,
             .shared_insertion_mutex = shared_insertion_mutex,
             .shared_candidates = shared_candidates,
             .shared_candidates_size = shared_candidates_size});
    }
};

/// Alpaka kernel functor for @c traccc::device::condense_tracks
struct condense_tracks {
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc, const device::condense_tracks_payload& payload) const {

        const device::global_index_t globalThreadIdx =
            ::alpaka::getIdx<::alpaka::Grid, ::alpaka::Threads>(acc)[0];
        device::condense_tracks(globalThreadIdx, payload);
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
        const device::fill_finding_propagation_sort_keys_payload& payload)
        const {

        const device::global_index_t globalThreadIdx =
            ::alpaka::getIdx<::alpaka::Grid, ::alpaka::Threads>(acc)[0];
        device::fill_finding_propagation_sort_keys(globalThreadIdx, payload);
    }
};

/// Alpaka kernel functor for @c traccc::device::progressive_kalman_filter
template <typename propagator_t>
struct progressive_kalman_filter {
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc, const finding_config& cfg,
        const typename propagator_t::detector_type::const_view_type* det_data,
        const typename propagator_t::stepper_type::magnetic_field_type&
            field_data,
        vecmem::data::jagged_vector_view<
            typename propagator_t::detector_type::surface_type>
            sf_sequences,
        const device::progressive_kalman_filter_payload& payload) const {

        device::global_index_t globalThreadIdx =
            ::alpaka::getIdx<::alpaka::Grid, ::alpaka::Threads>(acc)[0];
        device::progressive_kalman_filter<propagator_t>(
            globalThreadIdx, cfg, *det_data, field_data, sf_sequences, payload);
    }
};

/// Alpaka kernel functor for @c traccc::device::propagate_to_next_surface
template <typename propagator_t, typename bfield_t>
struct propagate_to_next_surface {
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc, const finding_config& cfg,
        const typename propagator_t::detector_type::const_view_type* det_data,
        const bfield_t& field_data,
        const device::propagate_to_next_surface_payload& payload) const {

        device::global_index_t globalThreadIdx =
            ::alpaka::getIdx<::alpaka::Grid, ::alpaka::Threads>(acc)[0];
        device::propagate_to_next_surface<propagator_t, bfield_t>(
            globalThreadIdx, cfg, *det_data, field_data, payload);
    }
};

/// Alpaka kernel functor for
/// @c traccc::device::fill_finding_propagation_sort_keys
struct gather_best_tips_per_measurement {
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc,
        const device::gather_best_tips_per_measurement_payload<default_algebra>&
            payload) const {

        const device::global_index_t globalThreadIdx =
            ::alpaka::getIdx<::alpaka::Grid, ::alpaka::Threads>(acc)[0];
        device::gather_best_tips_per_measurement(
            globalThreadIdx, alpaka::barrier<TAcc>{&acc}, payload);
    }
};

/// Alpaka kernel functor for @c traccc::device::gather_measurement_votes
struct gather_measurement_votes {
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc,
        const device::gather_measurement_votes_payload& payload) const {

        const device::global_index_t globalThreadIdx =
            ::alpaka::getIdx<::alpaka::Grid, ::alpaka::Threads>(acc)[0];
        device::gather_measurement_votes(globalThreadIdx, payload);
    }
};

/// Alpaka kernel functor for @c traccc::device::update_tip_length_buffer
struct update_tip_length_buffer {
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc,
        const device::update_tip_length_buffer_payload& payload) const {

        const device::global_index_t globalThreadIdx =
            ::alpaka::getIdx<::alpaka::Grid, ::alpaka::Threads>(acc)[0];
        device::update_tip_length_buffer(globalThreadIdx, payload);
    }
};

/// Alpaka kernel functor for @c traccc::device::build_tracks
struct build_tracks {
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc, bool run_mbf,
        const measurement_selector::config& calib_cfg,
        const device::build_tracks_payload& payload) const {

        device::global_index_t globalThreadIdx =
            ::alpaka::getIdx<::alpaka::Grid, ::alpaka::Threads>(acc)[0];
        device::build_tracks(globalThreadIdx, run_mbf, calib_cfg, payload);
    }
};

}  // namespace kernels

combinatorial_kalman_filter_algorithm::combinatorial_kalman_filter_algorithm(
    const config_type& config, const traccc::memory_resource& mr,
    const vecmem::copy& copy, alpaka::queue& q,
    std::unique_ptr<const Logger> logger,
    std::unique_ptr<traccc::alpaka::kalman_fitting_algorithm> kf_fitter)
    : device::combinatorial_kalman_filter_algorithm(
          config, mr, copy, std::move(logger), std::move(kf_fitter)),
      alpaka::algorithm_base(q) {}

bool combinatorial_kalman_filter_algorithm::input_is_valid(
    const edm::measurement_collection::const_view&) const {

    // TODO: Implement sanity check(s).
    return true;
}

vecmem::data::vector_buffer<edm::measurement_collection::const_view::size_type>
combinatorial_kalman_filter_algorithm::build_measurement_ranges_buffer(
    const detector_buffer& detector,
    const edm::measurement_collection::const_view::size_type n_measurements,
    const edm::measurement_collection::const_view& measurements) const {

    return detector_buffer_visitor<detector_type_list>(
        detector, [&]<typename detector_traits_t>(
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
            details::upper_bound(
                details::get_queue(queue()), mr(),
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

    // Launch the kernel for the appropriate detector and magnetic field type.
    detector_buffer_magnetic_field_visitor<detector_type_list,
                                           alpaka::bfield_type_list<scalar>>(
        detector, field,
        [&]<typename detector_traits_t, typename bfield_view_t>(
            const typename detector_traits_t::view& det,
            const bfield_view_t& bfield) {
            using detector_t = typename detector_traits_t::device;
            using surface_t = typename detector_t::surface_type;

            // Copy the detector data to device memory.
            vecmem::data::vector_buffer<typename detector_traits_t::view>
                device_det(1u, mr().main);
            copy().setup(device_det)->ignore();
            copy()({1u, &det}, device_det)->ignore();

            // If the Kalman smoother should be run, obtain the real allocation
            vecmem::data::jagged_vector_view<surface_t> sf_sequences;
            if (config.run_smoother == smoother_type::e_kalman) {
                sf_sequences =
                    smoothing_payload.surfaces
                        .as<vecmem::data::jagged_vector_buffer<surface_t>>();
            }

            // Launch the kernel to propagate all active tracks to the next
            // surface.
            ::alpaka::exec<Acc>(
                details::get_queue(queue()),
                makeWorkDiv<Acc>(deviceBlocks, deviceThreads),
                kernels::progressive_kalman_filter<
                    traccc::details::ckf_propagator_t<detector_t,
                                                      bfield_view_t>>{},
                config, device_det.ptr(), bfield, sf_sequences, payload);
        });
}

void combinatorial_kalman_filter_algorithm::find_tracks_kernel(
    unsigned int n_threads, const finding_config& config,
    const detector_buffer& detector,
    const device::find_tracks_payload& payload) const {

    // Establish the kernel launch parameters.
    const unsigned int deviceThreads = warp_size() * 2;
    const unsigned int deviceBlocks =
        (n_threads + deviceThreads - 1) / deviceThreads;

    // Launch the kernel for the appropriate detector type.
    detector_buffer_visitor<detector_type_list>(
        detector, [&]<typename detector_traits_t>(
                      const typename detector_traits_t::view& det) {
            // Copy the detector data to device memory.
            vecmem::data::vector_buffer<typename detector_traits_t::view>
                device_det(1u, mr().main);
            copy().setup(device_det)->ignore();
            copy()({1u, &det}, device_det)->ignore();

            // Submit the kernel to the queue.
            ::alpaka::exec<Acc>(
                details::get_queue(queue()),
                makeWorkDiv<Acc>(deviceBlocks, deviceThreads),
                kernels::find_tracks<typename detector_traits_t::device>{},
                config, device_det.ptr(), payload);
        });
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
    details::inclusive_scan(details::get_queue(queue()), mr(),
                            out_params_per_in_param_vector.begin(),
                            out_params_per_in_param_vector.end(),
                            params_index_vector.begin());

    // Establish the kernel launch parameters.
    const unsigned int deviceThreads = warp_size() * 8;
    const unsigned int deviceBlocks =
        (n_threads + deviceThreads - 1) / deviceThreads;

    // Launch the kernel.
    assert(params_index.ptr() == payload.in_params_index_view.ptr());
    ::alpaka::exec<Acc>(details::get_queue(queue()),
                        makeWorkDiv<Acc>(deviceBlocks, deviceThreads),
                        kernels::condense_tracks{}, payload);
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
    ::alpaka::exec<Acc>(details::get_queue(queue()),
                        makeWorkDiv<Acc>(deviceBlocks, deviceThreads),
                        kernels::fill_finding_duplicate_removal_sort_keys{},
                        payload);
}

void combinatorial_kalman_filter_algorithm::sort_param_ids_by_last_measurement(
    vecmem::data::vector_view<unsigned int>& link_last_measurement,
    vecmem::data::vector_view<unsigned int>& param_ids) const {

    assert(link_last_measurement.capacity() == param_ids.capacity());
    assert(link_last_measurement.size_ptr() == nullptr);
    assert(param_ids.size_ptr() == nullptr);
    details::sort_by_key(
        details::get_queue(queue()), mr(), link_last_measurement.ptr(),
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
    ::alpaka::exec<Acc>(details::get_queue(queue()),
                        makeWorkDiv<Acc>(deviceBlocks, deviceThreads),
                        kernels::remove_duplicates{}, config, payload);
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
    ::alpaka::exec<Acc>(details::get_queue(queue()),
                        makeWorkDiv<Acc>(deviceBlocks, deviceThreads),
                        kernels::fill_finding_propagation_sort_keys{}, payload);
}

void combinatorial_kalman_filter_algorithm::sort_param_ids_by_keys(
    vecmem::data::vector_view<device::sort_key>& keys,
    vecmem::data::vector_view<unsigned int>& param_ids) const {

    assert(keys.capacity() == param_ids.capacity());
    assert(keys.size_ptr() == nullptr);
    assert(param_ids.size_ptr() == nullptr);
    details::sort_by_key(details::get_queue(queue()), mr(), keys.ptr(),
                         keys.ptr() + keys.capacity(), param_ids.ptr());
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
                                           alpaka::bfield_type_list<scalar>>(
        detector, field,
        [&]<typename detector_traits_t, typename bfield_view_t>(
            const typename detector_traits_t::view& det,
            const bfield_view_t& bfield) {
            // Copy the detector data to device memory.
            vecmem::data::vector_buffer<typename detector_traits_t::view>
                device_det(1u, mr().main);
            copy().setup(device_det)->ignore();
            copy()({1u, &det}, device_det)->ignore();

            // Launch the kernel to propagate all active tracks to the next
            // surface.
            ::alpaka::exec<Acc>(
                details::get_queue(queue()),
                makeWorkDiv<Acc>(deviceBlocks, deviceThreads),
                kernels::propagate_to_next_surface<
                    traccc::details::ckf_propagator_t<
                        typename detector_traits_t::device, bfield_view_t>,
                    bfield_view_t>{},
                config, device_det.ptr(), bfield, payload);
        });
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
    ::alpaka::exec<Acc>(details::get_queue(queue()),
                        makeWorkDiv<Acc>(deviceBlocks, deviceThreads),
                        kernels::gather_best_tips_per_measurement{}, payload);
}

void combinatorial_kalman_filter_algorithm::gather_measurement_votes_kernel(
    unsigned int n_threads,
    const device::gather_measurement_votes_payload& payload) const {

    // Establish the kernel launch parameters.
    const unsigned int deviceThreads = warp_size() * 16;
    const unsigned int deviceBlocks =
        (n_threads + deviceThreads - 1) / deviceThreads;

    // Launch the kernel.
    ::alpaka::exec<Acc>(details::get_queue(queue()),
                        makeWorkDiv<Acc>(deviceBlocks, deviceThreads),
                        kernels::gather_measurement_votes{}, payload);
}

void combinatorial_kalman_filter_algorithm::update_tip_length_buffer_kernel(
    unsigned int n_threads,
    const device::update_tip_length_buffer_payload& payload) const {

    // Establish the kernel launch parameters.
    const unsigned int deviceThreads = warp_size() * 16;
    const unsigned int deviceBlocks =
        (n_threads + deviceThreads - 1) / deviceThreads;

    // Launch the kernel.
    ::alpaka::exec<Acc>(details::get_queue(queue()),
                        makeWorkDiv<Acc>(deviceBlocks, deviceThreads),
                        kernels::update_tip_length_buffer{}, payload);
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
    ::alpaka::exec<Acc>(details::get_queue(queue()),
                        makeWorkDiv<Acc>(deviceBlocks, deviceThreads),
                        kernels::build_tracks{}, run_mbf_smoother, calib_cfg,
                        payload);
}

}  // namespace traccc::alpaka

namespace alpaka::trait {

/// Specify how much dynamic shared memory is needed for the
/// @c traccc::alpaka::details::kernels::find_tracks kernel.
template <typename TAcc, typename detector_t>
struct BlockSharedMemDynSizeBytes<
    traccc::alpaka::kernels::find_tracks<detector_t>, TAcc> {
    template <typename TVec, typename... TArgs>
    ALPAKA_FN_HOST_ACC static auto getBlockSharedMemDynSizeBytes(
        traccc::alpaka::kernels::find_tracks<detector_t> const& /* kernel */,
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
