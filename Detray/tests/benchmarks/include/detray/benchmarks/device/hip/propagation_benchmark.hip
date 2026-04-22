/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "detray/benchmarks/device/hip/propagation_benchmark.hpp"
#include "detray/definitions/detail/hip_definitions.hpp"

namespace detray::benchmarks {

template <typename propagator_t, detray::benchmarks::propagation_opt kOPT>
__global__ void __launch_bounds__(256, 4) propagator_benchmark_kernel(
    propagation::config cfg,
    typename propagator_t::detector_type::view_type det_view,
    typename propagator_t::stepper_type::magnetic_field_type field_view,
    const typename propagator_t::actor_chain_type::state_tuple
        *device_actor_state_ptr,
    vecmem::data::vector_view<
        free_track_parameters<typename propagator_t::algebra_type>>
        tracks_view) {

    using detector_device_t =
        detector<typename propagator_t::detector_type::metadata,
                 device_container_types>;
    using algebra_t = typename detector_device_t::algebra_type;
    using actor_chain_t = typename propagator_t::actor_chain_type;
    using propagator_device_t =
        propagator<typename propagator_t::stepper_type,
                   caching_navigator<detector_device_t>, actor_chain_t>;

    const detector_device_t det(det_view);
    const vecmem::device_vector<free_track_parameters<algebra_t>> tracks(
        tracks_view);

    int gid = threadIdx.x + blockIdx.x * blockDim.x;
    if (gid >= static_cast<int>(tracks.size())) {
        return;
    }

    // Create propagator
    propagator_device_t p{cfg};

    // Create the actor states on a fresh copy
    typename actor_chain_t::state_tuple actor_states = *device_actor_state_ptr;
    auto actor_state_refs = actor_chain_t::setup_actor_states(actor_states);

    // Create the propagator state

    // The track gets copied into the stepper state, so that the
    // original track sample vector remains unchanged
    typename propagator_device_t::state p_state(tracks.at(gid), field_view,
                                                det);

    // Particle hypothesis
    auto &ptc = p_state.stepping().particle_hypothesis();
    p_state.set_particle(update_particle_hypothesis(ptc, tracks.at(gid)));

    // Run propagation
    if constexpr (kOPT == detray::benchmarks::propagation_opt::e_unsync) {
        p.propagate(p_state, actor_state_refs);
    } else if constexpr (kOPT == detray::benchmarks::propagation_opt::e_sync) {
        /* Do nothing for now */
    }
}

template <typename propagator_t>
typename propagator_t::actor_chain_type::state_tuple *setup_actor_states(
    typename propagator_t::actor_chain_type::state_tuple *input_actor_states) {

    // Copy the actor state blueprint to the device
    using actor_state_t = typename propagator_t::actor_chain_type::state_tuple;
    actor_state_t *device_actor_state_ptr{nullptr};

    [[maybe_unused]] hipError_t success =
        hipMalloc((void **)&device_actor_state_ptr, sizeof(actor_state_t));
    assert(success == hipSuccess);

    success = hipMemcpy(device_actor_state_ptr, input_actor_states,
                        sizeof(actor_state_t), hipMemcpyHostToDevice);
    assert(success == hipSuccess);

    return device_actor_state_ptr;
}

template <typename propagator_t>
void release_actor_states(typename propagator_t::actor_chain_type::state_tuple
                              *device_actor_state_ptr) {
    [[maybe_unused]] hipError_t success = hipFree(device_actor_state_ptr);
    assert(success == hipSuccess);
}

template <typename propagator_t, detray::benchmarks::propagation_opt kOPT>
void run_propagation_kernel(
    const propagation::config &cfg,
    typename propagator_t::detector_type::view_type det_view,
    typename propagator_t::stepper_type::magnetic_field_type field_view,
    typename propagator_t::actor_chain_type::state_tuple
        *device_actor_state_ptr,
    vecmem::data::vector_view<
        free_track_parameters<typename propagator_t::algebra_type>>
        tracks_view,
    const int n_samples) {

    constexpr int thread_dim = 256;
    int block_dim = (n_samples + thread_dim - 1) / thread_dim;

    // run the test kernel
    hipLaunchKernelGGL((propagator_benchmark_kernel<propagator_t, kOPT>),
                       dim3(block_dim), dim3(thread_dim), 0, 0, cfg, det_view,
                       field_view, device_actor_state_ptr, tracks_view);

    // hip error check
    DETRAY_HIP_ERROR_CHECK(hipGetLastError());
    DETRAY_HIP_ERROR_CHECK(hipDeviceSynchronize());
}

/// Macro declaring the template instantiations for the different detector types
#define DECLARE_PROPAGATION_BENCHMARK(METADATA, CHAIN, FIELD, OPT)            \
                                                                              \
    template void                                                             \
    run_propagation_kernel<hip_propagator_type<METADATA, FIELD, CHAIN>, OPT>( \
        const propagation::config &, detector<METADATA>::view_type,           \
        covfie::field_view<FIELD>,                                            \
        hip_propagator_type<METADATA, FIELD,                                  \
                            CHAIN>::actor_chain_type::state_tuple *,          \
        vecmem::data::vector_view<                                            \
            free_track_parameters<detector<METADATA>::algebra_type>>,         \
        const int);                                                           \
                                                                              \
    template hip_propagator_type<METADATA, FIELD,                             \
                                 CHAIN>::actor_chain_type::state_tuple *      \
    setup_actor_states<hip_propagator_type<METADATA, FIELD, CHAIN>>(          \
        hip_propagator_type<METADATA, FIELD,                                  \
                            CHAIN>::actor_chain_type::state_tuple *);         \
                                                                              \
    template void                                                             \
    release_actor_states<hip_propagator_type<METADATA, FIELD, CHAIN>>(        \
        hip_propagator_type<METADATA, FIELD,                                  \
                            CHAIN>::actor_chain_type::state_tuple *);

DECLARE_PROPAGATION_BENCHMARK(benchmarks::default_metadata, empty_chain,
                              const_field_t, propagation_opt::e_unsync)
DECLARE_PROPAGATION_BENCHMARK(benchmarks::default_metadata, default_chain,
                              const_field_t, propagation_opt::e_unsync)

DECLARE_PROPAGATION_BENCHMARK(benchmarks::toy_metadata, empty_chain,
                              const_field_t, propagation_opt::e_unsync)
DECLARE_PROPAGATION_BENCHMARK(benchmarks::toy_metadata, default_chain,
                              const_field_t, propagation_opt::e_unsync)

DECLARE_PROPAGATION_BENCHMARK(benchmarks::wire_chamber_metadata, empty_chain,
                              const_field_t, propagation_opt::e_unsync)
DECLARE_PROPAGATION_BENCHMARK(benchmarks::wire_chamber_metadata, default_chain,
                              const_field_t, propagation_opt::e_unsync)

}  // namespace detray::benchmarks
