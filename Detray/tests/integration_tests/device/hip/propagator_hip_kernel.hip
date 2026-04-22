/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/definitions/detail/hip_definitions.hpp"

// Detray test include(s)
#include "propagator_hip_kernel.hpp"

namespace detray {

template <typename bfield_bknd_t, typename detector_t>
__global__ void propagator_test_kernel(
    typename detector_t::view_type det_data, const propagation::config cfg,
    covfie::field_view<bfield_bknd_t> field_data,
    vecmem::data::vector_view<test_track> tracks_data,
    vecmem::data::jagged_vector_view<detail::step_data<test_algebra>>
        steps_data) {

    int gid = threadIdx.x + blockIdx.x * blockDim.x;
    using detector_device_t =
        detector<typename detector_t::metadata, device_container_types>;

    static_assert(std::is_same_v<typename detector_t::view_type,
                                 typename detector_device_t::view_type>,
                  "Host and device detector views do not match");

    detector_device_t det(det_data);
    vecmem::device_vector<test_track> tracks(tracks_data);
    vecmem::jagged_device_vector<detail::step_data<test_algebra>> steps(
        steps_data);

    if (gid >= static_cast<int>(tracks.size())) {
        return;
    }

    auto stepr = rk_stepper_t<covfie::field_view<bfield_bknd_t>>{};
    auto nav = navigator_t<detector_device_t>{};

    // Create propagator
    using propagator_device_t =
        propagator<decltype(stepr), decltype(nav), actor_chain_device_t>;

    propagator_device_t p{cfg};

    // Create actor states
    step_tracer_device_t::state tracer_state(steps.at(gid));
    pathlimit_aborter_t::state aborter_state{cfg.stepping.path_limit};
    actor::parameter_updater_state<test_algebra> updater_state{cfg};
    actor::pointwise_material_interactor<test_algebra>::state
        interactor_state{};

    // Create the actor states
    auto actor_states = ::detray::tie(tracer_state, aborter_state,
                                      updater_state, interactor_state);
    // Create the propagator state
    typename propagator_device_t::state state(tracks.at(gid), field_data, det);

    auto& ptc = state.stepping().particle_hypothesis();
    state.set_particle(update_particle_hypothesis(ptc, tracks.at(gid)));

    state.stepping().template set_constraint<step::constraint::e_accuracy>(
        cfg.stepping.step_constraint);

    // Run propagation
    p.propagate(state, actor_states);
}

/// Launch the device kernel
template <typename bfield_bknd_t, typename detector_t>
void propagator_test(
    typename detector_t::view_type det_view, const propagation::config& cfg,
    covfie::field_view<bfield_bknd_t> field_data,
    vecmem::data::vector_view<test_track>& tracks_data,
    vecmem::data::jagged_vector_view<detail::step_data<test_algebra>>&
        step_data) {

    constexpr int thread_dim = 2 * WARP_SIZE;
    int block_dim = tracks_data.size() / thread_dim + 1;

    // run the test kernel
    hipLaunchKernelGGL((propagator_test_kernel<bfield_bknd_t, detector_t>),
                       dim3(block_dim), dim3(thread_dim), 0, 0, det_view, cfg,
                       field_data, tracks_data, step_data);

    // hip error check
    DETRAY_HIP_ERROR_CHECK(hipGetLastError());
    DETRAY_HIP_ERROR_CHECK(hipDeviceSynchronize());
}

/// Explicit instantiation for a constant magnetic field
template void
propagator_test<bfield::const_bknd_t<dscalar<test_algebra>>,
                detector<toy_metadata<test_algebra>, host_container_types>>(
    detector<toy_metadata<test_algebra>, host_container_types>::view_type,
    const propagation::config&,
    covfie::field_view<bfield::const_bknd_t<dscalar<test_algebra>>>,
    vecmem::data::vector_view<test_track>&,
    vecmem::data::jagged_vector_view<detail::step_data<test_algebra>>&);

/// Explicit instantiation for an inhomogeneous magnetic field
/*
template void
propagator_test<bfield::hip::inhom_bknd_t,
                detector<toy_metadata<test_algebra>, host_container_types>>(
    detector<toy_metadata<test_algebra>, host_container_types>::view_type,
    const propagation::config&, covfie::field_view<bfield::hip::inhom_bknd_t>,
    vecmem::data::vector_view<test_track>&,
    vecmem::data::jagged_vector_view<detail::step_data<test_algebra>>&); */

}  // namespace detray
