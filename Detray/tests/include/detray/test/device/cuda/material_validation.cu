// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "detray/definitions/detail/cuda_definitions.hpp"
#include "detray/propagator/actors.hpp"
#include "detray/propagator/line_stepper.hpp"
#include "material_validation.hpp"

namespace detray::cuda {

template <typename detector_t>
__global__ void material_validation_kernel(
    typename detector_t::view_type det_data, const propagation::config cfg,
    vecmem::data::vector_view<
        free_track_parameters<typename detector_t::algebra_type>>
        tracks_view,
    vecmem::data::vector_view<
        material_validator::track_material<typename detector_t::scalar_type>>
        track_mat_view,
    vecmem::data::jagged_vector_view<
        material_record<typename detector_t::scalar_type>>
        mat_steps_view) {
  using detector_device_t =
      detector<typename detector_t::metadata, device_container_types>;
  using algebra_t = typename detector_device_t::algebra_type;
  using scalar_t = dscalar<algebra_t>;

  using stepper_t = line_stepper<algebra_t>;
  using navigator_t = caching_navigator<detector_device_t>;
  // Propagator with full covariance transport, pathlimit aborter and
  // material tracer
  using material_tracer_t =
      material_validator::material_tracer<scalar_t, vecmem::device_vector>;
  using pathlimit_aborter_t = actor::pathlimit_aborter<scalar_t>;
  using actor_chain_t = actor_chain<
      pathlimit_aborter_t,
      actor::parameter_updater<algebra_t,
                               actor::pointwise_material_interactor<algebra_t>,
                               material_tracer_t>>;
  using propagator_t = propagator<stepper_t, navigator_t, actor_chain_t>;

  detector_device_t det(det_data);

  vecmem::device_vector<free_track_parameters<algebra_t>> tracks(tracks_view);
  vecmem::device_vector<typename material_tracer_t::track_material_type>
      track_mat_vec(track_mat_view);
  vecmem::jagged_device_vector<typename material_tracer_t::material_record_type>
      mat_steps(mat_steps_view);

  int trk_id = threadIdx.x + blockIdx.x * blockDim.x;
  if (trk_id >= tracks.size()) {
    return;
  }

  propagator_t p{cfg};

  // Create the actor states
  typename pathlimit_aborter_t::state aborter_state{cfg.stepping.path_limit};
  actor::parameter_updater_state<algebra_t> updater_state{cfg};
  typename actor::pointwise_material_interactor<algebra_t>::state
      interactor_state{};
  typename material_tracer_t::state mat_tracer_state{mat_steps.at(trk_id)};

  auto actor_states = ::detray::tie(aborter_state, updater_state,
                                    interactor_state, mat_tracer_state);

  // Run propagation
  typename navigator_t::state::view_type nav_view{};
  typename propagator_t::state propagation(tracks[trk_id], det, nav_view);

  p.propagate(propagation, actor_states);

  // Record the accumulated material
  assert(track_mat_vec.size() == tracks.size());
  track_mat_vec.at(trk_id) = mat_tracer_state.get_track_material();
}

/// Launch the device kernel
template <typename detector_t>
void material_validation_device(
    typename detector_t::view_type det_view, const propagation::config &cfg,
    vecmem::data::vector_view<
        free_track_parameters<typename detector_t::algebra_type>> &tracks_view,
    vecmem::data::vector_view<
        material_validator::track_material<typename detector_t::scalar_type>>
        &track_mat_view,
    vecmem::data::jagged_vector_view<
        material_record<typename detector_t::scalar_type>> &mat_steps_view) {
  constexpr int thread_dim = 2 * WARP_SIZE;
  int block_dim = tracks_view.size() / thread_dim + 1;

  // run the test kernel
  material_validation_kernel<detector_t><<<block_dim, thread_dim>>>(
      det_view, cfg, tracks_view, track_mat_view, mat_steps_view);

  // cuda error check
  DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
  DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

/// Macro declaring the template instantiations for the different detector types
#define DECLARE_MATERIAL_VALIDATION(METADATA)                                  \
                                                                               \
  template void material_validation_device<detector<METADATA>>(                \
      typename detector<METADATA>::view_type, const propagation::config &,     \
      vecmem::data::vector_view<                                               \
          free_track_parameters<typename detector<METADATA>::algebra_type>> &, \
      vecmem::data::vector_view<material_validator::track_material<            \
          typename detector<METADATA>::scalar_type>> &,                        \
      vecmem::data::jagged_vector_view<                                        \
          material_record<typename detector<METADATA>::scalar_type>> &);

DECLARE_MATERIAL_VALIDATION(test::default_metadata)
DECLARE_MATERIAL_VALIDATION(test::toy_metadata)
DECLARE_MATERIAL_VALIDATION(test::default_telescope_metadata)
DECLARE_MATERIAL_VALIDATION(test::wire_chamber_metadata)

}  // namespace detray::cuda
