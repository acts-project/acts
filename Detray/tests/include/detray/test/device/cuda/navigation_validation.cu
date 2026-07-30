// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "detray/definitions/detail/cuda_definitions.hpp"
#include "navigation_validation.hpp"

namespace detray::cuda {

template <typename bfield_t, typename detector_t>
__global__ void navigation_validation_kernel(
    typename detector_t::view_type det_data, const propagation::config cfg,
    pdg_particle<typename detector_t::scalar_type> ptc_hypo,
    bfield_t field_data,
    vecmem::data::jagged_vector_view<const intersection_record<detector_t>>
        truth_intersection_traces_view,
    vecmem::data::jagged_vector_view<intersection_record<detector_t>>
        recorded_intersections_view,
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

  static_assert(std::is_same_v<typename detector_t::view_type,
                               typename detector_device_t::view_type>,
                "Host and device detector view types do not match");

  using hom_bfield_view_t = typename bfield::const_field_t<scalar_t>::view_t;
  using rk_stepper_t = rk_stepper<hom_bfield_view_t, algebra_t>;
  using line_stepper_t = line_stepper<algebra_t>;
  // Use RK-stepper when a non-empty b-field was passed
  static constexpr auto is_no_bfield{
      std::is_same_v<bfield_t, navigation_validator::empty_bfield>};
  using stepper_t =
      std::conditional_t<is_no_bfield, line_stepper_t, rk_stepper_t>;

  // Inspector that records all encountered surfaces
  using intersection_t =
      typename intersection_record<detector_t>::intersection_type;
  using object_tracer_t =
      navigation::object_tracer<detector_t, vecmem::device_vector,
                                navigation::status::e_on_object,
                                navigation::status::e_on_portal>;
  // Navigation with inspection
  using navigator_t =
      caching_navigator<detector_device_t, navigation::default_cache_size,
                        object_tracer_t, intersection_t>;

  // Propagator with pathlimit aborter
  using material_tracer_t =
      material_validator::material_tracer<scalar_t, vecmem::device_vector>;
  using pathlimit_aborter_t = actor::pathlimit_aborter<scalar_t>;
  using actor_chain_t = actor_chain<pathlimit_aborter_t, material_tracer_t>;
  using propagator_t = propagator<stepper_t, navigator_t, actor_chain_t>;

  detector_device_t det(det_data);

  vecmem::jagged_device_vector<const intersection_record<detector_t>>
      truth_intersection_traces(truth_intersection_traces_view);
  vecmem::jagged_device_vector<intersection_record<detector_t>>
      recorded_intersections(recorded_intersections_view);
  vecmem::device_vector<typename material_tracer_t::track_material_type>
      track_mat_vec(track_mat_view);
  vecmem::jagged_device_vector<typename material_tracer_t::material_record_type>
      mat_steps(mat_steps_view);

  // Check the memory setup
  assert(truth_intersection_traces.size() ==
         recorded_intersections_view.size());

  int trk_id = threadIdx.x + blockIdx.x * blockDim.x;
  if (trk_id >= truth_intersection_traces.size()) {
    return;
  }

  propagator_t p{cfg};

  // Create the actor states
  typename pathlimit_aborter_t::state aborter_state{cfg.stepping.path_limit};
  typename material_tracer_t::state mat_tracer_state{mat_steps.at(trk_id)};
  auto actor_states = ::detray::tie(aborter_state, mat_tracer_state);

  // Get the initial track parameters
  const auto &track = truth_intersection_traces[trk_id].front().track_param();

  // Save the initial intersection, since it is not recorded by the
  // object tracer
  assert(recorded_intersections.at(trk_id).empty());
  recorded_intersections.at(trk_id).push_back(
      {track.pos(), track.dir(),
       truth_intersection_traces[trk_id].front().intersection});
  // Did the insertion of an element work?
  assert(recorded_intersections.at(trk_id).size() == 1);

  // Run propagation
  if constexpr (is_no_bfield) {
    typename propagator_t::state propagation(
        track, det,
        typename navigator_t::state::view_type{
            recorded_intersections_view.ptr()[trk_id]});
    propagation.set_particle(update_particle_hypothesis(ptc_hypo, track));

    p.propagate(propagation, actor_states);
  } else {
    typename propagator_t::stepper_type::magnetic_field_type bfield_view(
        field_data);
    typename propagator_t::state propagation(
        track, bfield_view, det,
        typename navigator_t::state::view_type{
            recorded_intersections_view.ptr()[trk_id]});
    propagation.set_particle(update_particle_hypothesis(ptc_hypo, track));

    p.propagate(propagation, actor_states);
  }

  // Record the accumulated material
  assert(truth_intersection_traces.size() == track_mat_vec.size());
  track_mat_vec.at(trk_id) = mat_tracer_state.get_track_material();
}

/// Launch the device kernel
template <typename bfield_t, typename detector_t>
void navigation_validation_device(
    typename detector_t::view_type det_view, const propagation::config &cfg,
    pdg_particle<typename detector_t::scalar_type> ptc_hypo,
    bfield_t field_data,
    vecmem::data::jagged_vector_view<const intersection_record<detector_t>>
        &truth_intersection_traces_view,
    vecmem::data::jagged_vector_view<intersection_record<detector_t>>
        &recorded_intersections_view,
    vecmem::data::vector_view<
        material_validator::track_material<typename detector_t::scalar_type>>
        &track_mat_view,
    vecmem::data::jagged_vector_view<
        material_record<typename detector_t::scalar_type>> &mat_steps_view) {
  constexpr int thread_dim = 2 * WARP_SIZE;
  int block_dim = truth_intersection_traces_view.size() / thread_dim + 1;

  // run the test kernel
  navigation_validation_kernel<bfield_t, detector_t><<<block_dim, thread_dim>>>(
      det_view, cfg, ptc_hypo, field_data, truth_intersection_traces_view,
      recorded_intersections_view, track_mat_view, mat_steps_view);

  // cuda error check
  DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
  DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

/// Macro declaring the template instantiations for the different detector types
#define DECLARE_NAVIGATION_VALIDATION(METADATA)                            \
                                                                           \
  template void navigation_validation_device<                              \
      covfie::field_view<                                                  \
          bfield::const_bknd_t<dscalar<typename METADATA::algebra_type>>>, \
      detector<METADATA>>(                                                 \
      typename detector<METADATA>::view_type, const propagation::config &, \
      pdg_particle<typename detector<METADATA>::scalar_type>,              \
      covfie::field_view<                                                  \
          bfield::const_bknd_t<dscalar<typename METADATA::algebra_type>>>, \
      vecmem::data::jagged_vector_view<                                    \
          const detray::intersection_record<detector<METADATA>>> &,        \
      vecmem::data::jagged_vector_view<                                    \
          intersection_record<detector<METADATA>>> &,                      \
      vecmem::data::vector_view<material_validator::track_material<        \
          typename detector<METADATA>::scalar_type>> &,                    \
      vecmem::data::jagged_vector_view<                                    \
          material_record<typename detector<METADATA>::scalar_type>> &);   \
                                                                           \
  template void navigation_validation_device<                              \
      detray::navigation_validator::empty_bfield, detector<METADATA>>(     \
      typename detector<METADATA>::view_type, const propagation::config &, \
      pdg_particle<typename detector<METADATA>::scalar_type>,              \
      detray::navigation_validator::empty_bfield,                          \
      vecmem::data::jagged_vector_view<                                    \
          const detray::intersection_record<detector<METADATA>>> &,        \
      vecmem::data::jagged_vector_view<                                    \
          intersection_record<detector<METADATA>>> &,                      \
      vecmem::data::vector_view<material_validator::track_material<        \
          typename detector<METADATA>::scalar_type>> &,                    \
      vecmem::data::jagged_vector_view<                                    \
          material_record<typename detector<METADATA>::scalar_type>> &);

DECLARE_NAVIGATION_VALIDATION(test::default_metadata)
DECLARE_NAVIGATION_VALIDATION(test::toy_metadata)
DECLARE_NAVIGATION_VALIDATION(test::default_telescope_metadata)
DECLARE_NAVIGATION_VALIDATION(test::wire_chamber_metadata)

}  // namespace detray::cuda
