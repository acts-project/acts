/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/edm/track_container.hpp"
#include "traccc/finding/actors/expected_layer_pattern_collector.hpp"
#include "traccc/finding/details/combinatorial_kalman_filter_types.hpp"
#include "traccc/finding/finding_config.hpp"
#include "traccc/utils/logging.hpp"
#include "traccc/utils/particle.hpp"

// Detray include(s).
#include <detray/geometry/shapes/line.hpp>
#include <detray/geometry/tracking_surface.hpp>
#include <detray/material/material_rod.hpp>
#include <detray/navigation/direct_navigator.hpp>
#include <detray/navigation/external_surface.hpp>
#include <detray/propagator/actor_chain.hpp>
#include <detray/propagator/actors/parameter_updater.hpp>
#include <detray/propagator/constrained_step.hpp>
#include <detray/propagator/perigee_extrapolator.hpp>
#include <detray/propagator/propagator.hpp>
#include <detray/utils/ranges/single.hpp>

// System include(s).
#include <limits>
#include <vector>

namespace traccc::details {

/// Build an expected-layer pattern from a sequence of detector barcodes.
template <typename geo_id_range_t>
TRACCC_HOST_DEVICE inline expected_layer_pattern_type
collect_expected_layer_pattern(const barcode_range_t& geo_ids,
                               const expected_layer_table_mapper& mapper) {
  expected_layer_pattern_type pattern{};

  for (const auto& barcode : barcodes) {
    const auto mapping = mapper(barcode);
    if (!mapping.valid || mapping.pattern_index >= pattern.size() ||
        mapping.layer_index >=
            static_cast<unsigned int>(sizeof(unsigned int) * 8u)) {
      continue;
    }

    pattern[mapping.pattern_index] |= (1u << mapping.layer_index);
  }

  return pattern;
}

/// Build an expected-layer pattern for a track by extrapolating from the
/// perigee to the first sensitive surface.
template <typename detector_t, typename bfield_t, typename track_t>
TRACCC_HOST_DEVICE inline expected_layer_pattern_type
collect_expected_layer_pattern_from_perigee(
    const detector_t& det, const bfield_t& field, const track_t& track,
    const finding_config& config, const expected_layer_mapping_entry* map,
    std::size_t map_size) {
  expected_layer_pattern_type pattern{};

  if (map == nullptr || map_size == 0u) {
    TRACCC_INFO_HOST_DEVICE(
        "Expected-layer pattern collection skipped: map=%p size=%llu",
        static_cast<const void*>(map),
        static_cast<unsigned long long>(map_size));
    return pattern;
  }

  using stepper_t = ckf_stepper_t<bfield_t>;

  // Start from fitted/found bound parameters on the current surface.
  const auto bound_params = track.params();
  const detray::tracking_surface sf{det, bound_params.surface_link()};
  const auto free_seed = sf.bound_to_free_vector(
      typename detector_t::geometry_context{}, bound_params);

  // Extrapolate to perigee first, then propagate outward for layer tagging.
  using nav_link_t = typename detector_t::surface_type::navigation_link;
  using perigee_surface_t = detray::external_surface<
      typename detector_t::algebra_type, detray::line_circular,
      detray::material_rod<traccc::scalar>, nav_link_t>;
  const perigee_surface_t perigee_surface{
      detray::dtransform3D<traccc::default_algebra>{},
      typename perigee_surface_t::mask_type{
          0u, std::numeric_limits<traccc::scalar>::max(),
          -std::numeric_limits<traccc::scalar>::max()},
      typename perigee_surface_t::material_type{
          detray::vacuum<traccc::scalar>{},
          std::numeric_limits<traccc::scalar>::max()},
      0u, detray::surface_id::e_passive};
  detray::ranges::single_view<perigee_surface_t> perigee_view{perigee_surface};

  const auto perigee_free =
      [&]() -> typename stepper_t::free_track_parameters_type {
#if defined(__CUDA_ARCH__)
    using perigee_nav_t = detray::direct_navigator<
        detector_t, detray::ranges::single_view<perigee_surface_t>>;
    using perigee_extrapolation_propagator_t =
        detray::propagator<stepper_t, perigee_nav_t, detray::actor_chain<>>;

    perigee_extrapolation_propagator_t perigee_to_perigee{config.propagation};
    typename perigee_extrapolation_propagator_t::state perigee_state(
        free_seed, field, det, perigee_view, config.propagation.context);
    perigee_state.set_particle(traccc::detail::correct_particle_hypothesis(
        config.ptc_hypothesis, free_seed));
    perigee_state.stepping()
        .template set_constraint<detray::step::constraint::e_accuracy>(
            config.propagation.stepping.step_constraint);
    perigee_state.navigation().set_direction(
        detray::navigation::direction::e_backward);

    perigee_to_perigee.propagate(perigee_state);
    const bool perigee_finished = perigee_to_perigee.finished(perigee_state);
    const bool perigee_invalid = perigee_state.stepping()().is_invalid();
    if (!perigee_finished || perigee_invalid) {
      const auto nav_status =
          static_cast<unsigned int>(perigee_state.navigation().status());
      const auto nav_target =
          perigee_state.navigation().target().surface().identifier().value();
      const auto nav_current =
          perigee_state.navigation().current().surface().identifier().value();
      TRACCC_WARNING_HOST_DEVICE(
          "Perigee extrapolation (CUDA direct navigator) failed: finished=%u "
          "invalid=%u nav_alive=%u nav_status=%u path=%f target_id=%llu "
          "current_id=%llu",
          perigee_finished ? 1u : 0u, perigee_invalid ? 1u : 0u,
          perigee_state.navigation().is_alive() ? 1u : 0u, nav_status,
          static_cast<double>(perigee_state.stepping().path_length()),
          static_cast<unsigned long long>(nav_target),
          static_cast<unsigned long long>(nav_current));
      return typename stepper_t::free_track_parameters_type{};
    }

    return perigee_state.stepping()();
#else
    using perigee_extrapolator_t =
        detray::perigee_extrapolator<detector_t, stepper_t>;
    perigee_extrapolator_t perigee_extrapolator{config.propagation};

    typename perigee_extrapolator_t::state perigee_state(
        free_seed, field, det, perigee_view, config.propagation.context);
    perigee_state.set_particle(traccc::detail::correct_particle_hypothesis(
        config.ptc_hypothesis, free_seed));
    perigee_state.stepping()
        .template set_constraint<detray::step::constraint::e_accuracy>(
            config.propagation.stepping.step_constraint);

    (void)perigee_extrapolator.extrapolate(perigee_state);
    if (!perigee_extrapolator.finished(perigee_state) ||
        perigee_state.stepping()().is_invalid()) {
      TRACCC_WARNING_HOST_DEVICE(
          "Perigee extrapolation (host extrapolator path) failed: finished=%u "
          "invalid=%u",
          perigee_extrapolator.finished(perigee_state) ? 1u : 0u,
          perigee_state.stepping()().is_invalid() ? 1u : 0u);
      return typename stepper_t::free_track_parameters_type{};
    }

    return perigee_state.stepping()();
#endif
  }();

  if (perigee_free.is_invalid()) {
    TRACCC_WARNING_HOST_DEVICE(
        "Expected-layer pattern collection aborted: invalid perigee-free "
        "parameters");
    return pattern;
  }

  // Run a CKF-like actor chain with expected-layer collector enabled.
  using perigee_actor_chain_t = detray::actor_chain<
      detray::actor::pathlimit_aborter<traccc::scalar>,
      detray::actor::parameter_updater<traccc::default_algebra,
                                       ckf_interactor_t>,
      detray::actor::momentum_aborter<traccc::scalar>,
      expected_layer_pattern_collector<expected_layer_table_mapper>>;
  using perigee_propagator_t = detray::propagator<
      stepper_t, detray::caching_navigator<std::add_const_t<detector_t>>,
      perigee_actor_chain_t>;

  detray::propagation::config prop_cfg{config.propagation};
  perigee_propagator_t perigee_to_first(prop_cfg);
  typename perigee_propagator_t::state perigee_propagation(
      perigee_free, field, det, config.propagation.context);
  perigee_propagation.set_particle(traccc::detail::correct_particle_hypothesis(
      config.ptc_hypothesis, perigee_free));
  perigee_propagation.stepping()
      .template set_constraint<detray::step::constraint::e_accuracy>(
          config.propagation.stepping.step_constraint);

  typename detray::actor::pathlimit_aborter<traccc::scalar>::state
      aborter_state{};
  detray::actor::parameter_updater_state<typename detector_t::algebra_type>
      updater_state{prop_cfg, bound_params};
  traccc::details::ckf_interactor_t::state interactor_state{};
  typename detray::actor::momentum_aborter<traccc::scalar>::state
      momentum_aborter_state{};
  typename expected_layer_pattern_collector<expected_layer_table_mapper>::state
      expected_layer_collector_state{};

  updater_state.notify_on_initial(true);
  momentum_aborter_state.min_pT(static_cast<traccc::scalar>(config.min_pT));
  momentum_aborter_state.min_p(static_cast<traccc::scalar>(config.min_p));
  // Write collected layer bits into caller-provided output pattern.
  expected_layer_collector_state.pattern = &pattern;
  expected_layer_collector_state.mapper.entries = map;
  expected_layer_collector_state.mapper.size = map_size;

  perigee_to_first.propagate(
      perigee_propagation,
      detray::tie(aborter_state, updater_state, interactor_state,
                  momentum_aborter_state, expected_layer_collector_state));

  TRACCC_DEBUG_HOST_DEVICE(
      "Expected-layer pattern collection: %u updates, %u skipped, pattern={%u, "
      "%u, %u, %u}",
      expected_layer_collector_state.n_updates,
      expected_layer_collector_state.n_skipped, pattern[0], pattern[1],
      pattern[2], pattern[3]);

  return pattern;
}

/// Build expected-layer patterns for all tracks in a track container along with
/// perigee extrapolation.
template <typename detector_t, typename bfield_t>
TRACCC_HOST inline std::vector<expected_layer_pattern_type>
expected_layer_patterns_from_perigee(
    const detector_t& det, const bfield_t& field,
    const typename traccc::edm::track_container<
        typename detector_t::algebra_type>::const_device& tracks,
    const finding_config& config, const expected_layer_mapping_entry* map,
    std::size_t map_size) {
  std::vector<expected_layer_pattern_type> patterns;
  patterns.reserve(tracks.tracks.size());

  // One output pattern per input track.
  for (std::size_t track_index = 0u; track_index < tracks.tracks.size();
       ++track_index) {
    patterns.push_back(collect_expected_layer_pattern_from_perigee(
        det, field, tracks.tracks.at(track_index), config, map, map_size));
  }

  return patterns;
}

}  // namespace traccc::details
