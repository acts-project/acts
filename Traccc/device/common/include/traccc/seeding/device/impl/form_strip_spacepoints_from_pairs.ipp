/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/seeding/detail/spacepoint_formation.hpp"
#include "traccc/seeding/detail/strip_geometry.hpp"

namespace traccc::device {

template <detray::concepts::detector detector_t>
TRACCC_HOST_DEVICE inline void form_strip_spacepoints_from_pairs(
    const global_index_t globalIndex,
    const detray::detector_view_t<detector_t> det_view,
    const edm::measurement_collection::const_view& measurements_view,
    const strip_pair_collection_types::const_view& pairs_view,
    const strip_measurement_surface_info_collection_types::const_view&
        surface_infos_view,
    const point3& beam_spot,
    edm::spacepoint_collection::view spacepoints_view) {
  const edm::measurement_collection::const_device measurements(
      measurements_view);
  const strip_pair_collection_types::const_device pairs(pairs_view);
  const strip_measurement_surface_info_collection_types::const_device
      surface_infos(surface_infos_view);
  if (globalIndex >= pairs.size()) {
    return;
  }

  (void)det_view;
  edm::spacepoint_collection::device spacepoints(spacepoints_view);
  const strip_pair pair = pairs.at(globalIndex);
  const edm::measurement first_meas = measurements.at(pair.measurement_index_1);
  const edm::measurement second_meas =
      measurements.at(pair.measurement_index_2);

  strip_measurement_surface_info first_info{};
  strip_measurement_surface_info second_info{};
  if (!details::find_strip_surface_info(
          surface_infos, first_meas.surface_link().value(), first_info) ||
      !details::find_strip_surface_info(
          surface_infos, second_meas.surface_link().value(), second_info)) {
    return;
  }
  const details::strip_material first_material =
      details::make_strip_material(first_meas, first_info, beam_spot);
  const details::strip_material second_material =
      details::make_strip_material(second_meas, second_info, beam_spot);
  if ((first_material.valid == 0u) || (second_material.valid == 0u)) {
    return;
  }
  const scalar strip_length_gap_tolerance = pair.strip_length_gap_tolerance;
  point3 position{};
  const bool accepted = traccc::details::make_strip_spacepoint(
      position, first_material.center, first_material.direction,
      second_material.direction, first_material.trajectory,
      second_material.trajectory, first_material.normal, second_material.normal,
      first_material.half_length, second_material.half_length,
      strip_length_gap_tolerance, pair.strip_length_tolerance);
  if (!accepted) {
    return;
  }

  const edm::spacepoint_collection::device::size_type i =
      spacepoints.push_back_default();
  edm::spacepoint_collection::device::proxy_type sp = spacepoints.at(i);
  sp.x() = static_cast<float>(position[0]);
  sp.y() = static_cast<float>(position[1]);
  sp.z() = static_cast<float>(position[2]);
  sp.radius_variance() = static_cast<float>(first_info.spacepoint_variance_r);
  sp.z_variance() = static_cast<float>(first_info.spacepoint_variance_z);
  sp.measurement_index_1() = pair.measurement_index_1;
  sp.measurement_index_2() = pair.measurement_index_2;
}

}  // namespace traccc::device
