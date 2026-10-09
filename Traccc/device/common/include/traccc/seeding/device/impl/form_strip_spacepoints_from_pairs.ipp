// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>

// Project include(s).
#include "traccc/seeding/detail/spacepoint_formation.hpp"
#include "traccc/seeding/detail/strip_geometry.hpp"

namespace traccc::device {

TRACCC_HOST_DEVICE inline void form_strip_spacepoints_from_pairs(
    const global_index_t globalIndex,
    const edm::measurement_collection::const_view& measurements_view,
    const strip_pair_collection_types::const_view& pairs_view,
    const strip_measurement_surface_info_collection_types::const_view&
        surface_infos_view,
    const point3& beam_spot,
    vecmem::data::vector_view<unsigned int> accepted_view,
    edm::spacepoint_collection::view spacepoints_view) {
  const edm::measurement_collection::const_device measurements(
      measurements_view);
  const strip_pair_collection_types::const_device pairs(pairs_view);
  const strip_measurement_surface_info_collection_types::const_device
      surface_infos(surface_infos_view);
  if (globalIndex >= pairs.size()) {
    return;
  }

  vecmem::device_vector<unsigned int> accepted_flags(accepted_view);
  accepted_flags.at(globalIndex) = 0u;
  edm::spacepoint_collection::device spacepoints(spacepoints_view);
  const strip_pair pair = pairs.at(globalIndex);
  // A count pass may very rarely overestimate the number of pairs found by the
  // write pass because the two kernels can evaluate floating-point cut
  // boundaries differently. Pair buffers are pre-filled with sentinel indices;
  // reject those, as well as any other invalid index, before dereferencing the
  // measurement collection.
  if ((pair.measurement_index_1 >= measurements.size()) ||
      (pair.measurement_index_2 >= measurements.size())) {
    return;
  }
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

  const edm::spacepoint_collection::device::size_type i = globalIndex;
  edm::spacepoint_collection::device::proxy_type sp = spacepoints.at(i);
  sp.x() = static_cast<float>(position[0]);
  sp.y() = static_cast<float>(position[1]);
  sp.z() = static_cast<float>(position[2]);
  sp.radius_variance() = static_cast<float>(first_info.spacepoint_variance_r);
  sp.z_variance() = static_cast<float>(first_info.spacepoint_variance_z);
  sp.measurement_index_1() = pair.measurement_index_1;
  sp.measurement_index_2() = pair.measurement_index_2;
  accepted_flags.at(globalIndex) = 1u;
}

/// Compact accepted spacepoints in pair order using inclusive flag offsets.
TRACCC_HOST_DEVICE inline void gather_strip_spacepoints(
    const global_index_t globalIndex,
    const edm::spacepoint_collection::const_view& candidates_view,
    vecmem::data::vector_view<const unsigned int> offsets_view,
    edm::spacepoint_collection::view output_view) {
  const edm::spacepoint_collection::const_device candidates(candidates_view);
  if (globalIndex >= candidates.size()) {
    return;
  }
  const vecmem::device_vector<const unsigned int> offsets(offsets_view);
  const unsigned int end = offsets.at(globalIndex);
  const unsigned int start =
      globalIndex == 0u ? 0u : offsets.at(globalIndex - 1u);
  if (end == start) {
    return;
  }
  edm::spacepoint_collection::device output(output_view);
  const auto source = candidates.at(globalIndex);
  auto destination = output.at(end - 1u);
  destination.x() = source.x();
  destination.y() = source.y();
  destination.z() = source.z();
  destination.radius_variance() = source.radius_variance();
  destination.z_variance() = source.z_variance();
  destination.measurement_index_1() = source.measurement_index_1();
  destination.measurement_index_2() = source.measurement_index_2();
}

}  // namespace traccc::device
