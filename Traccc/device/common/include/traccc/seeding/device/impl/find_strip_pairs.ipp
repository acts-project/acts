// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <vecmem/containers/device_vector.hpp>

#include "traccc/seeding/detail/strip_geometry.hpp"
#include "traccc/seeding/detail/strip_pairing.hpp"

namespace traccc::device {

struct strip_pair_write_visitor {
  unsigned int reference_index;
  unsigned int standard_position;
  unsigned int overlap_position;
  unsigned int standard_end;
  unsigned int overlap_end;
  vecmem::device_vector<strip_pair>& standard;
  vecmem::device_vector<strip_pair>& overlap;
  TRACCC_HOST_DEVICE void operator()(unsigned int candidate_index,
                                     const strip_pairing_rule& rule) {
    const strip_pair pair{reference_index, candidate_index,
                          rule.strip_length_gap_tolerance,
                          rule.strip_length_tolerance};
    auto& output =
        rule.category == strip_pair_category::standard ? standard : overlap;
    auto& counter = rule.category == strip_pair_category::standard
                        ? standard_position
                        : overlap_position;
    const auto end = rule.category == strip_pair_category::standard
                         ? standard_end
                         : overlap_end;
    // Each measurement owns only the range assigned by the count/scan pass.
    // Never overwrite the following measurement if find produces more pairs.
    if ((counter < end) && (counter < output.size())) {
      const auto position = counter++;
      output.at(position) = pair;
    }
  }
};
TRACCC_HOST_DEVICE inline void find_strip_pairs(
    const global_index_t globalIndex,
    const edm::measurement_collection::const_view& measurements_view,
    const strip_measurement_surface_info_collection_types::const_view&
        surface_infos_view,
    const strip_pairing_rule_collection_types::const_view& rules_view,
    const point3& beam_spot,
    vecmem::data::vector_view<const unsigned int> standard_offsets_view,
    vecmem::data::vector_view<const unsigned int> overlap_offsets_view,
    strip_pair_collection_types::view opposite_pairs_view,
    strip_pair_collection_types::view overlap_pairs_view) {
  const edm::measurement_collection::const_device measurements(
      measurements_view);
  const strip_measurement_surface_info_collection_types::const_device
      surface_infos(surface_infos_view);
  const strip_pairing_rule_collection_types::const_device rules(rules_view);

  vecmem::device_vector<strip_pair> opposite_pairs(opposite_pairs_view);
  vecmem::device_vector<strip_pair> overlap_pairs(overlap_pairs_view);
  if (globalIndex >= measurements.size()) {
    return;
  }
  const vecmem::device_vector<const unsigned int> standard_offsets(
      standard_offsets_view);
  const vecmem::device_vector<const unsigned int> overlap_offsets(
      overlap_offsets_view);
  strip_pair_write_visitor visitor{
      static_cast<unsigned int>(globalIndex),
      globalIndex == 0u ? 0u : standard_offsets.at(globalIndex - 1u),
      globalIndex == 0u ? 0u : overlap_offsets.at(globalIndex - 1u),
      standard_offsets.at(globalIndex),
      overlap_offsets.at(globalIndex),
      opposite_pairs,
      overlap_pairs};

  details::visit_strip_pairs(static_cast<unsigned int>(globalIndex),
                             measurements, surface_infos, rules, beam_spot,
                             visitor);
}

}  // namespace traccc::device
