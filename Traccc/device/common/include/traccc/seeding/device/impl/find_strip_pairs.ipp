/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <vecmem/memory/device_atomic_ref.hpp>

#include "traccc/seeding/detail/strip_geometry.hpp"
#include "traccc/seeding/detail/strip_pairing.hpp"

namespace traccc::device {

struct strip_pair_write_visitor {
  unsigned int reference_index;
  unsigned int& standard_position;
  unsigned int& overlap_position;
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
    const auto position =
        vecmem::device_atomic_ref<unsigned int>(counter).fetch_add(1u);
    if (position < output.size()) {
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
    const point3& beam_spot, unsigned int& opposite_position,
    unsigned int& overlap_position,
    strip_pair_collection_types::view opposite_pairs_view,
    strip_pair_collection_types::view overlap_pairs_view) {
  const edm::measurement_collection::const_device measurements(
      measurements_view);
  const strip_measurement_surface_info_collection_types::const_device
      surface_infos(surface_infos_view);
  const strip_pairing_rule_collection_types::const_device rules(rules_view);

  vecmem::device_vector<strip_pair> opposite_pairs(opposite_pairs_view);
  vecmem::device_vector<strip_pair> overlap_pairs(overlap_pairs_view);
  strip_pair_write_visitor visitor{static_cast<unsigned int>(globalIndex),
                                   opposite_position, overlap_position,
                                   opposite_pairs, overlap_pairs};

  details::visit_strip_pairs(static_cast<unsigned int>(globalIndex),
                             measurements, surface_infos, rules, beam_spot,
                             visitor);
}

}  // namespace traccc::device
