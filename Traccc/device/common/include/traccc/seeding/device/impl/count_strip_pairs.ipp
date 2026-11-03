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

struct strip_pair_count_visitor {
  unsigned int standard{0u}, overlap{0u};
  TRACCC_HOST_DEVICE void operator()(unsigned int,
                                     const strip_pairing_rule& rule) {
    if (rule.category == strip_pair_category::standard) {
      ++standard;
    } else {
      ++overlap;
    }
  }
};
TRACCC_HOST_DEVICE inline void count_strip_pairs(
    const global_index_t globalIndex,
    const edm::measurement_collection::const_view& measurements_view,
    const strip_measurement_surface_info_collection_types::const_view&
        surface_infos_view,
    const strip_pairing_rule_collection_types::const_view& rules_view,
    const point3& beam_spot,
    vecmem::data::vector_view<unsigned int> standard_counts_view,
    vecmem::data::vector_view<unsigned int> overlap_counts_view) {
  const edm::measurement_collection::const_device measurements(
      measurements_view);
  const strip_measurement_surface_info_collection_types::const_device
      surface_infos(surface_infos_view);
  const strip_pairing_rule_collection_types::const_device rules(rules_view);

  strip_pair_count_visitor visitor{};

  details::visit_strip_pairs(static_cast<unsigned int>(globalIndex),
                             measurements, surface_infos, rules, beam_spot,
                             visitor);

  if (globalIndex < measurements.size()) {
    vecmem::device_vector<unsigned int> standard_counts(standard_counts_view);
    vecmem::device_vector<unsigned int> overlap_counts(overlap_counts_view);
    standard_counts.at(globalIndex) = visitor.standard;
    overlap_counts.at(globalIndex) = visitor.overlap;
  }
}

}  // namespace traccc::device
