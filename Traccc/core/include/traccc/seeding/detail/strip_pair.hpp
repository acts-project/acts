/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "traccc/seeding/strip_spacepoint_formation_data.hpp"

namespace traccc {

/// Event-local pair of indices into the sorted measurement collection.
/// Carry the selected rule's tolerances through to the formation pass.
struct strip_pair {
  unsigned int measurement_index_1;
  unsigned int measurement_index_2;
  scalar strip_length_gap_tolerance;
  scalar strip_length_tolerance;
};

using strip_pair_collection_types = collection_types<strip_pair>;

}  // namespace traccc
