/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/memory/device_atomic_ref.hpp>

// System include(s).
#include <bit>

namespace traccc::device {

TRACCC_HOST_DEVICE inline unsigned int encode_seed_weight(float weight) {
  const unsigned int bits = std::bit_cast<unsigned int>(weight);
  return (bits & 0x80000000u) ? ~bits : (bits | 0x80000000u);
}

TRACCC_HOST_DEVICE inline void mark_seed_weights(
    const global_index_t globalIndex,
    const edm::seed_collection::const_view& seeds_view,
    vecmem::data::vector_view<unsigned int> best_weight_view) {
  const edm::seed_collection::const_device seeds(seeds_view);
  if (globalIndex >= seeds.size()) {
    return;
  }
  vecmem::device_vector<unsigned int> best_weight(best_weight_view);

  const auto seed = seeds.at(globalIndex);
  const unsigned int weight = encode_seed_weight(seed.quality());
  for (const unsigned int sp :
       {seed.bottom_index(), seed.middle_index(), seed.top_index()}) {
    vecmem::device_atomic_ref<unsigned int> best(best_weight.at(sp));
    best.fetch_max(weight);
  }
}

TRACCC_HOST_DEVICE inline void suppress_seeds(
    const global_index_t globalIndex,
    const edm::seed_collection::const_view& seeds_view,
    const vecmem::data::vector_view<const unsigned int>& best_weight_view,
    edm::seed_collection::view out_view) {
  const edm::seed_collection::const_device seeds(seeds_view);
  if (globalIndex >= seeds.size()) {
    return;
  }
  const vecmem::device_vector<const unsigned int> best_weight(best_weight_view);
  edm::seed_collection::device out(out_view);

  const auto seed = seeds.at(globalIndex);
  const unsigned int weight = encode_seed_weight(seed.quality());
  if (weight < best_weight.at(seed.bottom_index()) &&
      weight < best_weight.at(seed.middle_index()) &&
      weight < best_weight.at(seed.top_index())) {
    return;
  }
  out.push_back({seed.bottom_index(), seed.middle_index(), seed.top_index(),
                 seed.quality()});
}

}  // namespace traccc::device
