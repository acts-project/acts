/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"

// Detray include(s).
#include <detray/geometry/identifier.hpp>
#include <detray/propagator/base_actor.hpp>

// System include(s).
#include <array>
#include <cstddef>
#include <utility>

namespace traccc {

/// Bitmask layout: [pixel barrel, pixel endcap, strip barrel, strip endcap].
using expected_layer_pattern_type = std::array<unsigned int, 4>;

/// Mapping from detector surface identifier to bitmask position.
/// pattern_index selects one of the 4 detector partitions in the pattern.
/// layer_index selects the bit index (layer) inside the selected detector
/// partition.
struct expected_layer_mapping_entry {
  detray::geometry::identifier::value_t geo_id{0u};
  unsigned int pattern_index{0u};
  unsigned int layer_index{0u};
};

/// Default mapper that disables pattern updates.
struct null_expected_layer_mapper {
  struct result {
    bool valid{false};
    unsigned int pattern_index{0u};
    unsigned int layer_index{0u};
  };

  TRACCC_HOST_DEVICE constexpr result operator()(
      const detray::geometry::identifier&) const {
    return {};
  }
};

/// Mapper that uses a flat lookup table of geometry identifier entries.
struct expected_layer_table_mapper {
  using entry_type = expected_layer_mapping_entry;

  struct result {
    bool valid{false};
    unsigned int pattern_index{0u};
    unsigned int layer_index{0u};
  };

  const entry_type* entries{nullptr};
  std::size_t size{0u};

  TRACCC_HOST_DEVICE result
  operator()(const detray::geometry::identifier& geo_id) const {
    if (entries == nullptr || size == 0u) {
      return {};
    }

    const auto value = barcode.value();
    for (std::size_t i = 0u; i < size; ++i) {
      const auto& entry = entries[i];
      if (entry.barcode == value) {
        return {true, entry.pattern_index, entry.layer_index};
      }
    }

    return {};
  }
};

/// Collector actor for expected-layer-pattern bitmasks.
///
/// The mapper must return:
/// - valid = true if the barcode can be mapped
/// - pattern_index in [0,3]
/// - layer_index in [0,31]
template <typename mapper_t = null_expected_layer_mapper>
struct expected_layer_pattern_collector : detray::base_actor {
  using pattern_type = std::array<unsigned int, 4>;
  using mapper_type = mapper_t;
  using mapping_result_type = typename mapper_type::result;

  struct state {
    /// Target bitmask to update.
    pattern_type* pattern{nullptr};
    /// Mapping from detray barcode to (pattern index, layer index).
    mapper_type mapper{};

    /// Optional guard against repeated updates on the same surface.
    bool deduplicate_consecutive_surfaces{true};
    detray::geometry::identifier last_surface{};
    bool has_last_surface{false};

    /// Lightweight counters for debugging.
    unsigned int n_updates{0u};
    unsigned int n_skipped{0u};
  };

  template <typename propagator_state_t>
  TRACCC_HOST_DEVICE void operator()(state& actor_state,
                                     propagator_state_t& propagation) const {
    const auto& navigation = propagation.navigation();

    // Keep behavior aligned with sensitive detector element collection.
    if (!navigation.is_on_sensitive()) {
      return;
    }

    const auto& sf_desc = std::as_const(navigation).current().surface();
    const auto sf_barcode = sf_desc.identifier();

    // Avoid double counting if navigator revisits the same sensitive.
    if (actor_state.deduplicate_consecutive_surfaces &&
        actor_state.has_last_surface &&
        sf_barcode == actor_state.last_surface) {
      ++actor_state.n_skipped;
      return;
    }

    actor_state.last_surface = sf_barcode;
    actor_state.has_last_surface = true;

    if (actor_state.pattern == nullptr) {
      ++actor_state.n_skipped;
      return;
    }

    const mapping_result_type mapping = actor_state.mapper(sf_barcode);
    if (!mapping.valid || mapping.pattern_index >= 4u ||
        mapping.layer_index >=
            static_cast<unsigned int>(sizeof(unsigned int) * 8u)) {
      ++actor_state.n_skipped;
      return;
    }

    (*actor_state.pattern)[mapping.pattern_index] |=
        (1u << mapping.layer_index);
    ++actor_state.n_updates;
  }
};

}  // namespace traccc
