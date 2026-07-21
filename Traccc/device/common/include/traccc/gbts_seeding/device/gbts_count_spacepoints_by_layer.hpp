/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/device/concepts/thread_id.hpp"
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/spacepoint_collection.hpp"
#include "traccc/gbts_seeding/gbts_types.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>

// System include(s).
#include <utility>

namespace traccc::device {

/// (Global Event Data) Payload for the @c
/// traccc::device::gbts_count_spacepoints_by_layer function
struct gbts_count_spacepoints_by_layer_payload {
  /// Number of spacepoints in the event
  unsigned int nSp;
  /// All spacepoints in the event
  edm::spacepoint_collection::const_view spacepoints;
  /// All measurements in the event (used to look up surface IDs)
  edm::measurement_collection::const_view measurements;
  /// Map from detector volume index to GBTS layer index
  vecmem::data::vector_view<const short> volumeToLayerMap;
  /// Map from (volume, surface) pair to GBTS layer index (optional)
  vecmem::data::vector_view<const std::pair<unsigned int, unsigned int>>
      surfaceToLayerMap;
  /// Per-layer type code (barrel/endcap/etc.) used for cluster-width cuts
  vecmem::data::vector_view<const char> layerType;
  /// Output: reduced (x, y, z, r) per spacepoint after filtering
  vecmem::data::vector_view<float4> reducedSP;
  /// Output: per-layer spacepoint counts (atomically incremented)
  vecmem::data::vector_view<unsigned int> layerCounts;
  /// Output: GBTS layer index assigned to each kept spacepoint
  vecmem::data::vector_view<unsigned short> spacepointsLayer;
  /// Size of the volume-to-layer map (for bounds checking)
  unsigned long int volumeMapSize;
  /// Size of the surface-to-layer map (for bounds checking)
  unsigned long int surfaceMapSize;
  /// Parameters for SP counting (passed through from config, used for tau
  /// cut if enabled)
  traccc::gbts_count_spacepoints_by_layer_params
      gbts_count_spacepoints_by_layer_params;
};

/// @brief Count and tag spacepoints by GBTS layer, producing the reduced SP
/// view.
///
/// Each thread inspects one spacepoint, looks up its GBTS layer via the
/// volume / surface map, optionally applies a cluster-width cut, and on
/// acceptance atomically increments the per-layer count and writes the
/// reduced (x, y, z, r) tuple plus the assigned layer index.
///
/// @param[in] thread_id Thread identifier for the kernel launch
/// @param[in] payload   The global memory payload
///
template <concepts::thread_id1 thread_id_t>
TRACCC_HOST_DEVICE inline void gbts_count_spacepoints_by_layer(
    const thread_id_t& thread_id,
    const gbts_count_spacepoints_by_layer_payload& payload);

}  // namespace traccc::device

#include "traccc/gbts_seeding/device/impl/gbts_count_spacepoints_by_layer.ipp"
