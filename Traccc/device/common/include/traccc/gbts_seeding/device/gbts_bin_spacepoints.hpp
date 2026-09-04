/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/device/concepts/thread_id.hpp"
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/spacepoint_collection.hpp"
#include "traccc/gbts_seeding/gbts_seeding_config.hpp"
#include "traccc/gbts_seeding/gbts_types.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>

// System include(s).
#include <bit>
#include <utility>

namespace traccc::device {

/// (Global Event Data) Payload for the @c traccc::device::gbts_bin_spacepoints
/// function
///
/// Each spacepoint is read once: it is assigned a GBTS layer (or rejected),
/// its reduced parameters are written, its eta bin's node count is bumped,
/// and its node sort key is appended to the compacted key array, all in a
/// single pass.
struct gbts_bin_spacepoints_payload {
  /// Number of spacepoints in the event
  unsigned int nSp;
  /// Number of eta bins
  unsigned int nEtaBins;
  /// All spacepoints in the event
  edm::spacepoint_collection::const_view spacepoints;
  /// All measurements in the event (used to look up surface IDs)
  edm::measurement_collection::const_view measurements;
  /// Map from detector volume index to GBTS layer index
  vecmem::data::vector_view<const short> volumeToLayerMap;
  /// Map from (volume, surface) pair to GBTS layer index (optional)
  vecmem::data::vector_view<const std::pair<unsigned int, unsigned int>>
      surfaceToLayerMap;
  /// Per-layer type code (barrel or endcap) used for cluster-width cuts
  vecmem::data::vector_view<const char> layerType;
  /// Per-layer (first eta bin, number of eta bins) pair
  vecmem::data::vector_view<const std::pair<unsigned int, unsigned int>>
      layer_info;
  /// Per-layer geometry pair used to compute eta (e.g. (rmin, zmax))
  vecmem::data::vector_view<const std::pair<float, float>> layer_geo;
  /// Output: reduced (x, y, z, cluster width) per spacepoint after filtering
  vecmem::data::vector_view<float4> reducedSP;
  /// Output: nEtaBins + 1 counters, all atomically incremented.
  vecmem::data::vector_view<unsigned int> eta_node_counter;
  /// Output: one 64-bit node sort key per accepted spacepoint
  /// Which is used to sort the nodes by (eta bin, phi, spacepoint index) before
  /// the node-making kernel runs.
  vecmem::data::vector_view<unsigned long long int> sort_keys;
  /// Output: the full spacepoint index matching each @c sort_keys slot;
  /// sorted alongside the keys by the gbts_sort_nodes launcher
  vecmem::data::vector_view<unsigned int> sort_values;
  /// Size of the volume-to-layer map (for bounds checking)
  unsigned long int volumeMapSize;
  /// Size of the surface-to-layer map (for bounds checking)
  unsigned long int surfaceMapSize;
  /// Parameters for SP filtering (passed through from config, used for tau
  /// cut if enabled)
  traccc::gbts_count_spacepoints_by_layer_params
      gbts_count_spacepoints_by_layer_params;
};

/// @brief
///
/// Per-spacepoint binning kernel: look up the GBTS layer via the
/// volume / surface map, optionally apply a cluster-width cut, and on
/// acceptance write the reduced (x, y, z, width) tuple, bump the node's eta
/// bin count and append its node sort key.
///
/// @param[in] thread_id Thread identifier for the kernel launch
/// @param[in] payload   The global memory payload
///
template <concepts::thread_id1 thread_id_t>
TRACCC_HOST_DEVICE inline void gbts_bin_spacepoints(
    const thread_id_t& thread_id, const gbts_bin_spacepoints_payload& payload);

}  // namespace traccc::device

#include "traccc/gbts_seeding/device/impl/gbts_bin_spacepoints.ipp"
