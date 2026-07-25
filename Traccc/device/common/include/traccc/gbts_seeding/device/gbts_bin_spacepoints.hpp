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
#include "traccc/gbts_seeding/gbts_types.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>

// System include(s).
#include <utility>

namespace traccc::device {

/// (Global Event Data) Payload for the @c traccc::device::gbts_bin_spacepoints
/// function
///
/// Each spacepoint is read once and scattered into its layer slot, its
/// eta/phi bin indices are computed, and the (eta, phi) histogram is
/// bumped, all in a single pass.
struct gbts_bin_spacepoints_payload {
  /// Number of spacepoints in the event
  unsigned int nSp;
  /// Number of phi bins per eta slice
  unsigned int nPhiBins;
  /// Output: per-spacepoint (x, y, z, r) bin parameters, layer-ordered
  vecmem::data::vector_view<float4> sp_params;
  /// Input: reduced (x, y, z, r) per spacepoint from
  /// gbts_count_spacepoints_by_layer
  vecmem::data::vector_view<const float4> reducedSP;
  /// In/out: per-layer running write cursors (decremented as each SP lands)
  vecmem::data::vector_view<unsigned int> layerCounts;
  /// GBTS layer assignment for each spacepoint
  vecmem::data::vector_view<const unsigned short> spacepointsLayer;
  /// Output: layer-ordered index back to the original spacepoint slot
  vecmem::data::vector_view<unsigned int> original_sp_idx;
  /// Per-layer (first eta bin, number of eta bins) pair
  vecmem::data::vector_view<const std::pair<unsigned int, unsigned int>>
      layer_info;
  /// Per-layer geometry pair used to compute eta (e.g. (rmin, zmax))
  vecmem::data::vector_view<const std::pair<float, float>> layer_geo;
  /// Output: global eta-bin index assigned to each node
  vecmem::data::vector_view<unsigned int> node_eta_index;
  /// Output: phi-bin index assigned to each node
  vecmem::data::vector_view<unsigned int> node_phi_index;
  /// Output: flat (eta, phi) histogram, atomically incremented per node
  vecmem::data::vector_view<unsigned int> eta_phi_histo;
};

/// @brief Per-spacepoint binning kernel: atomically claim a layer-ordered slot,
/// compute the node's eta- and phi-bin indices, and bump the (eta, phi)
/// histogram bucket -- all from a single read of the source spacepoint.
///
/// @param[in] thread_id Thread identifier for the kernel launch
/// @param[in] payload   The global memory payload
///
template <concepts::thread_id1 thread_id_t>
TRACCC_HOST_DEVICE inline void gbts_bin_spacepoints(
    const thread_id_t& thread_id, const gbts_bin_spacepoints_payload& payload);

}  // namespace traccc::device

#include "traccc/gbts_seeding/device/impl/gbts_bin_spacepoints.ipp"
