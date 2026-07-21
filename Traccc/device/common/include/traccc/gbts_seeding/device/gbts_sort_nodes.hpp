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
#include "traccc/gbts_seeding/gbts_seeding_config.hpp"
#include "traccc/gbts_seeding/gbts_types.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>

namespace traccc::device {

/// (Global Event Data) Payload for the @c traccc::device::gbts_sort_nodes
/// function
struct gbts_sort_nodes_payload {
  /// Total number of GBTS nodes
  unsigned int nNodes;
  /// Number of phi bins per eta slice
  unsigned int nPhiBins;
  /// Layer-ordered (x, y, z, r) per spacepoint
  vecmem::data::vector_view<const float4> sp_params;
  /// Eta-bin index per node
  vecmem::data::vector_view<const unsigned int> node_eta_index;
  /// Phi-bin index per node
  vecmem::data::vector_view<const unsigned int> node_phi_index;
  /// In/out: per-(eta, phi) write cursor (atomically advanced)
  vecmem::data::vector_view<unsigned int> phi_cusums;
  /// Output: per-node (tau_min, tau_max, r, z), written in sorted order
  vecmem::data::vector_view<float4> node_params;
  /// Output: per-node phi, written in sorted order
  vecmem::data::vector_view<float> node_phi;
  /// Output: per-sorted-slot original layer-ordered spacepoint index
  vecmem::data::vector_view<unsigned int> node_index;
  /// Map from layer-ordered SP index to the original SP slot
  vecmem::data::vector_view<const unsigned int> original_sp_idx;
  /// Optional tau lookup table (used iff gbts_sort_nodes_params.useTauLUT)
  vecmem::data::vector_view<const float> tau_lut;
  /// Tau-prediction cuts read by @c device::gbts_sort_nodes
  traccc::gbts_sort_nodes_params gbts_sort_nodes_params;
};

/// @brief Scatter nodes into (eta, phi)-sorted slots and pack their geometry
/// tuple.
///
/// Each thread picks one node, atomically increments its (eta, phi) write
/// cursor, packs the geometry / kinematic tuple into the node parameters at
/// that slot, and stores the original SP index in the node index map.
///
/// @param[in] thread_id Thread identifier for the kernel launch
/// @param[in] payload   The global memory payload
///
template <concepts::thread_id1 thread_id_t>
TRACCC_HOST_DEVICE inline void gbts_sort_nodes(
    const thread_id_t& thread_id, const gbts_sort_nodes_payload& payload);

}  // namespace traccc::device

#include "traccc/gbts_seeding/device/impl/gbts_sort_nodes.ipp"
