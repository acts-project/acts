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
  /// Reduced (x, y, z, cluster width) per spacepoint, in original order
  vecmem::data::vector_view<const float4> reducedSP;
  /// In/out: the nNodes node sort keys from gbts_bin_spacepoints.
  vecmem::data::vector_view<unsigned long long int> sort_keys;
  /// In/out: the spacepoint index belonging to each key, sorted alongside
  /// @c sort_keys by the kernel launcher
  vecmem::data::vector_view<unsigned int> sort_values;
  /// Output: per-node (tau_min, tau_max, r, z), written in sorted order
  vecmem::data::vector_view<float4> node_params;
  /// Output: per-node phi, written in sorted order
  vecmem::data::vector_view<float> node_phi;
  /// Output: per-sorted-slot original spacepoint index
  vecmem::data::vector_view<unsigned int> node_index;
  /// Optional tau lookup table (used iff gbts_sort_nodes_params.useTauLUT)
  vecmem::data::vector_view<const float> tau_lut;
  /// Tau-prediction cuts read by @c device::gbts_sort_nodes
  traccc::gbts_sort_nodes_params gbts_sort_nodes_params;
};

/// @brief
///
/// Gather nodes into their (eta bin, phi, spacepoint)-sorted slots and
/// pack their geometry tuple.
///
/// The kernel launcher first sorts the (sort_keys, sort_values) pairs
/// built in gbts_bin_spacepoints; thread i then reads the spacepoint index
/// at sort_values[i] and writes that spacepoint's node data at rank i
/// directly.
///
/// @param[in] thread_id Thread identifier for the kernel launch
/// @param[in] payload   The global memory payload
///
template <concepts::thread_id1 thread_id_t>
TRACCC_HOST_DEVICE inline void gbts_sort_nodes(
    const thread_id_t& thread_id, const gbts_sort_nodes_payload& payload);

}  // namespace traccc::device

#include "traccc/gbts_seeding/device/impl/gbts_sort_nodes.ipp"
