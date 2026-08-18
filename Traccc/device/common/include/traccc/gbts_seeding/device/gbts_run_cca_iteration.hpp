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
#include "traccc/gbts_seeding/gbts_types.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>

// System include(s).
#include <cstdint>

namespace traccc::device {

/// (Global Event Data) Payload for the @c
/// traccc::device::gbts_run_cca_iteration function
struct gbts_run_cca_iteration_payload {
  /// Number of edges in the compacted graph
  unsigned int nConnectedEdges;
  /// Maximum number of neighbours retained per edge
  unsigned int max_num_neighbours;
  /// Minimum path length required for an edge to be considered active
  unsigned char minLevel;
  /// Compacted graph from gbts_compress_graph
  vecmem::data::vector_view<const unsigned int> output_graph;
  /// In/out: per-edge level ping-pong buffer (2 * nConnectedEdges bytes)
  vecmem::data::vector_view<unsigned char> levels;
  /// In/out: per-edge active-flag (holds the next iter index, or -1
  /// once the edge is no longer active).
  vecmem::data::vector_view<char> active_edges;
  /// Output: longest outgoing-path summary per edge (length, next-edge)
  vecmem::data::vector_view<short2> outgoing_paths;
  /// Iteration index (0-based)
  unsigned char iter;
};

/// @brief One iteration of the cellular-automaton "longest path" relaxation.
///
/// Threads cooperatively process the current active-edge list, propagate
/// levels along the compact graph, and write the next iteration's active set
/// into the opposite ping-pong buffer (selected by iter parity).  The
/// final block to finish records the longest outgoing path summary per edge.
///
/// @param[in] thread_id Thread identifier for the kernel launch
/// @param[in] payload   The global memory payload
///
template <concepts::thread_id1 thread_id_t>
TRACCC_HOST_DEVICE inline void gbts_run_cca_iteration(
    const thread_id_t& thread_id,
    const gbts_run_cca_iteration_payload& payload);

}  // namespace traccc::device

#include "traccc/gbts_seeding/device/impl/gbts_run_cca_iteration.ipp"
