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

namespace traccc::device {

/// (Global Event Data) Payload for the @c traccc::device::gbts_compress_graph
/// function
struct gbts_compress_graph_payload {
  /// Number of original (uncompressed) edges
  unsigned int nEdges;
  /// Maximum number of neighbours retained per edge
  unsigned int nMaxNei;
  /// Sorted-slot to original spacepoint index map
  vecmem::data::vector_view<const unsigned int> orig_node_index;
  /// (src, dst) node indices per edge
  vecmem::data::vector_view<const uint2> edge_nodes;
  /// Number of accepted neighbours per edge
  vecmem::data::vector_view<const unsigned char> num_neighbours;
  /// Neighbour edge indices per edge (nMaxNei per edge, flat)
  vecmem::data::vector_view<const unsigned int> neighbours;
  /// Old-edge to compacted-edge index map
  vecmem::data::vector_view<const int> reIndexer;
  /// Output: compacted graph in row-major layout; each edge owns a block
  /// of edge_size = 2 + 1 + nMaxNei ints (node1, node2, nNei,
  /// nei0..neiN-1).
  vecmem::data::vector_view<unsigned int> output_graph;
};

/// @brief Pack kept edges into the compact "output graph" layout.
///
/// Each thread processes one original edge; if it survived re-indexing, the
/// thread writes a record at its compact slot containing the
/// source/destination original-SP indices, the neighbour count, and up to
/// nMaxNei remapped neighbour indices.
///
/// @param[in] thread_id Thread identifier for the kernel launch
/// @param[in] payload   The global memory payload
///
template <concepts::thread_id1 thread_id_t>
TRACCC_HOST_DEVICE inline void gbts_compress_graph(
    const thread_id_t& thread_id, const gbts_compress_graph_payload& payload);

}  // namespace traccc::device

#include "traccc/gbts_seeding/device/impl/gbts_compress_graph.ipp"
