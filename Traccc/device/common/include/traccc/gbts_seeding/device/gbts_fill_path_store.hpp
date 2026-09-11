/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/device/concepts/barrier.hpp"
#include "traccc/device/concepts/thread_id.hpp"
#include "traccc/gbts_seeding/gbts_types.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>

namespace traccc::device {

/// (Global Event Data) Payload for the @c traccc::device::gbts_fill_path_store
/// function
struct gbts_fill_path_store_payload {
  /// Number of terminus edges
  unsigned int nTerminusEdges;
  /// Maximum number of neighbours retained per edge
  unsigned int max_num_neighbours;
  /// Total number of paths
  unsigned int nPaths;
  /// In/out: per-path (edge index, parent path-store index or -1) entries
  vecmem::data::vector_view<int2> path_store;
  /// Compacted graph (read for per-edge neighbour lookup)
  vecmem::data::vector_view<const unsigned int> output_graph;
  /// Per-edge CCA level array
  vecmem::data::vector_view<const unsigned char> levels;
  /// In/out: global atomic write cursor into path_store
  unsigned int* nPathStoreSizeCounter;
};

/// (Shared Event Data) Payload for the @c traccc::device::gbts_fill_path_store
/// function
///
/// Shared-memory scratch for gbts_fill_path_store: the block-local stack of
/// live paths being walked and its running length.
struct gbts_fill_path_store_shared_payload {
  /// Shared-mem frontier of in-flight paths
  vecmem::data::vector_view<traccc::uint2> live_paths;
  /// Shared-mem frontier size
  int& n_live_paths;
};

/// @brief Walk each terminus edge backwards along live levels, growing the path
/// store.
///
/// One block expands a fixed number of terminus seeds at once using a shared
/// "live paths" frontier.  At each step the block reads neighbour edges that
/// match the next-lower level, atomically reserves slots in the path store,
/// and continues until all paths reach the graph boundary or until the
/// frontier is empty.
///
/// @param[in] thread_id          Thread/block identifier (one block/task)
/// @param[in] barrier            Block-wide barrier
/// @param[in,out] payload        The global memory payload
/// @param[in,out] shared_payload The shared memory payload
///
template <concepts::thread_id1 thread_id_t, concepts::barrier barrier_t>
TRACCC_HOST_DEVICE inline void gbts_fill_path_store(
    const thread_id_t& thread_id, const barrier_t& barrier,
    const gbts_fill_path_store_payload& payload,
    const gbts_fill_path_store_shared_payload& shared_payload);

}  // namespace traccc::device

#include "traccc/gbts_seeding/device/impl/gbts_fill_path_store.ipp"
