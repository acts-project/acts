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

/// (Global Event Data) Payload for the @c traccc::device::gbts_fill_path_store
/// function
struct gbts_fill_path_store_payload {
  /// Capacity of the path store (maximum number of paths)
  unsigned int nPathsMax;
  /// Expected number of paths, only used to size the kernel launch
  unsigned int nPathsGrid;
  /// Device-side number of paths (clamped to nPathsMax by the kernel)
  vecmem::data::vector_view<const unsigned int> path_count;
  /// Number of edges in the compacted graph
  unsigned int nConnectedEdges;
  /// Maximum number of neighbours retained per edge
  unsigned int max_num_neighbours;
  /// Output: per-path (edge index, parent path-store index or -1) entries
  vecmem::data::vector_view<int2> path_store;
  /// Compacted graph (read for per-edge neighbour lookup)
  vecmem::data::vector_view<const unsigned int> output_graph;
  /// Per-edge CCA level array
  vecmem::data::vector_view<const unsigned char> levels;
  /// Per-edge (subtree path count, terminus flag) from CCA
  vecmem::data::vector_view<const int2> outgoing_paths;
  /// Inclusive prefix sum of the per-edge path counts
  vecmem::data::vector_view<const unsigned int> path_counts;
  /// Output: seed proposals, initialised to (-1, -1) per path
  vecmem::data::vector_view<int2> seed_proposals;
  /// Output: seed ambiguity flags, initialised to 0 per path
  vecmem::data::vector_view<char> seed_ambiguity;
  /// @name Segment fit of the path
  /// @{
  /// Minimum number of edges a path must have to be fit
  unsigned char minLevel;
  /// Reduced (x, y, z, w) per original spacepoint
  vecmem::data::vector_view<const float4> reducedSP;
  /// Curvature / pT / chi-squared cut parameters
  traccc::gbts_fit_segments_params gbts_fit_segments_params;
  /// Maximum |z0| at the beamline for extrapolation cuts
  float max_z0;
  /// In/out: per-edge highest seed bid (zeroed on entry)
  vecmem::data::vector_view<unsigned long long int> edge_bids;
  /// @}
};

/// @brief Fill the path store and fit every path.
///
/// The paths of a terminus edge are stored in preorder. One thread per path
/// finds its terminus edge by binary search on the path offsets, descends to
/// the edge the path ends at, fits the path and bids for that edge.
///
/// @param[in] thread_id Thread identifier for the kernel launch
/// @param[in] payload   The global memory payload
///
template <concepts::thread_id1 thread_id_t>
TRACCC_HOST_DEVICE inline void gbts_fill_path_store(
    const thread_id_t& thread_id, const gbts_fill_path_store_payload& payload);

}  // namespace traccc::device

#include "traccc/gbts_seeding/device/impl/gbts_fill_path_store.ipp"
