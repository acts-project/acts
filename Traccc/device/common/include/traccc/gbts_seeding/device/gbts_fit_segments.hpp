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

// System include(s).
#include <cstdint>

namespace traccc::device {

/// (Global Event Data) Payload for the @c traccc::device::gbts_fit_segments
/// function
struct gbts_fit_segments_payload {
  /// Upper bound on the number of paths (== path_store size minus
  /// terminus prefix)
  unsigned int nPaths;
  /// Number of terminus edges (path-store offset for fittable paths)
  unsigned int nTerminusEdges;
  /// Maximum number of neighbours retained per edge
  unsigned int max_num_neighbours;
  /// Minimum number of edges a path must have to be fit
  unsigned char minLevel;
  /// Reduced (x, y, z, r) per original spacepoint
  vecmem::data::vector_view<const float4> reducedSP;
  /// Compacted graph from gbts_compress_graph
  vecmem::data::vector_view<const unsigned int> output_graph;
  /// Per-path (edge index, parent path-store index or -1) entries
  vecmem::data::vector_view<const int2> path_store;
  /// Output: per-accepted-path (path_store index, level) seed proposal
  vecmem::data::vector_view<int2> seed_proposals;
  /// In/out: per-edge highest-bidder seed proposal (packed 64-bit)
  vecmem::data::vector_view<unsigned long long int> edge_bids;
  /// Output: per-seed-proposal ambiguity tag (multi-bid resolution flag)
  vecmem::data::vector_view<char> seed_ambiguity;
  /// Read-only upper bound on path indices (set by earlier kernels)
  unsigned int* nPathStoreSize;
  /// In/out: global atomic counter of accepted seed proposals
  unsigned int* nPropsCounter;
  /// Curvature / pT / chi-squared cut parameters
  traccc::gbts_fit_segments_params gbts_fit_segments_params;
  /// Maximum |z0| at the beamline for extrapolation cuts
  float max_z0;
};

/// @brief Fit each candidate path and emit seed proposals that pass quality
/// cuts.
///
/// One thread per path-store entry walks backwards from a leaf, gathers the
/// involved spacepoints, runs the helix / chi-squared fit, and on success
/// atomically claims a seed-proposal slot, bids for its leaf edge, and tags
/// ambiguity.
///
/// @param[in] thread_id Thread identifier for the kernel launch
/// @param[in] payload     The global memory payload
///
template <concepts::thread_id1 thread_id_t>
TRACCC_HOST_DEVICE inline void gbts_fit_segments(
    const thread_id_t& thread_id, const gbts_fit_segments_payload& payload);

}  // namespace traccc::device

#include "traccc/gbts_seeding/device/impl/gbts_fit_segments.ipp"
