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

/// Number of relaxation sweeps (and of per-sweep change flags)
inline constexpr unsigned int gbts_run_cca_max_sweeps =
    gbts_consts::max_seed_candidate_length + 1u;

/// (Global Event Data) Payload for the @c
/// traccc::device::gbts_run_cca_iteration function
struct gbts_run_cca_iteration_payload {
  /// Number of edges in the compacted graph
  unsigned int nConnectedEdges;
  /// Maximum number of neighbours retained per edge
  unsigned int max_num_neighbours;
  /// Compacted graph from gbts_compress_graph
  vecmem::data::vector_view<const unsigned int> output_graph;
  /// In/out: per-edge level, the longest path from the edge in edges,
  /// initialised to 1
  vecmem::data::vector_view<unsigned char> levels;
  /// Output: per-edge (subtree path count, terminus flag). The flag is
  /// written by gbts_finish_cca
  vecmem::data::vector_view<int2> outgoing_paths;
  /// Sweep index
  unsigned char iter;
  /// Per-sweep "something changed" flags (gbts_run_cca_max_sweeps, zeroed
  /// per event)
  vecmem::data::vector_view<unsigned int> changed;
};

/// @brief One sweep of the longest-path relaxation (the "CCA").
///
/// A sweep recomputes every edge's level (1 + the maximum neighbour level)
/// and subtree path count in place. A sweep after one without changes
/// returns at once. gbts_finish_cca runs after the last sweep.
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
