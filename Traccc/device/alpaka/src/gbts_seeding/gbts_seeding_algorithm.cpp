/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/alpaka/gbts_seeding/gbts_seeding_algorithm.hpp"

#include "../utils/barrier.hpp"
#include "../utils/get_queue.hpp"
#include "../utils/parallel_algorithms.hpp"
#include "../utils/thread_id.hpp"
#include "../utils/utils.hpp"

// Project include(s).
#include "traccc/gbts_seeding/device/gbts_bid_seeds_for_hits.hpp"
#include "traccc/gbts_seeding/device/gbts_bin_spacepoints.hpp"
#include "traccc/gbts_seeding/device/gbts_build_edge_work_list.hpp"
#include "traccc/gbts_seeding/device/gbts_compress_graph.hpp"
#include "traccc/gbts_seeding/device/gbts_convert_seeds.hpp"
#include "traccc/gbts_seeding/device/gbts_count_graph_edges.hpp"
#include "traccc/gbts_seeding/device/gbts_count_paths.hpp"
#include "traccc/gbts_seeding/device/gbts_fill_graph_edges.hpp"
#include "traccc/gbts_seeding/device/gbts_fill_path_store.hpp"
#include "traccc/gbts_seeding/device/gbts_find_minmax_radius.hpp"
#include "traccc/gbts_seeding/device/gbts_finish_cca.hpp"
#include "traccc/gbts_seeding/device/gbts_match_graph_edges.hpp"
#include "traccc/gbts_seeding/device/gbts_run_cca_iteration.hpp"
#include "traccc/gbts_seeding/device/gbts_sort_nodes.hpp"
#include "traccc/gbts_seeding/gbts_types.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>
#include <vecmem/containers/device_vector.hpp>

// System include(s).
#include <algorithm>

namespace traccc::alpaka {

namespace kernels {

// ---------------------------------------------------------------------------
// Stage 1 — nodes-making kernels
// ---------------------------------------------------------------------------

/// Alpaka kernel for running @c traccc::device::gbts_bin_spacepoints
struct gbts_bin_spacepoints {
  template <typename TAcc>
  ALPAKA_FN_ACC void operator()(
      TAcc const& acc,
      const device::gbts_bin_spacepoints_payload payload) const {
    device::gbts_bin_spacepoints(details::thread_id1{acc}, payload);
  }
};

/// Alpaka kernel for running @c traccc::device::gbts_sort_nodes
struct gbts_sort_nodes {
  template <typename TAcc>
  ALPAKA_FN_ACC void operator()(
      TAcc const& acc, const device::gbts_sort_nodes_payload payload) const {
    device::gbts_sort_nodes(details::thread_id1{acc}, payload);
  }
};

/// Alpaka kernel for running @c traccc::device::gbts_find_minmax_radius
struct gbts_find_minmax_radius {
  template <typename TAcc>
  ALPAKA_FN_ACC void operator()(
      TAcc const& acc,
      const device::gbts_find_minmax_radius_payload payload) const {
    device::gbts_find_minmax_radius(details::thread_id1{acc}, payload);
  }
};

// ---------------------------------------------------------------------------
// Stage 2 — graph-making kernels
// ---------------------------------------------------------------------------

/// Alpaka kernel for running @c traccc::device::gbts_build_edge_work_list
struct gbts_build_edge_work_list {
  template <typename TAcc>
  ALPAKA_FN_ACC void operator()(
      TAcc const& acc,
      const device::gbts_build_edge_work_list_payload payload) const {
    auto& scratch = ::alpaka::declareSharedVar<
        unsigned int[device::gbts_build_edge_work_list_block_size],
        __COUNTER__>(acc);
    const alpaka::barrier<TAcc> barrier(&acc);
    device::gbts_build_edge_work_list(
        details::thread_id1{acc}, barrier, payload,
        {vecmem::data::vector_view<unsigned int>(
            device::gbts_build_edge_work_list_block_size, &scratch[0])});
  }
};

/// Alpaka kernel for running @c traccc::device::gbts_count_graph_edges
struct gbts_count_graph_edges {
  template <typename TAcc>
  ALPAKA_FN_ACC void operator()(
      TAcc const& acc,
      const device::gbts_count_graph_edges_payload payload) const {
    device::gbts_count_graph_edges(details::thread_id1{acc}, payload);
  }
};

/// Alpaka kernel for running @c traccc::device::gbts_fill_graph_edges
struct gbts_fill_graph_edges {
  template <typename TAcc>
  ALPAKA_FN_ACC void operator()(
      TAcc const& acc,
      const device::gbts_fill_graph_edges_payload payload) const {
    device::gbts_fill_graph_edges(details::thread_id1{acc}, payload);
  }
};

/// Alpaka kernel for running @c traccc::device::gbts_match_graph_edges
struct gbts_match_graph_edges {
  template <typename TAcc>
  ALPAKA_FN_ACC void operator()(
      TAcc const& acc,
      const device::gbts_match_graph_edges_payload payload) const {
    device::gbts_match_graph_edges(details::thread_id1{acc}, payload);
  }
};

/// Alpaka kernel for running @c traccc::device::gbts_compress_graph
struct gbts_compress_graph {
  template <typename TAcc>
  ALPAKA_FN_ACC void operator()(
      TAcc const& acc,
      const device::gbts_compress_graph_payload payload) const {
    device::gbts_compress_graph(details::thread_id1{acc}, payload);
  }
};

// ---------------------------------------------------------------------------
// Stage 3 — graph-processing kernels
// ---------------------------------------------------------------------------

/// Alpaka kernel for running @c traccc::device::gbts_run_cca_iteration
struct gbts_run_cca_iteration {
  template <typename TAcc>
  ALPAKA_FN_ACC void operator()(
      TAcc const& acc,
      const device::gbts_run_cca_iteration_payload payload) const {
    device::gbts_run_cca_iteration(details::thread_id1{acc}, payload);
  }
};

/// Alpaka kernel for running @c traccc::device::gbts_finish_cca
struct gbts_finish_cca {
  template <typename TAcc>
  ALPAKA_FN_ACC void operator()(
      TAcc const& acc, const device::gbts_finish_cca_payload payload) const {
    device::gbts_finish_cca(details::thread_id1{acc}, payload);
  }
};

/// Alpaka kernel for running @c traccc::device::gbts_count_paths
struct gbts_count_paths {
  template <typename TAcc>
  ALPAKA_FN_ACC void operator()(
      TAcc const& acc, const device::gbts_count_paths_payload payload) const {
    device::gbts_count_paths(details::thread_id1{acc}, payload);
  }
};

/// Alpaka kernel for running @c traccc::device::gbts_fill_path_store
struct gbts_fill_path_store {
  template <typename TAcc>
  ALPAKA_FN_ACC void operator()(
      TAcc const& acc,
      const device::gbts_fill_path_store_payload payload) const {
    device::gbts_fill_path_store(details::thread_id1{acc}, payload);
  }
};

/// Alpaka kernel for running @c traccc::device::gbts_bid_seeds_for_hits
struct gbts_bid_seeds_for_hits {
  template <typename TAcc>
  ALPAKA_FN_ACC void operator()(
      TAcc const& acc,
      const device::gbts_bid_seeds_for_hits_payload payload) const {
    device::gbts_bid_seeds_for_hits(details::thread_id1{acc}, payload);
  }
};

/// Alpaka kernel for running @c traccc::device::gbts_convert_seeds
struct gbts_convert_seeds {
  template <typename TAcc>
  ALPAKA_FN_ACC void operator()(
      TAcc const& acc, const device::gbts_convert_seeds_payload payload) const {
    device::gbts_convert_seeds(details::thread_id1{acc}, payload);
  }
};

}  // namespace kernels

// ===========================================================================
// gbts_seeding_algorithm: kernel launchers
// ===========================================================================

gbts_seeding_algorithm::gbts_seeding_algorithm(
    const gbts_seedfinder_config& cfg, const memory_resource& mr,
    const vecmem::copy& copy, alpaka::queue& q,
    std::unique_ptr<const Logger> logger)
    : device::gbts_seeding_algorithm(cfg, mr, copy, std::move(logger)),
      alpaka::algorithm_base{q} {}

void gbts_seeding_algorithm::gbts_bin_spacepoints_kernel(
    const device::gbts_bin_spacepoints_payload& payload) const {
  const unsigned int n_threads = 128;
  const unsigned int n_blocks = 1 + (payload.nSp - 1) / n_threads;
  ::alpaka::exec<Acc>(details::get_queue(queue()),
                      makeWorkDiv<Acc>(n_blocks, n_threads),
                      kernels::gbts_bin_spacepoints{}, payload);

  // Turn the per-bin node counts into the node offsets.
  vecmem::device_vector<unsigned int> d_eta_node_counter(
      payload.eta_node_counter);
  details::exclusive_scan(details::get_queue(queue()), mr(),
                          d_eta_node_counter.begin(), d_eta_node_counter.end(),
                          d_eta_node_counter.begin());
}

void gbts_seeding_algorithm::gbts_sort_nodes_kernel(
    const device::gbts_sort_nodes_payload& payload) const {
  // Order the nodes by their (eta bin, phi, spacepoint index bits) keys,
  // carrying the full spacepoint index along as the value.
  details::sort_by_key(
      details::get_queue(queue()), mr(), payload.sort_keys.ptr(),
      payload.sort_keys.ptr() + payload.nNodes, payload.sort_values.ptr());

  const unsigned int n_threads = 256;
  const unsigned int n_blocks = 1 + (payload.nNodes - 1) / n_threads;
  ::alpaka::exec<Acc>(details::get_queue(queue()),
                      makeWorkDiv<Acc>(n_blocks, n_threads),
                      kernels::gbts_sort_nodes{}, payload);
}

void gbts_seeding_algorithm::gbts_find_minmax_radius_kernel(
    const device::gbts_find_minmax_radius_payload& payload) const {
  const unsigned int n_threads = 128;
  const unsigned int n_blocks = 1 + (payload.nEtaBins - 1) / n_threads;
  ::alpaka::exec<Acc>(details::get_queue(queue()),
                      makeWorkDiv<Acc>(n_blocks, n_threads),
                      kernels::gbts_find_minmax_radius{}, payload);
}

void gbts_seeding_algorithm::gbts_build_edge_work_list_kernel(
    const device::gbts_build_edge_work_list_payload& payload) const {
  ::alpaka::exec<Acc>(
      details::get_queue(queue()),
      makeWorkDiv<Acc>(1u, device::gbts_build_edge_work_list_block_size),
      kernels::gbts_build_edge_work_list{}, payload);
}

void gbts_seeding_algorithm::gbts_count_graph_edges_kernel(
    const device::gbts_count_graph_edges_payload& payload) const {
  // One thread per inner node of a chunk; the blocks stride over the work
  // items.
  const unsigned int n_threads = device::gbts_consts::edge_chunk_size;
  const unsigned int n_blocks =
      std::min(payload.nWorkMax, device::gbts_count_graph_edges_max_blocks);
  ::alpaka::exec<Acc>(details::get_queue(queue()),
                      makeWorkDiv<Acc>(n_blocks, n_threads),
                      kernels::gbts_count_graph_edges{}, payload);

  // Turn the per-node counts into the edge buckets.
  vecmem::device_vector<unsigned int> d_num_outgoing_edges(
      payload.num_outgoing_edges);
  details::inclusive_scan(
      details::get_queue(queue()), mr(), d_num_outgoing_edges.begin(),
      d_num_outgoing_edges.end(), d_num_outgoing_edges.begin());
}

void gbts_seeding_algorithm::gbts_fill_graph_edges_kernel(
    const device::gbts_fill_graph_edges_payload& payload) const {
  const unsigned int n_threads = device::gbts_consts::edge_chunk_size;
  const unsigned int n_blocks =
      std::min(payload.nWorkMax, device::gbts_fill_graph_edges_max_blocks);
  ::alpaka::exec<Acc>(details::get_queue(queue()),
                      makeWorkDiv<Acc>(n_blocks, n_threads),
                      kernels::gbts_fill_graph_edges{}, payload);
}

void gbts_seeding_algorithm::gbts_match_graph_edges_kernel(
    const device::gbts_match_graph_edges_payload& payload) const {
  const unsigned int n_threads = 256;
  const unsigned int n_blocks = 1u + (payload.nEdgesMax - 1u) / n_threads;
  ::alpaka::exec<Acc>(details::get_queue(queue()),
                      makeWorkDiv<Acc>(n_blocks, n_threads),
                      kernels::gbts_match_graph_edges{}, payload);

  // Compact the kept edges with a prefix sum over their 0/1 flags.
  details::inclusive_scan(details::get_queue(queue()), mr(), payload.kept.ptr(),
                          payload.kept.ptr() + payload.nEdgesMax,
                          payload.reIndexer.ptr());
}

void gbts_seeding_algorithm::gbts_compress_graph_kernel(
    const device::gbts_compress_graph_payload& payload) const {
  const unsigned int n_threads = 256;
  const unsigned int n_blocks = 1u + (payload.nEdgesMax - 1u) / n_threads;
  ::alpaka::exec<Acc>(details::get_queue(queue()),
                      makeWorkDiv<Acc>(n_blocks, n_threads),
                      kernels::gbts_compress_graph{}, payload);
}

void gbts_seeding_algorithm::gbts_run_cca_iteration_kernel(
    const device::gbts_run_cca_iteration_payload& payload) const {
  const unsigned int n_threads = 256;
  const unsigned int n_blocks = 1 + (payload.nConnectedEdges - 1) / n_threads;
  ::alpaka::exec<Acc>(details::get_queue(queue()),
                      makeWorkDiv<Acc>(n_blocks, n_threads),
                      kernels::gbts_run_cca_iteration{}, payload);
}

void gbts_seeding_algorithm::gbts_finish_cca_kernel(
    const device::gbts_finish_cca_payload& payload) const {
  const unsigned int n_threads = 256;
  const unsigned int n_blocks = 1 + (payload.nConnectedEdges - 1) / n_threads;
  ::alpaka::exec<Acc>(details::get_queue(queue()),
                      makeWorkDiv<Acc>(n_blocks, n_threads),
                      kernels::gbts_finish_cca{}, payload);
}

void gbts_seeding_algorithm::gbts_count_paths_kernel(
    const device::gbts_count_paths_payload& payload) const {
  const unsigned int n_threads = 128;
  const unsigned int n_blocks = 1 + (payload.nConnectedEdges - 1) / n_threads;
  ::alpaka::exec<Acc>(details::get_queue(queue()),
                      makeWorkDiv<Acc>(n_blocks, n_threads),
                      kernels::gbts_count_paths{}, payload);
  // Path offsets of the path store.
  vecmem::device_vector<unsigned int> d_path_counts(payload.path_counts);
  details::inclusive_scan(
      details::get_queue(queue()), mr(), d_path_counts.begin(),
      d_path_counts.begin() + payload.nConnectedEdges, d_path_counts.begin());
}

void gbts_seeding_algorithm::gbts_fill_path_store_kernel(
    const device::gbts_fill_path_store_payload& payload) const {
  const unsigned int n_threads = 128;
  const unsigned int n_blocks = 1 + (payload.nPathsGrid - 1) / n_threads;
  ::alpaka::exec<Acc>(details::get_queue(queue()),
                      makeWorkDiv<Acc>(n_blocks, n_threads),
                      kernels::gbts_fill_path_store{}, payload);
}

void gbts_seeding_algorithm::gbts_bid_seeds_for_hits_kernel(
    const device::gbts_bid_seeds_for_hits_payload& payload) const {
  const unsigned int n_threads = 128;
  const unsigned int n_blocks = 1 + (payload.nPathsGrid - 1) / n_threads;
  ::alpaka::exec<Acc>(details::get_queue(queue()),
                      makeWorkDiv<Acc>(n_blocks, n_threads),
                      kernels::gbts_bid_seeds_for_hits{}, payload);
}

void gbts_seeding_algorithm::gbts_convert_seeds_kernel(
    const device::gbts_convert_seeds_payload& payload) const {
  const unsigned int n_threads = 128;
  const unsigned int n_blocks = 1 + (payload.nPathsGrid - 1) / n_threads;
  ::alpaka::exec<Acc>(details::get_queue(queue()),
                      makeWorkDiv<Acc>(n_blocks, n_threads),
                      kernels::gbts_convert_seeds{}, payload);
}

}  // namespace traccc::alpaka
