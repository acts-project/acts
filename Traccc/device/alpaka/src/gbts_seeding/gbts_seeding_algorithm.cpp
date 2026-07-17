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
#include "traccc/gbts_seeding/device/gbts_add_terminus_to_path_store.hpp"
#include "traccc/gbts_seeding/device/gbts_bid_seeds_for_hits.hpp"
#include "traccc/gbts_seeding/device/gbts_bin_spacepoints.hpp"
#include "traccc/gbts_seeding/device/gbts_compress_graph.hpp"
#include "traccc/gbts_seeding/device/gbts_convert_seeds.hpp"
#include "traccc/gbts_seeding/device/gbts_count_eta_phi_bins.hpp"
#include "traccc/gbts_seeding/device/gbts_count_spacepoints_by_layer.hpp"
#include "traccc/gbts_seeding/device/gbts_count_terminus_edges.hpp"
#include "traccc/gbts_seeding/device/gbts_fill_path_store.hpp"
#include "traccc/gbts_seeding/device/gbts_find_minmax_radius.hpp"
#include "traccc/gbts_seeding/device/gbts_fit_segments.hpp"
#include "traccc/gbts_seeding/device/gbts_link_graph_edges.hpp"
#include "traccc/gbts_seeding/device/gbts_make_graph_edges.hpp"
#include "traccc/gbts_seeding/device/gbts_match_graph_edges.hpp"
#include "traccc/gbts_seeding/device/gbts_prefix_sum_eta_phi_bins.hpp"
#include "traccc/gbts_seeding/device/gbts_rebid_seeds_for_edges.hpp"
#include "traccc/gbts_seeding/device/gbts_reindex_edges.hpp"
#include "traccc/gbts_seeding/device/gbts_reset_edge_bids.hpp"
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

/// Alpaka kernel for running @c traccc::device::gbts_count_spacepoints_by_layer
struct gbts_count_spacepoints_by_layer {
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc,
        const device::gbts_count_spacepoints_by_layer_payload payload) const {

        device::gbts_count_spacepoints_by_layer(details::thread_id1{acc},
                                                payload);
    }
};

/// Alpaka kernel for running @c traccc::device::gbts_bin_spacepoints
struct gbts_bin_spacepoints {
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc,
        const device::gbts_bin_spacepoints_payload payload) const {

        device::gbts_bin_spacepoints(details::thread_id1{acc}, payload);
    }
};

/// Alpaka kernel for running @c traccc::device::gbts_count_eta_phi_bins
struct gbts_count_eta_phi_bins {
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc,
        const device::gbts_count_eta_phi_bins_payload payload) const {

        device::gbts_count_eta_phi_bins(details::thread_id1{acc}, payload);
    }
};

/// Alpaka kernel for running @c traccc::device::gbts_prefix_sum_eta_phi_bins
struct gbts_prefix_sum_eta_phi_bins {
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc,
        const device::gbts_prefix_sum_eta_phi_bins_payload payload) const {

        device::gbts_prefix_sum_eta_phi_bins(details::thread_id1{acc}, payload);
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

/// Alpaka kernel for running @c traccc::device::gbts_make_graph_edges
struct gbts_make_graph_edges {
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc,
        const device::gbts_make_graph_edges_payload payload) const {

        auto& phi = ::alpaka::declareSharedVar<
            float[traccc::device::gbts_consts::node_buffer_length],
            __COUNTER__>(acc);
        auto& node_pack = ::alpaka::declareSharedVar<
            traccc::float4[traccc::device::gbts_consts::node_buffer_length],
            __COUNTER__>(acc);
        const alpaka::barrier<TAcc> barrier(&acc);

        device::gbts_make_graph_edges(
            details::thread_id1{acc}, barrier, payload,
            {vecmem::data::vector_view<float>(
                 traccc::device::gbts_consts::node_buffer_length, &phi[0]),
             vecmem::data::vector_view<traccc::float4>(
                 traccc::device::gbts_consts::node_buffer_length,
                 &node_pack[0])});
    }
};

/// Alpaka kernel for running @c traccc::device::gbts_link_graph_edges
struct gbts_link_graph_edges {
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc,
        const device::gbts_link_graph_edges_payload payload) const {

        device::gbts_link_graph_edges(details::thread_id1{acc}, payload);
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

/// Alpaka kernel for running @c traccc::device::gbts_reindex_edges
struct gbts_reindex_edges {
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc,
        const device::gbts_reindex_edges_payload payload) const {

        device::gbts_reindex_edges(details::thread_id1{acc}, payload);
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

/// Alpaka kernel for running @c traccc::device::gbts_count_terminus_edges
struct gbts_count_terminus_edges {
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc,
        const device::gbts_count_terminus_edges_payload payload) const {

        device::gbts_count_terminus_edges(details::thread_id1{acc}, payload);
    }
};

/// Alpaka kernel for running @c traccc::device::gbts_add_terminus_to_path_store
struct gbts_add_terminus_to_path_store {
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc,
        const device::gbts_add_terminus_to_path_store_payload payload) const {

        device::gbts_add_terminus_to_path_store(details::thread_id1{acc},
                                                payload);
    }
};

/// Alpaka kernel for running @c traccc::device::gbts_fill_path_store
struct gbts_fill_path_store {
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc,
        const device::gbts_fill_path_store_payload payload) const {

        auto& live_paths = ::alpaka::declareSharedVar<
            traccc::uint2[traccc::device::gbts_consts::live_path_buffer],
            __COUNTER__>(acc);
        auto& n_live_paths = ::alpaka::declareSharedVar<int, __COUNTER__>(acc);
        const alpaka::barrier<TAcc> barrier(&acc);

        device::gbts_fill_path_store(
            details::thread_id1{acc}, barrier, payload,
            {vecmem::data::vector_view<traccc::uint2>(
                 traccc::device::gbts_consts::live_path_buffer, &live_paths[0]),
             n_live_paths});
    }
};

/// Alpaka kernel for running @c traccc::device::gbts_fit_segments
struct gbts_fit_segments {
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc,
        const device::gbts_fit_segments_payload payload) const {

        device::gbts_fit_segments(details::thread_id1{acc}, payload);
    }
};

/// Alpaka kernel for running @c traccc::device::gbts_reset_edge_bids
struct gbts_reset_edge_bids {
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc,
        const device::gbts_reset_edge_bids_payload payload) const {

        device::gbts_reset_edge_bids(details::thread_id1{acc}, payload);
    }
};

/// Alpaka kernel for running @c traccc::device::gbts_rebid_seeds_for_edges
struct gbts_rebid_seeds_for_edges {
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc,
        const device::gbts_rebid_seeds_for_edges_payload payload) const {

        device::gbts_rebid_seeds_for_edges(details::thread_id1{acc}, payload);
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
        TAcc const& acc,
        const device::gbts_convert_seeds_payload payload) const {

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

void gbts_seeding_algorithm::gbts_count_spacepoints_by_layer_kernel(
    const device::gbts_count_spacepoints_by_layer_payload& payload) const {

    const unsigned int n_threads = 128;
    const unsigned int n_blocks = 1 + (payload.nSp - 1) / n_threads;
    ::alpaka::exec<Acc>(details::get_queue(queue()),
                        makeWorkDiv<Acc>(n_blocks, n_threads),
                        kernels::gbts_count_spacepoints_by_layer{}, payload);
    vecmem::device_vector<unsigned int> d_layer_sums(payload.layerCounts);
    details::inclusive_scan(details::get_queue(queue()), mr(),
                            d_layer_sums.begin(), d_layer_sums.end(),
                            d_layer_sums.begin());
}

void gbts_seeding_algorithm::gbts_bin_spacepoints_kernel(
    const device::gbts_bin_spacepoints_payload& payload) const {

    const unsigned int n_threads = 128;
    const unsigned int n_blocks = 1 + (payload.nSp - 1) / n_threads;
    ::alpaka::exec<Acc>(details::get_queue(queue()),
                        makeWorkDiv<Acc>(n_blocks, n_threads),
                        kernels::gbts_bin_spacepoints{}, payload);
}

void gbts_seeding_algorithm::gbts_count_eta_phi_bins_kernel(
    const device::gbts_count_eta_phi_bins_payload& payload) const {

    const unsigned int n_threads = 128;
    const unsigned int n_blocks = 1 + (payload.nEtaBins - 1) / n_threads;
    ::alpaka::exec<Acc>(details::get_queue(queue()),
                        makeWorkDiv<Acc>(n_blocks, n_threads),
                        kernels::gbts_count_eta_phi_bins{}, payload);
    vecmem::device_vector<unsigned int> d_eta_sums(payload.eta_node_counter);
    details::inclusive_scan(details::get_queue(queue()), mr(),
                            d_eta_sums.begin(), d_eta_sums.end(),
                            d_eta_sums.begin());
}

void gbts_seeding_algorithm::gbts_prefix_sum_eta_phi_bins_kernel(
    const device::gbts_prefix_sum_eta_phi_bins_payload& payload) const {

    const unsigned int n_threads = 128;
    const unsigned int n_blocks = 1 + (payload.nEtaBins - 1) / n_threads;
    ::alpaka::exec<Acc>(details::get_queue(queue()),
                        makeWorkDiv<Acc>(n_blocks, n_threads),
                        kernels::gbts_prefix_sum_eta_phi_bins{}, payload);
}

void gbts_seeding_algorithm::gbts_sort_nodes_kernel(
    const device::gbts_sort_nodes_payload& payload) const {

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

void gbts_seeding_algorithm::gbts_make_graph_edges_kernel(
    const device::gbts_make_graph_edges_payload& payload) const {

    const unsigned int n_threads = 128;
    const unsigned int n_blocks = payload.nUsedBinPairs;
    ::alpaka::exec<Acc>(details::get_queue(queue()),
                        makeWorkDiv<Acc>(n_blocks, n_threads),
                        kernels::gbts_make_graph_edges{}, payload);
    vecmem::device_vector<unsigned int> d_num_outgoing_edges(
        payload.num_outgoing_edges);
    details::inclusive_scan(
        details::get_queue(queue()), mr(), d_num_outgoing_edges.begin(),
        d_num_outgoing_edges.end(), d_num_outgoing_edges.begin());
}

void gbts_seeding_algorithm::gbts_link_graph_edges_kernel(
    const device::gbts_link_graph_edges_payload& payload) const {

    const unsigned int n_threads = 256;
    const unsigned int n_blocks = 1 + (payload.nEdges - 1) / n_threads;
    ::alpaka::exec<Acc>(details::get_queue(queue()),
                        makeWorkDiv<Acc>(n_blocks, n_threads),
                        kernels::gbts_link_graph_edges{}, payload);
}

void gbts_seeding_algorithm::gbts_match_graph_edges_kernel(
    const device::gbts_match_graph_edges_payload& payload) const {

    const unsigned int n_threads = 256;
    const unsigned int n_blocks = 1 + (payload.nEdges - 1) / n_threads;
    ::alpaka::exec<Acc>(details::get_queue(queue()),
                        makeWorkDiv<Acc>(n_blocks, n_threads),
                        kernels::gbts_match_graph_edges{}, payload);
}

void gbts_seeding_algorithm::gbts_reindex_edges_kernel(
    const device::gbts_reindex_edges_payload& payload) const {

    const unsigned int n_threads = 256;
    const unsigned int n_blocks = 1 + (payload.nEdges - 1) / n_threads;
    ::alpaka::exec<Acc>(details::get_queue(queue()),
                        makeWorkDiv<Acc>(n_blocks, n_threads),
                        kernels::gbts_reindex_edges{}, payload);
}

void gbts_seeding_algorithm::gbts_compress_graph_kernel(
    const device::gbts_compress_graph_payload& payload) const {

    const unsigned int n_threads = 256;
    const unsigned int n_blocks = 1 + (payload.nEdges - 1) / n_threads;
    ::alpaka::exec<Acc>(details::get_queue(queue()),
                        makeWorkDiv<Acc>(n_blocks, n_threads),
                        kernels::gbts_compress_graph{}, payload);
}

void gbts_seeding_algorithm::gbts_run_cca_iteration_kernel(
    const device::gbts_run_cca_iteration_payload& payload) const {

    const unsigned int n_threads = 128;
    const unsigned int n_blocks = 1 + (payload.nConnectedEdges - 1) / n_threads;
    ::alpaka::exec<Acc>(details::get_queue(queue()),
                        makeWorkDiv<Acc>(n_blocks, n_threads),
                        kernels::gbts_run_cca_iteration{}, payload);
}

void gbts_seeding_algorithm::gbts_count_terminus_edges_kernel(
    const device::gbts_count_terminus_edges_payload& payload) const {

    const unsigned int n_threads = 128;
    const unsigned int n_blocks = 1 + (payload.nConnectedEdges - 1) / n_threads;
    ::alpaka::exec<Acc>(details::get_queue(queue()),
                        makeWorkDiv<Acc>(n_blocks, n_threads),
                        kernels::gbts_count_terminus_edges{}, payload);
}

void gbts_seeding_algorithm::gbts_add_terminus_to_path_store_kernel(
    const device::gbts_add_terminus_to_path_store_payload& payload) const {

    const unsigned int n_threads = 128;
    const unsigned int n_blocks = 1 + (payload.nConnectedEdges - 1) / n_threads;
    ::alpaka::exec<Acc>(details::get_queue(queue()),
                        makeWorkDiv<Acc>(n_blocks, n_threads),
                        kernels::gbts_add_terminus_to_path_store{}, payload);
}

void gbts_seeding_algorithm::gbts_fill_path_store_kernel(
    const device::gbts_fill_path_store_payload& payload) const {

    const unsigned int n_threads = 128;
    const unsigned int pathsPerTerminus =
        1 + (payload.nPaths - 1) / payload.nTerminusEdges;
    const unsigned int terminusPerBlock = std::min(
        n_threads, 1 + (traccc::device::gbts_consts::live_path_buffer - 1) /
                           pathsPerTerminus);
    const unsigned int n_blocks =
        1 + (payload.nTerminusEdges - 1) / terminusPerBlock;
    ::alpaka::exec<Acc>(details::get_queue(queue()),
                        makeWorkDiv<Acc>(n_blocks, n_threads),
                        kernels::gbts_fill_path_store{}, payload);
}

void gbts_seeding_algorithm::gbts_fit_segments_kernel(
    const device::gbts_fit_segments_payload& payload) const {

    const unsigned int n_threads = 128;
    const unsigned int n_blocks = 1 + (payload.nPaths - 1) / n_threads;
    ::alpaka::exec<Acc>(details::get_queue(queue()),
                        makeWorkDiv<Acc>(n_blocks, n_threads),
                        kernels::gbts_fit_segments{}, payload);
}

void gbts_seeding_algorithm::gbts_reset_edge_bids_kernel(
    const device::gbts_reset_edge_bids_payload& payload) const {

    const unsigned int n_threads = 128;
    const unsigned int n_blocks = 1 + (payload.nProps - 1) / n_threads;
    ::alpaka::exec<Acc>(details::get_queue(queue()),
                        makeWorkDiv<Acc>(n_blocks, n_threads),
                        kernels::gbts_reset_edge_bids{}, payload);
}

void gbts_seeding_algorithm::gbts_rebid_seeds_for_edges_kernel(
    const device::gbts_rebid_seeds_for_edges_payload& payload) const {

    const unsigned int n_threads = 128;
    const unsigned int n_blocks = 1 + (payload.nProps - 1) / n_threads;
    ::alpaka::exec<Acc>(details::get_queue(queue()),
                        makeWorkDiv<Acc>(n_blocks, n_threads),
                        kernels::gbts_rebid_seeds_for_edges{}, payload);
}

void gbts_seeding_algorithm::gbts_bid_seeds_for_hits_kernel(
    const device::gbts_bid_seeds_for_hits_payload& payload) const {

    const unsigned int n_threads = 128;
    const unsigned int n_blocks = 1 + (payload.nProps - 1) / n_threads;
    ::alpaka::exec<Acc>(details::get_queue(queue()),
                        makeWorkDiv<Acc>(n_blocks, n_threads),
                        kernels::gbts_bid_seeds_for_hits{}, payload);
}

void gbts_seeding_algorithm::gbts_convert_seeds_kernel(
    const device::gbts_convert_seeds_payload& payload) const {

    const unsigned int n_threads = 128;
    const unsigned int n_blocks = 1 + (payload.nProps - 1) / n_threads;
    ::alpaka::exec<Acc>(details::get_queue(queue()),
                        makeWorkDiv<Acc>(n_blocks, n_threads),
                        kernels::gbts_convert_seeds{}, payload);
}

}  // namespace traccc::alpaka
