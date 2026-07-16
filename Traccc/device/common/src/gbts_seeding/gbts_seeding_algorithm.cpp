/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/gbts_seeding/device/gbts_seeding_algorithm.hpp"

#include "traccc/gbts_seeding/gbts_seeding_config.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_buffer.hpp>
#include <vecmem/containers/vector.hpp>

// System include(s).
#include <algorithm>
#include <cmath>
#include <utility>

namespace traccc::device {

// Stage 1:
// Bin the spacepoints by layer in eta and phi,
// compute the node parameters (x, y, z, w)
// and the bin-wise min/max radius for the graph-building cuts.
auto gbts_seeding_algorithm::make_nodes(
    const edm::spacepoint_collection::const_view& spacepoints,
    const edm::measurement_collection::const_view& measurements) const
    -> node_making_output {

    const gbts_seedfinder_config& cfg = m_config;
    const unsigned int nSp = copy().get_size(spacepoints);

    // 0. Bin spacepoints by the mapping supplied to config.surfaceToLayerMap.
    vecmem::data::vector_buffer<unsigned int> layerCounts_buf(cfg.nLayers + 1,
                                                              mr().main);
    copy().memset(layerCounts_buf, 0)->ignore();

    vecmem::data::vector_buffer<float4> reducedSP_buf(nSp, mr().main);
    copy().setup(reducedSP_buf)->ignore();

    vecmem::data::vector_buffer<unsigned short> spacepointsLayer_buf(nSp,
                                                                     mr().main);
    copy().setup(spacepointsLayer_buf)->ignore();

    vecmem::data::vector_buffer<short> volumeToLayerMap_buf(
        static_cast<unsigned int>(cfg.volumeToLayerMap.size()), mr().main);
    copy().setup(volumeToLayerMap_buf)->ignore();
    copy()(vecmem::get_data(cfg.volumeToLayerMap), volumeToLayerMap_buf)
        ->ignore();

    vecmem::data::vector_buffer<std::pair<unsigned int, unsigned int>>
        surfaceToLayerMap_buf;
    if (!cfg.surfaceToLayerMap.empty()) {
        surfaceToLayerMap_buf =
            vecmem::data::vector_buffer<std::pair<unsigned int, unsigned int>>(
                static_cast<unsigned int>(cfg.surfaceToLayerMap.size()),
                mr().main);
        copy().setup(surfaceToLayerMap_buf)->ignore();
        copy()(vecmem::get_data(cfg.surfaceToLayerMap), surfaceToLayerMap_buf)
            ->ignore();
    }

    vecmem::data::vector_buffer<char> layerType_buf(cfg.nLayers, mr().main);
    copy().setup(layerType_buf)->ignore();
    copy()(vecmem::get_data(cfg.layerInfo.type), layerType_buf)->ignore();
    gbts_count_spacepoints_by_layer_kernel(
        {nSp, spacepoints, measurements, volumeToLayerMap_buf,
         surfaceToLayerMap_buf, layerType_buf, reducedSP_buf, layerCounts_buf,
         spacepointsLayer_buf, cfg.volumeToLayerMap.size(),
         cfg.surfaceToLayerMap.size(),
         cfg.gbts_count_spacepoints_by_layer_params});

    vecmem::vector<unsigned int> layerCounts(cfg.nLayers + 1, mr().host);
    copy()(vecmem::get_data(layerCounts_buf), layerCounts)->wait();

    const unsigned int nNodes =
        static_cast<unsigned int>(layerCounts[cfg.nLayers]);
    TRACCC_DEBUG("nNodes " << nNodes);
    if (nNodes == 0) {
        TRACCC_WARNING("No nodes were found after spacepoint counting");
        return node_making_output{};
    }

    vecmem::data::vector_buffer<float4> sp_params_buf(nSp, mr().main);
    copy().setup(sp_params_buf)->ignore();
    vecmem::data::vector_buffer<unsigned int> original_sp_idx_buf(nSp,
                                                                  mr().main);
    copy().setup(original_sp_idx_buf)->ignore();

    // 1. Fused binning: scatter spacepoints into layer-ordered slots, compute
    //    their eta/phi bin indices and fill the (eta, phi) histogram, all in a
    //    single pass.
    vecmem::data::vector_buffer<std::pair<unsigned int, unsigned int>>
        layer_info_buf(cfg.nLayers, mr().main);
    copy().setup(layer_info_buf)->ignore();
    copy()(vecmem::get_data(cfg.layerInfo.info), layer_info_buf)->ignore();

    vecmem::data::vector_buffer<std::pair<float, float>> layer_geo_buf(
        cfg.nLayers, mr().main);
    copy().setup(layer_geo_buf)->ignore();
    copy()(vecmem::get_data(cfg.layerInfo.geo), layer_geo_buf)->ignore();

    vecmem::data::vector_buffer<unsigned int> node_phi_index_buf(nNodes,
                                                                 mr().main);
    copy().setup(node_phi_index_buf)->ignore();

    vecmem::data::vector_buffer<unsigned int> node_eta_index_buf(nNodes,
                                                                 mr().main);
    copy().setup(node_eta_index_buf)->ignore();

    const unsigned int hist_size = cfg.n_eta_bins * cfg.n_phi_bins;
    vecmem::data::vector_buffer<unsigned int> eta_phi_histo_buf(hist_size,
                                                                mr().main);
    copy().setup(eta_phi_histo_buf)->ignore();
    copy().memset(eta_phi_histo_buf, 0)->ignore();
    vecmem::data::vector_buffer<unsigned int> phi_cusums_buf(hist_size,
                                                             mr().main);
    copy().setup(phi_cusums_buf)->ignore();

    gbts_bin_spacepoints_kernel(
        {nSp, cfg.n_phi_bins, sp_params_buf, reducedSP_buf, layerCounts_buf,
         spacepointsLayer_buf, original_sp_idx_buf, layer_info_buf,
         layer_geo_buf, node_eta_index_buf, node_phi_index_buf,
         eta_phi_histo_buf});

    vecmem::data::vector_buffer<unsigned int> eta_node_counter_buf(
        cfg.n_eta_bins, mr().main);
    copy().setup(eta_node_counter_buf)->ignore();

    gbts_count_eta_phi_bins_kernel({cfg.n_eta_bins, cfg.n_phi_bins,
                                    eta_phi_histo_buf, eta_node_counter_buf,
                                    phi_cusums_buf});

    vecmem::vector<unsigned int> eta_sums(cfg.n_eta_bins, mr().host);
    copy()(vecmem::get_data(eta_node_counter_buf), eta_sums)->wait();

    vecmem::vector<unsigned int> eta_bin_views(2 * cfg.n_eta_bins, mr().host);
    for (unsigned int view_idx = 0; view_idx < cfg.n_eta_bins; view_idx++) {
        const unsigned int pos = 2 * view_idx;
        eta_bin_views[pos] = (view_idx == 0) ? 0 : eta_sums[view_idx - 1];
        eta_bin_views[pos + 1] = eta_sums[view_idx];
    }

    gbts_prefix_sum_eta_phi_bins_kernel(
        {cfg.n_eta_bins, cfg.n_phi_bins, eta_node_counter_buf, phi_cusums_buf});

    vecmem::data::vector_buffer<float4> node_params_buf(nNodes, mr().main);
    copy().setup(node_params_buf)->ignore();
    vecmem::data::vector_buffer<float> node_phi_buf(nNodes, mr().main);
    copy().setup(node_phi_buf)->ignore();
    vecmem::data::vector_buffer<unsigned int> node_index_buf(nNodes, mr().main);
    copy().setup(node_index_buf)->ignore();

    // Optional tau LUT consumed by device::gbts_sort_nodes when
    // cfg.gbts_sort_nodes_params.useTauLUT is set. A size-1 dummy is allocated
    // when the LUT is unused so the kernel always receives a valid (never-read)
    // view.
    const unsigned int tau_lut_size = std::max<unsigned int>(
        1u, static_cast<unsigned int>(cfg.tau_lut.size()));
    vecmem::data::vector_buffer<float> tau_lut_buf(tau_lut_size, mr().main);
    copy().setup(tau_lut_buf)->ignore();
    if (!cfg.tau_lut.empty()) {
        copy()(vecmem::get_data(cfg.tau_lut), tau_lut_buf)->ignore();
    }

    gbts_sort_nodes_kernel({nNodes, cfg.n_phi_bins, sp_params_buf,
                            node_eta_index_buf, node_phi_index_buf,
                            phi_cusums_buf, node_params_buf, node_phi_buf,
                            node_index_buf, original_sp_idx_buf, tau_lut_buf,
                            cfg.gbts_sort_nodes_params});

    vecmem::data::vector_buffer<unsigned int> eta_bin_views_buf(
        2 * cfg.n_eta_bins, mr().main);
    copy().setup(eta_bin_views_buf)->ignore();
    copy()(vecmem::get_data(eta_bin_views), eta_bin_views_buf)->ignore();

    vecmem::data::vector_buffer<float> bin_rads_buf(2 * cfg.n_eta_bins,
                                                    mr().main);
    copy().setup(bin_rads_buf)->ignore();

    gbts_find_minmax_radius_kernel(
        {cfg.n_eta_bins, eta_bin_views_buf, node_params_buf, bin_rads_buf});

    vecmem::vector<float> bin_rads(2 * cfg.n_eta_bins, mr().host);
    copy()(vecmem::get_data(bin_rads_buf), bin_rads)->wait();

    return node_making_output{std::move(reducedSP_buf),
                              std::move(node_params_buf),
                              std::move(node_phi_buf),
                              std::move(node_index_buf),
                              std::move(bin_rads),
                              std::move(eta_bin_views),
                              nNodes};
}

// Stage 2:
// Find edges between compatible nodes
// The main output is a graph in
// the form of an edge list (array of node index pairs)
// and an accompanying array of edge parameters
// (exp(-eta), curvature, extrapolated phi at node1, extrapolated phi at node2).
auto gbts_seeding_algorithm::create_edges(
    vecmem::data::vector_buffer<float4> node_params,
    vecmem::data::vector_buffer<float> node_phi,
    vecmem::data::vector_buffer<unsigned int> node_index,
    const vecmem::vector<float>& bin_rads,
    const vecmem::vector<unsigned int>& eta_bin_views,
    const unsigned int nNodes,
    vecmem::data::vector_buffer<unsigned int>& counters_buf,
    vecmem::vector<unsigned int>& h_counters) const -> graph_making_output {

    const gbts_seedfinder_config& cfg = m_config;
    unsigned int* d_counters = counters_buf.ptr();

    // CPU: build the per-bin-pair work list (begin/end node ranges + phi search
    // window) on the host from the eta-bin views, splitting large bins into
    // node_buffer_length-sized chunks. Two passes: count then fill.
    unsigned int nBinPairs = 0;
    for (const std::pair<unsigned int, unsigned int>& binPair : cfg.binTables) {
        const unsigned int bin1_begin = eta_bin_views[2 * binPair.first];
        const unsigned int bin1_end = eta_bin_views[2 * binPair.first + 1];
        unsigned int nNodesInBin1 = bin1_end - bin1_begin;
        if (bin1_begin > bin1_end) {
            nNodesInBin1 = bin1_begin - bin1_end;
        }
        nBinPairs += 1 + (nNodesInBin1 - 1) / gbts_consts::node_buffer_length;
    }

    vecmem::vector<unsigned int> bin_pair_views(4 * nBinPairs, mr().host);
    vecmem::vector<float> bin_pair_dphi(nBinPairs, mr().host);

    unsigned int pairIdx = 0;
    for (const std::pair<unsigned int, unsigned int>& binPair : cfg.binTables) {
        const float rb1 = bin_rads[2 * binPair.first];

        const unsigned int begin_bin1 = eta_bin_views[2 * binPair.first];
        const unsigned int end_bin1 = eta_bin_views[2 * binPair.first + 1];
        if (begin_bin1 == end_bin1) {
            continue;
        }
        if (eta_bin_views[2 * binPair.second] ==
            eta_bin_views[2 * binPair.second + 1]) {
            continue;
        }

        const float rb2 = bin_rads[2 * binPair.second + 1];
        const float maxDeltaR = std::fabs(rb2 - rb1);

        float deltaPhi = cfg.gbts_dphi_window_params.min_delta_phi +
                         cfg.gbts_dphi_window_params.dphi_coeff * maxDeltaR;
        if (maxDeltaR < cfg.gbts_dphi_window_params.low_dr_threshold) {
            deltaPhi =
                cfg.gbts_dphi_window_params.min_delta_phi_low_dr +
                cfg.gbts_dphi_window_params.dphi_coeff_low_dr * maxDeltaR;
        }

        unsigned int currBegin_bin1 = begin_bin1;
        unsigned int currEnd_bin1 =
            end_bin1 < gbts_consts::node_buffer_length
                ? end_bin1
                : begin_bin1 + gbts_consts::node_buffer_length;

        for (; currEnd_bin1 < end_bin1;
             currEnd_bin1 += gbts_consts::node_buffer_length, pairIdx++) {
            const unsigned int offset = 4 * pairIdx;
            bin_pair_views[offset] = currBegin_bin1;
            bin_pair_views[1 + offset] = currEnd_bin1;
            bin_pair_views[2 + offset] = eta_bin_views[2 * binPair.second];
            bin_pair_views[3 + offset] = eta_bin_views[2 * binPair.second + 1];
            bin_pair_dphi[pairIdx] = deltaPhi;
            currBegin_bin1 = currEnd_bin1;
        }
        currEnd_bin1 = end_bin1;

        const unsigned int offset = 4 * pairIdx;
        bin_pair_views[offset] = currBegin_bin1;
        bin_pair_views[1 + offset] = currEnd_bin1;
        bin_pair_views[2 + offset] = eta_bin_views[2 * binPair.second];
        bin_pair_views[3 + offset] = eta_bin_views[2 * binPair.second + 1];
        bin_pair_dphi[pairIdx] = deltaPhi;
        pairIdx++;
    }
    const unsigned int nUsedBinPairs = pairIdx;
    TRACCC_DEBUG("nUsedBinPairs " << nUsedBinPairs);
    if (nUsedBinPairs == 0) {
        TRACCC_WARNING("No bin pairs were used for edge finding");
        return graph_making_output{};
    }

    vecmem::data::vector_buffer<unsigned int> bin_pair_views_buf(
        4 * nUsedBinPairs, mr().main);
    copy().setup(bin_pair_views_buf)->ignore();
    copy()(vecmem::get_data(bin_pair_views), bin_pair_views_buf)->ignore();

    vecmem::data::vector_buffer<float> bin_pair_dphi_buf(nUsedBinPairs,
                                                         mr().main);
    copy().setup(bin_pair_dphi_buf)->ignore();
    copy()(vecmem::get_data(bin_pair_dphi), bin_pair_dphi_buf)->ignore();

    // 2. Find edges between spacepoint pairs.
    const unsigned int nMaxEdges = cfg.max_edges_factor * nNodes;
    // Packed per-edge parameter buffer ([exp(-eta), curv, phi_z, phi_w]).
    vecmem::data::vector_buffer<float4> edge_params_buf(nMaxEdges, mr().main);
    copy().setup(edge_params_buf)->ignore();
    vecmem::data::vector_buffer<uint2> edge_nodes_buf(nMaxEdges, mr().main);
    copy().setup(edge_nodes_buf)->ignore();
    vecmem::data::vector_buffer<unsigned int> num_incoming_edges_buf(nNodes + 1,
                                                                     mr().main);
    copy().setup(num_incoming_edges_buf)->ignore();
    copy().memset(num_incoming_edges_buf, 0)->ignore();

    gbts_make_graph_edges_kernel(
        {nUsedBinPairs, nMaxEdges, cfg.n_phi_bins, bin_pair_views_buf,
         bin_pair_dphi_buf, node_params, node_phi,
         cfg.gbts_make_graph_edges_params, d_counters + gbts_counter::nEdges,
         edge_nodes_buf, edge_params_buf, num_incoming_edges_buf});

    // Read back the number of edges produced.
    copy()(counters_buf, h_counters)->wait();

    unsigned int nEdges = h_counters[gbts_counter::nEdges];
    TRACCC_DEBUG("Created " << nEdges << " edges with a cap of " << nMaxEdges);
    if (nEdges > nMaxEdges) {
        TRACCC_WARNING("Number of edges exceeds the maximum allowed, Removing "
                       << nEdges - nMaxEdges << " edges");
        nEdges = nMaxEdges;
    } else if (nEdges == 0) {
        TRACCC_WARNING("No edges were found");
        return graph_making_output{};
    }

    // 3. Link edges and nodes.
    vecmem::data::vector_buffer<unsigned int> edge_links_buf(nEdges, mr().main);
    copy().setup(edge_links_buf)->ignore();

    gbts_link_graph_edges_kernel(
        {nEdges, edge_nodes_buf, edge_links_buf, num_incoming_edges_buf});

    // 4. Edge matching to create edge-to-edge connections.
    vecmem::data::vector_buffer<unsigned char> num_neighbours_buf(nEdges,
                                                                  mr().main);
    copy().setup(num_neighbours_buf)->ignore();
    copy().memset(num_neighbours_buf, 0)->ignore();

    vecmem::data::vector_buffer<int> reIndexer_buf(nEdges, mr().main);
    copy().setup(reIndexer_buf)->ignore();
    // Byte-fill 0xFF -> int -1, the "edge not kept" sentinel checked by
    // gbts_reindex_edges / gbts_compress_graph.
    copy().memset(reIndexer_buf, 0xFF)->ignore();

    vecmem::data::vector_buffer<unsigned int> neighbours_buf(
        cfg.max_num_neighbours * nEdges, mr().main);
    copy().setup(neighbours_buf)->ignore();
    copy().memset(neighbours_buf, 0)->ignore();

    gbts_match_graph_edges_kernel(
        {nEdges, cfg.max_num_neighbours, cfg.gbts_match_graph_edges_params,
         edge_params_buf, edge_nodes_buf, num_incoming_edges_buf,
         edge_links_buf, num_neighbours_buf, neighbours_buf, reIndexer_buf,
         d_counters + gbts_counter::nConnections});

    // 5. Edge re-indexing to keep only edges involved in any connection.
    gbts_reindex_edges_kernel(
        {nEdges, reIndexer_buf, d_counters + gbts_counter::nConnectedEdges});

    copy()(counters_buf, h_counters)->wait();

    const unsigned int nConnections = h_counters[gbts_counter::nConnections];
    const unsigned int nConnectedEdges =
        h_counters[gbts_counter::nConnectedEdges];
    TRACCC_DEBUG("created " << nConnections << " edge links, found "
                            << nConnectedEdges
                            << " connected edges for seed extraction");
    if (nConnectedEdges == 0) {
        TRACCC_WARNING("No connected edges were found");
        return graph_making_output{};
    }

    const unsigned int nIntsPerEdge = 2 + 1 + cfg.max_num_neighbours;
    vecmem::data::vector_buffer<unsigned int> output_graph_buf(
        nConnectedEdges * nIntsPerEdge, mr().main);
    copy().setup(output_graph_buf)->ignore();

    gbts_compress_graph_kernel(
        {nEdges, cfg.max_num_neighbours, node_index, edge_nodes_buf,
         num_neighbours_buf, neighbours_buf, reIndexer_buf, output_graph_buf});

    return graph_making_output{std::move(output_graph_buf), nConnectedEdges};
}

// Stage 3:
// Find seed candidates as long chains of connected edges using a CCA
// Then fit the potential seeds (eta, phi, curvature).
// Finally, disambiguate them by repeated seed-vs-edge bidding rounds.
auto gbts_seeding_algorithm::extract_seeds(
    vecmem::data::vector_buffer<unsigned int>& output_graph,
    vecmem::data::vector_buffer<float4>& reducedSP,
    const unsigned int nConnectedEdges, const unsigned int nSp,
    vecmem::data::vector_buffer<unsigned int>& counters_buf,
    vecmem::vector<unsigned int>& h_counters) const
    -> edm::seed_collection::buffer {

    const gbts_seedfinder_config& cfg = m_config;
    unsigned int* d_counters = counters_buf.ptr();

    // 6. Find longest segments with CCA.
    // active_edges is the per-edge "next iter index" flag: it holds `iter`
    // while the edge is active in iteration `iter`, and -1 once it settles.
    // Iteration 0 writes every entry before any later iteration reads it, so
    // no initialisation is required.
    vecmem::data::vector_buffer<char> active_edges_buf(nConnectedEdges,
                                                       mr().main);
    copy().setup(active_edges_buf)->ignore();

    vecmem::data::vector_buffer<unsigned char> levels_buf(2 * nConnectedEdges,
                                                          mr().main);
    copy().setup(levels_buf)->ignore();
    // Initialise to 1 so a level counts the maximum number of edge segments
    // for a seed originating at the edge.
    copy().memset(levels_buf, 0x1)->ignore();

    vecmem::data::vector_buffer<short2> outgoing_paths_buf(nConnectedEdges,
                                                           mr().main);
    copy().setup(outgoing_paths_buf)->ignore();

    for (unsigned char iter = 0;
         iter < traccc::device::gbts_consts::max_cca_iter; ++iter) {
        gbts_run_cca_iteration_kernel({nConnectedEdges, cfg.max_num_neighbours,
                                       cfg.minLevel, output_graph, levels_buf,
                                       active_edges_buf, outgoing_paths_buf,
                                       iter});
    }

    gbts_count_terminus_edges_kernel(
        {nConnectedEdges, outgoing_paths_buf, d_counters + gbts_counter::nPaths,
         d_counters + gbts_counter::nTerminusEdges});

    copy()(counters_buf, h_counters)->wait();

    const unsigned int nPaths = h_counters[gbts_counter::nPaths];
    const unsigned int nTerminusEdges =
        h_counters[gbts_counter::nTerminusEdges];
    if (nTerminusEdges == 0) {
        TRACCC_WARNING("No terminus edges were found");
        return {0, mr().main};
    }

    TRACCC_DEBUG(nPaths << " size of path store | nTerminusEdges "
                        << nTerminusEdges);

    vecmem::data::vector_buffer<int2> path_store_buf(nPaths + nTerminusEdges,
                                                     mr().main);
    copy().setup(path_store_buf)->ignore();
    vecmem::data::vector_buffer<int2> seed_proposals_buf(nPaths, mr().main);
    copy().setup(seed_proposals_buf)->ignore();
    vecmem::data::vector_buffer<char> seed_ambiguity_buf(nPaths, mr().main);
    copy().setup(seed_ambiguity_buf)->ignore();

    vecmem::data::vector_buffer<unsigned long long int> edge_bids_buf(
        nConnectedEdges, mr().main);
    copy().setup(edge_bids_buf)->ignore();
    copy().memset(edge_bids_buf, 0)->ignore();

    gbts_add_terminus_to_path_store_kernel(
        {nConnectedEdges, path_store_buf, outgoing_paths_buf});

    gbts_fill_path_store_kernel({nTerminusEdges, cfg.max_num_neighbours, nPaths,
                                 path_store_buf, output_graph, levels_buf,
                                 d_counters + gbts_counter::nTerminusEdges});

    gbts_fit_segments_kernel(
        {nPaths, nTerminusEdges, cfg.max_num_neighbours, cfg.minLevel,
         reducedSP, output_graph, path_store_buf, seed_proposals_buf,
         edge_bids_buf, seed_ambiguity_buf,
         d_counters + gbts_counter::nTerminusEdges,
         d_counters + gbts_counter::nProps, cfg.gbts_fit_segments_params,
         cfg.gbts_make_graph_edges_params.max_z0});

    copy()(counters_buf, h_counters)->wait();

    const unsigned int nProps = h_counters[gbts_counter::nProps];
    TRACCC_DEBUG("nProps " << nProps);
    if (nProps == 0) {
        TRACCC_WARNING("No seed proposals were found");
        return {0, mr().main};
    }

    // 7. Disambiguate seeds through repeated seed-vs-edge bidding rounds.
    for (unsigned int round = 0; round < cfg.edge_bidding_rounds; ++round) {
        copy().memset(edge_bids_buf, 0)->ignore();

        gbts_rebid_seeds_for_edges_kernel(
            {nProps, path_store_buf, seed_proposals_buf, edge_bids_buf,
             seed_ambiguity_buf, d_counters + gbts_counter::nRejected,
             round == 0u});

        gbts_reset_edge_bids_kernel({nProps, path_store_buf, seed_proposals_buf,
                                     edge_bids_buf, seed_ambiguity_buf,
                                     d_counters + gbts_counter::nRejected});
    }

    copy()(counters_buf, h_counters)->wait();
    const unsigned int nRejectedProps = h_counters[gbts_counter::nRejected];
    const unsigned int nSeeds =
        (nRejectedProps >= nProps) ? 0u : nProps - nRejectedProps;

    TRACCC_DEBUG("Rejected " << nRejectedProps << " out of " << nProps
                             << " seed proposals");
    if (nSeeds == 0) {
        TRACCC_WARNING("All seed proposals were rejected");
        return {0, mr().main};
    }

    // 8. Convert to 3sp seeds and make output buffer.
    edm::seed_collection::buffer output_seeds(
        2 * nSeeds, mr().main, vecmem::data::buffer_type::resizable);
    copy().setup(output_seeds)->ignore();

    vecmem::data::vector_buffer<unsigned long long int> hit_bids_buf(nSp,
                                                                     mr().main);
    copy().setup(hit_bids_buf)->ignore();
    copy().memset(hit_bids_buf, 0)->ignore();

    const unsigned int edge_size = 1u + 2u + cfg.max_num_neighbours;
    gbts_bid_seeds_for_hits_kernel({nProps, nSeeds, edge_size, output_graph,
                                    seed_proposals_buf, path_store_buf,
                                    seed_ambiguity_buf, hit_bids_buf});

    gbts_convert_seeds_kernel(
        {nProps, nSeeds, cfg.max_num_neighbours, seed_proposals_buf,
         seed_ambiguity_buf, path_store_buf, output_graph, reducedSP,
         output_seeds, hit_bids_buf, cfg.gbts_convert_seeds_params});

    const unsigned int outputSeeds = copy().get_size(output_seeds);
    TRACCC_DEBUG("GBTS found " << outputSeeds << " seeds");
    return output_seeds;
}

gbts_seeding_algorithm::gbts_seeding_algorithm(
    const gbts_seedfinder_config& cfg, const memory_resource& mr,
    const vecmem::copy& copy, std::unique_ptr<const Logger> logger)
    : messaging(std::move(logger)), algorithm_base{mr, copy}, m_config{cfg} {}

auto gbts_seeding_algorithm::operator()(
    const edm::spacepoint_collection::const_view& spacepoints,
    const edm::measurement_collection::const_view& measurements) const
    -> output_type {

    const unsigned int nSp = copy().get_size(spacepoints);
    TRACCC_DEBUG("nSp " << nSp);
    if (nSp == 0) {
        TRACCC_WARNING("No spacepoints were found in the event");
        return {0, mr().main};
    }

    // Stage 1: bin spacepoints and create nodes with the parameters (eta, phi,
    // r, z).
    node_making_output nodes = make_nodes(spacepoints, measurements);
    if (nodes.nNodes == 0) {
        // No nodes survived spacepoint counting -> no seeds.
        return {0, mr().main};
    }

    // Named counters shared by the graph-making and seed-extraction stages.
    vecmem::data::vector_buffer<unsigned int> counters_buf(
        gbts_counter::nCounters, mr().main);
    copy().setup(counters_buf)->ignore();
    copy().memset(counters_buf, 0)->ignore();
    vecmem::vector<unsigned int> h_counters(
        gbts_counter::nCounters, mr().host ? mr().host : &(mr().main));

    // Stage 2: graph. The per-node buffers are moved in so they are released
    // when create_gbts_edges_from_nodes returns, along with all the edge/link
    // transients.
    graph_making_output graph = create_edges(
        std::move(nodes.node_params), std::move(nodes.node_phi),
        std::move(nodes.node_index), nodes.bin_rads, nodes.eta_bin_views,
        nodes.nNodes, counters_buf, h_counters);
    if (graph.nConnectedEdges == 0) {
        // No connected edges survived graph making -> no seeds.
        return {0, mr().main};
    }

    // Stage 3: Create seeds from the graph edges.
    return extract_seeds(graph.output_graph, nodes.reducedSP,
                         graph.nConnectedEdges, nSp, counters_buf, h_counters);
}

}  // namespace traccc::device
