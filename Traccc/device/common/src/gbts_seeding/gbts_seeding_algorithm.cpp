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
#include <vector>

namespace traccc::device {

// Stage 1:
// Bin the spacepoints by layer in eta and phi,
// compute the node parameters (x, y, z, w)
// and the bin-wise min/max radius for the graph-building cuts.
auto gbts_seeding_algorithm::make_nodes(
    const edm::spacepoint_collection::const_view& spacepoints,
    const edm::measurement_collection::const_view& measurements,
    const unsigned int nSp) const -> node_making_output {
  const gbts_seedfinder_config& cfg = m_config;

  // 1. Fused binning: assign each spacepoint a layer (or reject it), write
  //    its reduced parameters, count its eta bin and append its node sort
  //    key, all in a single pass. The key array is sized for the worst case
  //    (every spacepoint accepted); only the first nNodes slots get used.

  vecmem::data::vector_buffer<float4> reducedSP_buf(nSp, mr().main);
  copy().setup(reducedSP_buf)->ignore();

  vecmem::data::vector_buffer<unsigned long long int> sort_keys_buf(nSp,
                                                                    mr().main);
  copy().setup(sort_keys_buf)->ignore();
  vecmem::data::vector_buffer<unsigned int> sort_values_buf(nSp, mr().main);
  copy().setup(sort_values_buf)->ignore();

  // Per eta bin node counts, scanned in place by the kernel launcher into
  // the node offsets.
  vecmem::data::vector_buffer<unsigned int> eta_bin_offsets_buf(
      cfg.n_eta_bins + 1, mr().main);
  copy().setup(eta_bin_offsets_buf)->ignore();
  copy().memset(eta_bin_offsets_buf, 0)->ignore();

  gbts_bin_spacepoints_kernel(
      {nSp, cfg.n_eta_bins, spacepoints, measurements,
       m_volume_to_layer_map_buffer, m_surface_to_layer_map_buffer,
       m_layer_type_buffer, m_layer_info_buffer, m_layer_geo_buffer,
       reducedSP_buf, eta_bin_offsets_buf, sort_keys_buf, sort_values_buf,
       cfg.volumeToLayerMap.size(), cfg.surfaceToLayerMap.size(),
       cfg.gbts_count_spacepoints_by_layer_params});

  // The node count sizes the node buffers: one readback of the last offset.
  vecmem::vector<unsigned int> h_nNodes(1,
                                        mr().host ? mr().host : &(mr().main));
  copy()(vecmem::data::vector_view<unsigned int>{1u, eta_bin_offsets_buf.ptr() +
                                                         cfg.n_eta_bins},
         vecmem::get_data(h_nNodes))
      ->wait();
  const unsigned int nNodes = h_nNodes[0];
  TRACCC_DEBUG("nNodes " << nNodes);
  if (nNodes == 0) {
    TRACCC_WARNING("No nodes were found after spacepoint binning");
    return node_making_output{};
  }

  vecmem::data::vector_buffer<float4> node_params_buf(nNodes, mr().main);
  copy().setup(node_params_buf)->ignore();
  vecmem::data::vector_buffer<float> node_phi_buf(nNodes, mr().main);
  copy().setup(node_phi_buf)->ignore();
  vecmem::data::vector_buffer<unsigned int> node_index_buf(nNodes, mr().main);
  copy().setup(node_index_buf)->ignore();

  gbts_sort_nodes_kernel({nNodes, reducedSP_buf, sort_keys_buf, sort_values_buf,
                          node_params_buf, node_phi_buf, node_index_buf,
                          m_tau_lut_buffer, cfg.gbts_sort_nodes_params});

  vecmem::data::vector_buffer<float> bin_rads_buf(2 * cfg.n_eta_bins,
                                                  mr().main);
  copy().setup(bin_rads_buf)->ignore();

  gbts_find_minmax_radius_kernel(
      {cfg.n_eta_bins, eta_bin_offsets_buf, node_params_buf, bin_rads_buf});

  return node_making_output{std::move(reducedSP_buf),
                            std::move(node_params_buf),
                            std::move(node_phi_buf),
                            std::move(node_index_buf),
                            std::move(bin_rads_buf),
                            std::move(eta_bin_offsets_buf),
                            nNodes};
}

// Stage 2:
// Find edges between compatible nodes. The output is the compacted graph
// (an edge list with the accepted neighbour edges of every edge).
auto gbts_seeding_algorithm::create_edges(
    vecmem::data::vector_buffer<float4> node_params,
    vecmem::data::vector_buffer<float> node_phi,
    vecmem::data::vector_buffer<unsigned int> node_index,
    vecmem::data::vector_buffer<float> bin_rads,
    vecmem::data::vector_buffer<unsigned int> eta_bin_offsets,
    const unsigned int nNodes, const unsigned int nSp,
    vecmem::data::vector_buffer<unsigned int>& counters_buf,
    vecmem::vector<unsigned int>& h_counters) const -> graph_making_output {
  const gbts_seedfinder_config& cfg = m_config;

  // Per bin pair tables.
  vecmem::data::vector_buffer<uint2> bin_pairs_buf(m_nBinPairs, mr().main);
  copy().setup(bin_pairs_buf)->ignore();
  copy()(vecmem::get_data(m_bin_pairs), bin_pairs_buf)->ignore();
  vecmem::data::vector_buffer<unsigned int> pair_group_begin_buf(m_nBinPairs,
                                                                 mr().main);
  copy().setup(pair_group_begin_buf)->ignore();
  copy()(vecmem::get_data(m_pair_group_begin), pair_group_begin_buf)->ignore();

  // 1. Work list: one item per (bin pair, chunk of the inner bin).
  constexpr unsigned int chunk_size = gbts_consts::edge_chunk_size;
  const unsigned int nWorkMax =
      m_nBinPairs + m_maxPairsPerBin1 * (nNodes / chunk_size + 1u);
  vecmem::data::vector_buffer<unsigned int> pair_work_begin_buf(m_nBinPairs + 1,
                                                                mr().main);
  copy().setup(pair_work_begin_buf)->ignore();
  vecmem::data::vector_buffer<gbts_graph_edges_work_item> work_items_buf(
      nWorkMax, mr().main);
  copy().setup(work_items_buf)->ignore();
  gbts_build_edge_work_list_kernel({m_nBinPairs, bin_pairs_buf, eta_bin_offsets,
                                    bin_rads, cfg.gbts_dphi_window_params,
                                    pair_work_begin_buf, work_items_buf});

  // 2. Count the edges per inner node.
  const float max_Kappa =
      std::max(cfg.gbts_make_graph_edges_params.max_Kappa_low_tau,
               cfg.gbts_make_graph_edges_params.max_Kappa_high_tau);
  const edge_params_converter edge_params_maker(
      max_Kappa, cfg.gbts_sort_nodes_params.maxTau);
  vecmem::data::vector_buffer<unsigned int> edge_counts_buf(
      nWorkMax * chunk_size, mr().main);
  copy().setup(edge_counts_buf)->ignore();
  // Per node edge counts at [node + 1]; the scan turns them into the edge
  // buckets and leaves the edge count at [nNodes].
  vecmem::data::vector_buffer<unsigned int> num_outgoing_edges_buf(nNodes + 1,
                                                                   mr().main);
  copy().setup(num_outgoing_edges_buf)->ignore();
  copy().memset(num_outgoing_edges_buf, 0)->ignore();
  gbts_count_graph_edges_kernel({nWorkMax, pair_work_begin_buf, work_items_buf,
                                 node_params, node_phi,
                                 cfg.gbts_make_graph_edges_params,
                                 edge_counts_buf, num_outgoing_edges_buf});

  // 3. Write the edges into their slots.
  const unsigned int nEdgesMax = cfg.max_edges_per_spacepoint * nSp;
  vecmem::data::vector_buffer<uint2> edge_nodes_buf(nEdgesMax, mr().main);
  copy().setup(edge_nodes_buf)->ignore();
  vecmem::data::vector_buffer<short4> edge_params_buf(nEdgesMax, mr().main);
  copy().setup(edge_params_buf)->ignore();

  gbts_fill_graph_edges_kernel(
      {nWorkMax, pair_work_begin_buf, work_items_buf, node_params, node_phi,
       cfg.gbts_make_graph_edges_params, pair_group_begin_buf, edge_counts_buf,
       num_outgoing_edges_buf, edge_params_maker, nEdgesMax, edge_nodes_buf,
       edge_params_buf});

  // 4. Match every edge with the edges entering its outer node, then
  //    compact the kept edges.
  vecmem::data::vector_buffer<unsigned char> num_neighbours_buf(nEdgesMax,
                                                                mr().main);
  copy().setup(num_neighbours_buf)->ignore();
  vecmem::data::vector_buffer<unsigned int> neighbours_buf(
      cfg.max_num_neighbours * nEdgesMax, mr().main);
  copy().setup(neighbours_buf)->ignore();
  vecmem::data::vector_buffer<unsigned int> edge_kept_buf(nEdgesMax, mr().main);
  copy().setup(edge_kept_buf)->ignore();
  copy().memset(edge_kept_buf, 0)->ignore();
  vecmem::data::vector_buffer<unsigned int> reIndexer_buf(nEdgesMax, mr().main);
  copy().setup(reIndexer_buf)->ignore();

  gbts_match_graph_edges_kernel(
      {nEdgesMax, cfg.max_num_neighbours, cfg.gbts_match_graph_edges_params,
       edge_params_buf, edge_nodes_buf, num_outgoing_edges_buf,
       num_neighbours_buf, neighbours_buf, edge_kept_buf, reIndexer_buf,
       edge_params_maker});

  // 5. Compress the kept edges into the output graph.
  const unsigned int nConnectedEdgesMax =
      cfg.max_connected_edges_per_spacepoint * nSp;
  const unsigned int nIntsPerEdge =
      gbts_consts::nei_start + cfg.max_num_neighbours;
  vecmem::data::vector_buffer<unsigned int> output_graph_buf(
      nConnectedEdgesMax * nIntsPerEdge, mr().main);
  copy().setup(output_graph_buf)->ignore();

  gbts_compress_graph_kernel({nEdgesMax, num_outgoing_edges_buf,
                              nConnectedEdgesMax, cfg.max_num_neighbours,
                              node_index, edge_nodes_buf, num_neighbours_buf,
                              neighbours_buf, reIndexer_buf, output_graph_buf});

  // The seed extraction needs the kept-edge count on the host.
  copy()(
      vecmem::data::vector_view<unsigned int>{
          1u, num_outgoing_edges_buf.ptr() + nNodes},
      vecmem::data::vector_view<unsigned int>{
          1u, counters_buf.ptr() + gbts_counter::nEdgesTotal})
      ->ignore();
  copy()(vecmem::data::vector_view<unsigned int>{1u, reIndexer_buf.ptr() +
                                                         (nEdgesMax - 1u)},
         vecmem::data::vector_view<unsigned int>{
             1u, counters_buf.ptr() + gbts_counter::nConnectedEdges})
      ->ignore();
  copy()(counters_buf, h_counters)->wait();
  if (h_counters[gbts_counter::nEdgesTotal] > nEdgesMax) {
    TRACCC_WARNING("Edge buffer capacity ("
                   << nEdgesMax << ") exceeded, "
                   << h_counters[gbts_counter::nEdgesTotal] - nEdgesMax
                   << " edges were dropped. "
                      "Increase max_edges_per_spacepoint");
  }
  unsigned int nConnectedEdges = h_counters[gbts_counter::nConnectedEdges];
  if (nConnectedEdges > nConnectedEdgesMax) {
    TRACCC_WARNING("Compacted graph capacity ("
                   << nConnectedEdgesMax << ") exceeded, "
                   << nConnectedEdges - nConnectedEdgesMax
                   << " connected edges were dropped. "
                      "Increase max_connected_edges_per_spacepoint");
    nConnectedEdges = nConnectedEdgesMax;
  }
  TRACCC_DEBUG("found " << nConnectedEdges
                        << " connected edges for seed extraction");
  if (nConnectedEdges == 0) {
    TRACCC_WARNING("No connected edges were found");
    return graph_making_output{};
  }
  return graph_making_output{std::move(output_graph_buf), nConnectedEdges};
}

// Stage 3:
// Find seed candidates as long chains of connected edges using a CCA
// Then fit the potential seeds (eta, phi, curvature).
// Finally, disambiguate them by seed-vs-edge and seed-vs-hit bidding.
auto gbts_seeding_algorithm::extract_seeds(
    vecmem::data::vector_buffer<unsigned int>& output_graph,
    vecmem::data::vector_buffer<float4>& reducedSP,
    const unsigned int nConnectedEdges, const unsigned int nSp,
    vecmem::vector<unsigned int>& h_counters) const
    -> edm::seed_collection::buffer {
  const gbts_seedfinder_config& cfg = m_config;

  // 6. Find longest segments with CCA.
  vecmem::data::vector_buffer<unsigned char> levels_buf(nConnectedEdges,
                                                        mr().main);
  copy().setup(levels_buf)->ignore();
  // Initialise to 1 so a level counts the maximum number of edge segments
  // for a seed originating at the edge.
  copy().memset(levels_buf, 1)->ignore();
  vecmem::data::vector_buffer<unsigned char> has_parent_buf(nConnectedEdges,
                                                            mr().main);
  copy().setup(has_parent_buf)->ignore();
  copy().memset(has_parent_buf, 0)->ignore();
  vecmem::data::vector_buffer<int2> outgoing_paths_buf(nConnectedEdges,
                                                       mr().main);
  copy().setup(outgoing_paths_buf)->ignore();
  vecmem::data::vector_buffer<unsigned int> cca_changed_buf(
      gbts_run_cca_max_sweeps, mr().main);
  copy().setup(cca_changed_buf)->ignore();
  copy().memset(cca_changed_buf, 0)->ignore();

  gbts_run_cca_iteration_payload cca{nConnectedEdges,    cfg.max_num_neighbours,
                                     output_graph,       levels_buf,
                                     outgoing_paths_buf, 0u,
                                     cca_changed_buf};
  for (unsigned char iter = 0; iter < gbts_run_cca_max_sweeps; ++iter) {
    cca.iter = iter;
    gbts_run_cca_iteration_kernel(cca);
  }
  gbts_finish_cca_kernel({nConnectedEdges, cfg.max_num_neighbours, cfg.minLevel,
                          output_graph, levels_buf, outgoing_paths_buf,
                          has_parent_buf});

  // 7. Lay out the path store, the paths of a terminus edge contiguous.
  const unsigned int nPathsMax = nSp;
  const unsigned int nPathsGrid = nSp / 2u;
  vecmem::data::vector_buffer<unsigned int> path_counts_buf(nConnectedEdges,
                                                            mr().main);
  copy().setup(path_counts_buf)->ignore();
  // The last path offset is the total number of paths.
  vecmem::data::vector_view<unsigned int> path_count{
      1u, path_counts_buf.ptr() + (nConnectedEdges - 1u)};
  vecmem::data::vector_buffer<unsigned long long int> edge_bids_buf(
      nConnectedEdges, mr().main);
  copy().setup(edge_bids_buf)->ignore();
  copy().memset(edge_bids_buf, 0)->ignore();
  vecmem::data::vector_buffer<unsigned long long int> hit_bids_buf(nSp,
                                                                   mr().main);
  copy().setup(hit_bids_buf)->ignore();
  copy().memset(hit_bids_buf, 0)->ignore();
  vecmem::data::vector_buffer<char> seed_ambiguity_buf(nPathsMax, mr().main);
  copy().setup(seed_ambiguity_buf)->ignore();
  copy().memset(seed_ambiguity_buf, 0)->ignore();

  gbts_count_paths_kernel(
      {nConnectedEdges, outgoing_paths_buf, has_parent_buf, path_counts_buf});

  // 8. Fill the path store, fit every path and bid for its edge.
  vecmem::data::vector_buffer<int2> path_store_buf(nPathsMax, mr().main);
  copy().setup(path_store_buf)->ignore();
  vecmem::data::vector_buffer<int2> seed_proposals_buf(nPathsMax, mr().main);
  copy().setup(seed_proposals_buf)->ignore();

  gbts_fill_path_store_kernel(
      {nPathsMax, nPathsGrid, path_count, nConnectedEdges,
       cfg.max_num_neighbours, path_store_buf, output_graph, levels_buf,
       outgoing_paths_buf, path_counts_buf, seed_proposals_buf,
       seed_ambiguity_buf, cfg.minLevel, reducedSP,
       cfg.gbts_fit_segments_params, cfg.gbts_make_graph_edges_params.max_z0,
       edge_bids_buf});

  // 9. Bid for the hits and convert the winners to 3sp seeds.
  edm::seed_collection::buffer output_seeds(
      2 * nPathsMax, mr().main, vecmem::data::buffer_type::resizable);
  copy().setup(output_seeds)->ignore();

  const unsigned int edge_size =
      gbts_consts::nei_start + cfg.max_num_neighbours;
  gbts_bid_seeds_for_hits_kernel(
      {nPathsMax, nPathsGrid, path_count, edge_size, output_graph,
       seed_proposals_buf, path_store_buf, seed_ambiguity_buf, hit_bids_buf});

  gbts_convert_seeds_kernel(
      {nPathsMax, nPathsGrid, path_count, cfg.max_num_neighbours,
       seed_proposals_buf, seed_ambiguity_buf, path_store_buf, output_graph,
       reducedSP, output_seeds, hit_bids_buf, cfg.gbts_convert_seeds_params});

  // The caller synchronises for the seed count right after this.
  copy()(path_count,
         vecmem::data::vector_view<unsigned int>{
             1u, h_counters.data() + gbts_counter::nPaths})
      ->wait();
  if (h_counters[gbts_counter::nPaths] > nPathsMax) {
    TRACCC_WARNING("Path store capacity ("
                   << nPathsMax << ") exceeded, "
                   << h_counters[gbts_counter::nPaths] - nPathsMax
                   << " paths were dropped");
  }
  return output_seeds;
}

gbts_seeding_algorithm::gbts_seeding_algorithm(
    const gbts_seedfinder_config& cfg, const memory_resource& mr,
    const vecmem::copy& copy, std::unique_ptr<const Logger> logger)
    : messaging(std::move(logger)),
      algorithm_base{mr, copy},
      m_config{cfg},
      m_volume_to_layer_map_buffer{
          static_cast<unsigned int>(cfg.volumeToLayerMap.size()), mr.main},
      m_layer_type_buffer{cfg.nLayers, mr.main},
      m_layer_info_buffer{cfg.nLayers, mr.main},
      m_layer_geo_buffer{cfg.nLayers, mr.main},
      m_tau_lut_buffer{std::max<unsigned int>(
                           1u, static_cast<unsigned int>(cfg.tau_lut.size())),
                       mr.main} {
  // The copies below may be asynchronous, so they read from m_config (which
  // lives as long as the buffers) rather than from the cfg argument.
  copy.setup(m_volume_to_layer_map_buffer)->ignore();
  copy(vecmem::get_data(m_config.volumeToLayerMap),
       m_volume_to_layer_map_buffer)
      ->ignore();
  if (!m_config.surfaceToLayerMap.empty()) {
    m_surface_to_layer_map_buffer =
        vecmem::data::vector_buffer<std::pair<unsigned int, unsigned int>>(
            static_cast<unsigned int>(m_config.surfaceToLayerMap.size()),
            mr.main);
    copy.setup(m_surface_to_layer_map_buffer)->ignore();
    copy(vecmem::get_data(m_config.surfaceToLayerMap),
         m_surface_to_layer_map_buffer)
        ->ignore();
  }
  copy.setup(m_layer_type_buffer)->ignore();
  copy(vecmem::get_data(m_config.layerInfo.type), m_layer_type_buffer)
      ->ignore();
  copy.setup(m_layer_info_buffer)->ignore();
  copy(vecmem::get_data(m_config.layerInfo.info), m_layer_info_buffer)
      ->ignore();
  copy.setup(m_layer_geo_buffer)->ignore();
  copy(vecmem::get_data(m_config.layerInfo.geo), m_layer_geo_buffer)->ignore();
  // Optional tau LUT consumed by device::gbts_sort_nodes when
  // cfg.gbts_sort_nodes_params.useTauLUT is set. A size-1 dummy is allocated
  // when the LUT is unused so the kernel always receives a valid (never-read)
  // view.
  copy.setup(m_tau_lut_buffer)->ignore();
  if (!m_config.tau_lut.empty()) {
    copy(vecmem::get_data(m_config.tau_lut), m_tau_lut_buffer)->ignore();
  }

  prepare_bin_pairs();
}

void gbts_seeding_algorithm::prepare_bin_pairs() {
  // The fill kernel relies on the bin pairs being sorted by (bin1, bin2)
  // without duplicates.
  std::vector<std::pair<unsigned int, unsigned int>>& binTables =
      m_config.binTables;
  const std::size_t nInput = binTables.size();
  std::erase_if(binTables,
                [this](const std::pair<unsigned int, unsigned int>& p) {
                  return (p.first >= m_config.n_eta_bins) ||
                         (p.second >= m_config.n_eta_bins);
                });
  if (binTables.size() != nInput) {
    TRACCC_ERROR("Dropped " << nInput - binTables.size()
                            << " bin pairs referring to eta bins >= "
                            << m_config.n_eta_bins);
  }
  std::ranges::sort(binTables);
  const auto duplicates = std::ranges::unique(binTables);
  if (!duplicates.empty()) {
    TRACCC_WARNING("Removed " << duplicates.size()
                              << " duplicate bin pairs from binTables");
    binTables.erase(duplicates.begin(), duplicates.end());
  }
  m_nBinPairs = static_cast<unsigned int>(binTables.size());
  m_maxPairsPerBin1 = 0;
  m_bin_pairs.resize(m_nBinPairs);
  m_pair_group_begin.resize(m_nBinPairs);
  unsigned int run = 0;
  for (unsigned int i = 0; i < m_nBinPairs; i++) {
    const bool same_bin1 =
        (i > 0) && (binTables[i - 1].first == binTables[i].first);
    run = same_bin1 ? run + 1 : 1;
    m_maxPairsPerBin1 = std::max(m_maxPairsPerBin1, run);
    m_bin_pairs[i] = uint2{binTables[i].first, binTables[i].second};
    m_pair_group_begin[i] = same_bin1 ? m_pair_group_begin[i - 1] : i;
  }
}

auto gbts_seeding_algorithm::operator()(
    const edm::spacepoint_collection::const_view& spacepoints,
    const edm::measurement_collection::const_view& measurements) const
    -> output_type {
  unsigned int nSp;
  if (mr().host) {
    vecmem::async_size size = copy().get_size(spacepoints, *(mr().host));
    // Here we could give control back to the caller, once our
    // code allows for it. (coroutines...)
    nSp = size.get();
  } else {
    nSp = copy().get_size(spacepoints);
  }
  TRACCC_DEBUG("nSp " << nSp);
  if (nSp == 0) {
    TRACCC_WARNING("No spacepoints were found in the event");
    return {0, mr().main};
  }

  // Stage 1: bin spacepoints and create nodes with the parameters (eta, phi,
  // r, z).
  node_making_output nodes = make_nodes(spacepoints, measurements, nSp);
  if (nodes.nNodes == 0) {
    // No nodes survived spacepoint counting -> no seeds.
    return {0, mr().main};
  }

  // Named counters shared by the graph-making and seed-extraction stages.
  vecmem::data::vector_buffer<unsigned int> counters_buf(
      gbts_counter::nCounters, mr().main);
  copy().setup(counters_buf)->ignore();
  copy().memset(counters_buf, 0)->ignore();
  vecmem::vector<unsigned int> h_counters(gbts_counter::nCounters,
                                          mr().host ? mr().host : &(mr().main));

  // Stage 2: graph. The per-node buffers are moved in so they are released
  // when create_edges returns, along with the edge transients.
  graph_making_output graph =
      create_edges(std::move(nodes.node_params), std::move(nodes.node_phi),
                   std::move(nodes.node_index), std::move(nodes.bin_rads),
                   std::move(nodes.eta_bin_offsets), nodes.nNodes, nSp,
                   counters_buf, h_counters);
  if (graph.nConnectedEdges == 0) {
    // No connected edges survived graph making -> no seeds.
    return {0, mr().main};
  }

  // Stage 3: Create seeds from the graph edges.
  return extract_seeds(graph.output_graph, nodes.reducedSP,
                       graph.nConnectedEdges, nSp, h_counters);
}

}  // namespace traccc::device
