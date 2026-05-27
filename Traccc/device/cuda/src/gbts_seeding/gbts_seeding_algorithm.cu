/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "../utils/cuda_error_handling.hpp"
#include "../utils/utils.hpp"
#include "./kernels/GbtsGraphMakingKernels.cuh"
#include "./kernels/GbtsGraphProcessingKernels.cuh"
#include "./kernels/GbtsNodesMakingKernels.cuh"
#include "traccc/cuda/gbts_seeding/gbts_seeding_algorithm.hpp"

// C++ include(s)
#include <ranges>

namespace traccc::cuda {

struct gbts_ctx {
    // counters
    unsigned int nSp{};
    unsigned int nNodes{};
    unsigned int nUsedBinPairs{};
    unsigned int nMaxEdges{};

    unsigned int nEdges{};
    unsigned int nConnections{};
    unsigned int nConnectedEdges{};
    unsigned int nSeeds{};
    // nEdges, nConnections, nConnectedEdges, .., nSeeds
    unsigned int* d_counters{};

    // device side graph building cuts
    gbts_graph_building_params* d_graph_building_params;

    // node making and binning
    int* d_layerCounts{};
    short* d_spacepointsLayer{};
    // begin_idx + 1 for the surfaceToLayerMap or -layerBin if one to one
    short* d_volumeToLayerMap{};
    uint2* d_surfaceToLayerMap{};  // surface_index, layerBin
    char* d_layerType{};
    // conversion to original sp from post layer binning index
    int* d_original_sp_idx{};
    // conversion to orignal sp/node index from post binning index
    int* d_node_index{};

    // x,y,z,cluster width in eta
    float4* d_reducedSP{};
    // layer binned reducedSP
    float4* d_sp_params{};

    int2* d_layer_info{};
    float2* d_layer_geo{};

    int* d_node_eta_index{};
    int* d_node_phi_index{};

    int* d_eta_phi_histo{};     // for data binning
    int* d_phi_cusums{};        // for data binning
    int* d_eta_node_counter{};  // for data binning

    int2* d_eta_bin_views{};  // views of the nodes
    // eta-bin views of the node_params array
    std::unique_ptr<int[]> h_eta_bin_views{};

    float2* d_bin_rads{};  // minimum and maximum r of nodes inside an eta-bin
    std::unique_ptr<float[]> h_bin_rads{};

    uint4* d_bin_pair_views{};
    std::unique_ptr<int[]> h_bin_pair_views{};

    std::unique_ptr<float[]> h_bin_pair_dphi{};
    float* d_bin_pair_dphi{};
    // node making output
    float* d_node_params{};

    // GraphMaking
    int2* d_edge_nodes{};
    kernels::half4* d_edge_params{};

    int* d_num_incoming_edges{};
    int* d_edge_links{};

    unsigned char* d_num_neighbours{};
    int* d_reIndexer{};
    int* d_neighbours{};
    // offload this for CPU-side seed extraction
    int* d_output_graph{};

    // message-passing CCA
    // holds flags for edges that need more CCA iterations
    char* d_active_edges{};
    // d_levels[edge_idx] = the maxium length of seeds starting with this edge
    char* d_levels{};
    // #paths, is terminus
    short2* d_outgoing_paths{};

    // seed-extraction walkthrough

    // edge_idx and prev path_store idx forms a uniuqe path through the graph
    int2* d_path_store{};
    int2* d_seed_proposals{};  // int quality and final mini_state_idx
    // first 32 bits are seed quality second 32 bits are seed_proposals index
    unsigned long long int* d_edge_bids{};
    unsigned long long int* d_hit_bids{};
    // 0 as default/is real seed, 1 as maybe seed,
    //-1 as maybe fake seed, -2 as fake
    char* d_seed_ambiguity{};
};

gbts_seeding_algorithm::gbts_seeding_algorithm(
    const gbts_seedfinder_config& cfg, const traccc::memory_resource& mr,
    vecmem::copy& copy, stream& str, std::unique_ptr<const Logger> logger)
    : messaging(logger->clone()),
      m_config(cfg),
      m_mr(mr),
      m_copy(copy),
      m_stream(str) {}

gbts_seeding_algorithm::output_type gbts_seeding_algorithm::operator()(
    const traccc::edm::spacepoint_collection::const_view& spacepoints,
    const edm::measurement_collection::const_view& measurements) const {

    gbts_ctx ctx;

    cudaStream_t stream = details::get_stream(m_stream);

    cudaMalloc(&ctx.d_graph_building_params,
               sizeof(m_config.graph_building_params));
    cudaMemcpyAsync(
        ctx.d_graph_building_params, &m_config.graph_building_params,
        sizeof(m_config.graph_building_params), cudaMemcpyHostToDevice);

    // 0. bin spacepoints by the maping supplied to config.m_surfaceToLayerMap
    ctx.nSp = m_copy.get().get_size(spacepoints);
    if (ctx.nSp == 0) {
        return {0, m_mr.main};
    }

    unsigned int nThreads = 128;
    unsigned int nBlocks = 1 + (ctx.nSp - 1) / nThreads;

    cudaMalloc(&ctx.d_layerCounts, (m_config.nLayers + 1) * sizeof(int));
    cudaMemsetAsync(ctx.d_layerCounts, 0, (m_config.nLayers + 1) * sizeof(int),
                    stream);

    cudaMalloc(&ctx.d_spacepointsLayer, ctx.nSp * sizeof(short));
    cudaMalloc(&ctx.d_reducedSP, ctx.nSp * sizeof(float4));

    cudaMalloc(&ctx.d_volumeToLayerMap,
               sizeof(short) * m_config.volumeToLayerMap.size());

    cudaMemcpyAsync(ctx.d_volumeToLayerMap, m_config.volumeToLayerMap.data(),
                    sizeof(short) * m_config.volumeToLayerMap.size(),
                    cudaMemcpyHostToDevice, stream);

    if (m_config.surfaceToLayerMap.size() != 0) {
        cudaMalloc(&ctx.d_surfaceToLayerMap,
                   sizeof(uint2) * m_config.surfaceToLayerMap.size());

        cudaMemcpyAsync(ctx.d_surfaceToLayerMap,
                        m_config.surfaceToLayerMap.data(),
                        sizeof(uint2) * m_config.surfaceToLayerMap.size(),
                        cudaMemcpyHostToDevice, stream);
    }  // may be zero and correct, volumeMapSize, nLayers are checked at config

    cudaMalloc(&ctx.d_layerType, sizeof(char) * m_config.nLayers);
    cudaMemcpyAsync(ctx.d_layerType, m_config.layerInfo.type.data(),
                    sizeof(char) * m_config.nLayers, cudaMemcpyHostToDevice,
                    stream);

    kernels::count_sp_by_layer<<<nBlocks, nThreads, 0, stream>>>(
        spacepoints, measurements, ctx.d_volumeToLayerMap,
        ctx.d_surfaceToLayerMap, ctx.d_layerType, ctx.d_reducedSP,
        ctx.d_layerCounts, ctx.d_spacepointsLayer,
        m_config.graph_building_params.type1_max_width, ctx.nSp,
        m_config.volumeToLayerMap.size(), m_config.surfaceToLayerMap.size());

    cudaStreamSynchronize(stream);

    cudaFree(ctx.d_volumeToLayerMap);
    cudaFree(ctx.d_surfaceToLayerMap);
    cudaFree(ctx.d_layerType);

    // prefix sum layerCounts
    std::unique_ptr<int[]> layerCounts =
        std::make_unique<int[]>(m_config.nLayers + 1);

    cudaMemcpyAsync(layerCounts.get(), ctx.d_layerCounts,
                    (m_config.nLayers + 1) * sizeof(int),
                    cudaMemcpyDeviceToHost, stream);

    for (size_t layer = 0; layer < m_config.nLayers; layer++) {
        layerCounts[layer + 1] += layerCounts[layer];
    }
    cudaMemcpyAsync(ctx.d_layerCounts, layerCounts.get(),
                    m_config.nLayers * sizeof(int), cudaMemcpyHostToDevice,
                    stream);

    ctx.nNodes = static_cast<unsigned int>(layerCounts[m_config.nLayers]);
    if (ctx.nNodes == 0)
        return {0, m_mr.main};
    layerCounts.reset();

    cudaMalloc(&ctx.d_sp_params, ctx.nSp * sizeof(float4));
    cudaMalloc(&ctx.d_original_sp_idx, ctx.nSp * sizeof(int));

    kernels::bin_sp_by_layer<<<nBlocks, nThreads, 0, stream>>>(
        ctx.d_sp_params, ctx.d_reducedSP, ctx.d_layerCounts,
        ctx.d_spacepointsLayer, ctx.d_original_sp_idx, ctx.nSp);

    cudaStreamSynchronize(stream);
    cudaError_t error = cudaGetLastError();

    if (error != cudaSuccess) {
        TRACCC_ERROR("spacepoint layer binning: CUDA error: "
                     << cudaGetErrorString(error));
        return {0, m_mr.main};
    }

    cudaFree(ctx.d_spacepointsLayer);

    // 1. histogram spacepoints by layer->eta->phi and convert to nodes
    // do this in config setup?
    cudaMalloc(&ctx.d_layer_info, sizeof(int2) * m_config.nLayers);
    cudaMemcpyAsync(ctx.d_layer_info, m_config.layerInfo.info.data(),
                    sizeof(int2) * m_config.nLayers, cudaMemcpyHostToDevice,
                    stream);

    cudaMalloc(&ctx.d_layer_geo, sizeof(float2) * m_config.nLayers);
    cudaMemcpyAsync(ctx.d_layer_geo, m_config.layerInfo.geo.data(),
                    sizeof(float2) * m_config.nLayers, cudaMemcpyHostToDevice,
                    stream);

    cudaMalloc(&ctx.d_node_phi_index, sizeof(int) * ctx.nNodes);

    nThreads = 256;
    unsigned int nNodesPerBlock = nThreads * 64;

    nBlocks = 1 + (ctx.nNodes - 1) / nNodesPerBlock;

    kernels::node_phi_binning_kernel<<<nBlocks, nThreads, 0, stream>>>(
        ctx.d_sp_params, ctx.d_node_phi_index, nNodesPerBlock, ctx.nNodes,
        m_config.n_phi_bins);

    cudaStreamSynchronize(stream);

    cudaMalloc(&ctx.d_node_eta_index, sizeof(int) * ctx.nNodes);

    nBlocks = m_config.nLayers;

    kernels::node_eta_binning_kernel<<<nBlocks, nThreads, 0, stream>>>(
        ctx.d_sp_params, ctx.d_layer_info, ctx.d_layer_geo,
        ctx.d_node_eta_index, ctx.d_layerCounts, m_config.nLayers);

    cudaStreamSynchronize(stream);

    cudaFree(ctx.d_layerCounts);
    cudaFree(ctx.d_layer_info);
    cudaFree(ctx.d_layer_geo);

    error = cudaGetLastError();

    if (error != cudaSuccess) {
        TRACCC_ERROR(
            "eta-phi binning: CUDA error: " << cudaGetErrorString(error));
        return {0, m_mr.main};
    }
    size_t hist_size = sizeof(int) * m_config.n_eta_bins * m_config.n_phi_bins;
    cudaMalloc(&ctx.d_eta_phi_histo, hist_size);
    cudaMemsetAsync(ctx.d_eta_phi_histo, 0, hist_size, stream);
    cudaMalloc(&ctx.d_phi_cusums, hist_size);

    nBlocks = 1 + (ctx.nNodes - 1) / nNodesPerBlock;

    kernels::eta_phi_histo_kernel<<<nBlocks, nThreads, 0, stream>>>(
        ctx.d_node_phi_index, ctx.d_node_eta_index, ctx.d_eta_phi_histo,
        nNodesPerBlock, ctx.nNodes, m_config.n_phi_bins);

    cudaStreamSynchronize(stream);
    error = cudaGetLastError();

    if (error != cudaSuccess) {
        TRACCC_ERROR(
            "eta-phi histo: CUDA error: " << cudaGetErrorString(error));
        return {0, m_mr.main};
    }

    cudaMalloc(&ctx.d_eta_node_counter, sizeof(int) * m_config.n_eta_bins);

    unsigned int nBinsPerBlock = 128;

    nThreads = nBinsPerBlock;

    nBlocks = 1 + (m_config.n_eta_bins - 1) / nBinsPerBlock;

    kernels::eta_phi_counting_kernel<<<nBlocks, nThreads, 0, stream>>>(
        ctx.d_eta_phi_histo, ctx.d_eta_node_counter, ctx.d_phi_cusums,
        nBinsPerBlock, m_config.n_eta_bins, m_config.n_phi_bins);

    cudaStreamSynchronize(stream);
    cudaFree(ctx.d_eta_phi_histo);

    error = cudaGetLastError();

    if (error != cudaSuccess) {
        TRACCC_ERROR(
            "eta-phi counting: CUDA error: " << cudaGetErrorString(error));
        return {0, m_mr.main};
    }

    std::unique_ptr<int[]> eta_sums =
        std::make_unique<int[]>(m_config.n_eta_bins);

    cudaMemcpyAsync(&eta_sums[0], &ctx.d_eta_node_counter[0],
                    sizeof(int) * m_config.n_eta_bins, cudaMemcpyDeviceToHost,
                    stream);

    cudaStreamSynchronize(stream);

    for (unsigned int k = 0; k < m_config.n_eta_bins; k++) {
        eta_sums[k + 1] += eta_sums[k];
    }
    // send back
    cudaMemcpyAsync(&ctx.d_eta_node_counter[0], &eta_sums[0],
                    sizeof(int) * m_config.n_eta_bins, cudaMemcpyHostToDevice,
                    stream);

    ctx.h_eta_bin_views = std::make_unique<int[]>(2 * m_config.n_eta_bins);

    for (unsigned view_idx = 0; view_idx < m_config.n_eta_bins; view_idx++) {
        unsigned int pos = 2 * view_idx;
        ctx.h_eta_bin_views[pos] = (view_idx == 0) ? 0 : eta_sums[view_idx - 1];
        ctx.h_eta_bin_views[pos + 1] = eta_sums[view_idx];
    }
    eta_sums.reset();

    cudaStreamSynchronize(stream);

    kernels::eta_phi_prefix_sum_kernel<<<nBlocks, nThreads, 0, stream>>>(
        ctx.d_eta_node_counter, ctx.d_phi_cusums, nBinsPerBlock,
        m_config.n_eta_bins, m_config.n_phi_bins);

    cudaStreamSynchronize(stream);
    cudaFree(ctx.d_eta_node_counter);

    error = cudaGetLastError();

    if (error != cudaSuccess) {
        TRACCC_ERROR(
            "eta-phi cusum: CUDA error: " << cudaGetErrorString(error));
        return {0, m_mr.main};
    }

    cudaMalloc(&ctx.d_node_params, 5 * sizeof(float) * ctx.nNodes);
    cudaMalloc(&ctx.d_node_index, sizeof(int) * ctx.nNodes);

    nThreads = 256;
    nNodesPerBlock = nThreads * 64;

    nBlocks = 1 + (ctx.nNodes - 1) / nNodesPerBlock;

    kernels::node_sorting_kernel<<<nBlocks, nThreads, 0, stream>>>(
        ctx.d_sp_params, ctx.d_node_eta_index, ctx.d_node_phi_index,
        ctx.d_phi_cusums, ctx.d_node_params, ctx.d_node_index,
        ctx.d_original_sp_idx, ctx.d_graph_building_params, nNodesPerBlock,
        ctx.nNodes, m_config.n_phi_bins);

    cudaStreamSynchronize(stream);
    cudaFree(ctx.d_sp_params);
    cudaFree(ctx.d_original_sp_idx);
    cudaFree(ctx.d_phi_cusums);
    cudaFree(ctx.d_node_eta_index);
    cudaFree(ctx.d_node_phi_index);

    error = cudaGetLastError();

    if (error != cudaSuccess) {
        TRACCC_ERROR("node sorting: CUDA error: " << cudaGetErrorString(error));
        return {0, m_mr.main};
    }

    cudaMalloc(&ctx.d_eta_bin_views, sizeof(int2) * m_config.n_eta_bins);
    cudaMalloc(&ctx.d_bin_rads, sizeof(float2) * m_config.n_eta_bins);

    cudaMemcpyAsync(&ctx.d_eta_bin_views[0], ctx.h_eta_bin_views.get(),
                    2 * m_config.n_eta_bins * sizeof(int),
                    cudaMemcpyHostToDevice, stream);

    cudaStreamSynchronize(stream);

    nBinsPerBlock = 128;

    nThreads = nBinsPerBlock;

    nBlocks = 1 + (m_config.n_eta_bins - 1) / nBinsPerBlock;

    kernels::minmax_rad_kernel<<<nBlocks, nThreads, 0, stream>>>(
        ctx.d_eta_bin_views, ctx.d_node_params, ctx.d_bin_rads, nBinsPerBlock,
        m_config.n_eta_bins);

    cudaStreamSynchronize(stream);
    cudaFree(ctx.d_eta_bin_views);

    error = cudaGetLastError();

    if (error != cudaSuccess) {
        TRACCC_ERROR("node sorting: CUDA error: " << cudaGetErrorString(error));
        return {0, m_mr.main};
    }
    ctx.h_bin_rads = std::make_unique<float[]>(2 * m_config.n_eta_bins);

    cudaMemcpyAsync(ctx.h_bin_rads.get(), &ctx.d_bin_rads[0],
                    2 * sizeof(float) * m_config.n_eta_bins,
                    cudaMemcpyDeviceToHost, stream);

    cudaStreamSynchronize(stream);

    cudaFree(ctx.d_bin_rads);

    // 2. prepare input for the graph making part of the code:

    int int_nBinPairs = 0;  // the number of eta bin pairs

    for (std::pair<unsigned int, unsigned int> binPair : m_config.binTables) {
        // loop over bin pairs defined by the layer
        // connection table and geometry settings

        int bin1_begin = ctx.h_eta_bin_views[2 * binPair.first];
        int bin1_end = ctx.h_eta_bin_views[2 * binPair.first + 1];
        // large bins will be split into smaller sub-views

        int nNodesInBin1 = bin1_end - bin1_begin;
        if (bin1_begin > bin1_end) {
            nNodesInBin1 = bin1_begin - bin1_end;
        }
        int_nBinPairs +=
            1 + (nNodesInBin1 - 1) /
                    traccc::device::gbts_consts::node_buffer_length;
    }
    unsigned int nBinPairs = static_cast<unsigned int>(int_nBinPairs);

    ctx.h_bin_pair_views = std::make_unique<int[]>(4 * nBinPairs);
    ctx.h_bin_pair_dphi = std::make_unique<float[]>(nBinPairs);

    unsigned int pairIdx = 0;
    for (std::pair<unsigned int, unsigned int> binPair : m_config.binTables) {
        float rb1 = ctx.h_bin_rads[2 * binPair.first];  // min radius

        int begin_bin1 = ctx.h_eta_bin_views[2 * binPair.first];
        int end_bin1 = ctx.h_eta_bin_views[2 * binPair.first + 1];
        // skip empty pairs
        if (begin_bin1 == end_bin1)
            continue;
        if (ctx.h_eta_bin_views[2 * binPair.second] ==
            ctx.h_eta_bin_views[2 * binPair.second + 1])
            continue;

        float rb2 = ctx.h_bin_rads[2 * binPair.second + 1];  // max radius

        // max radius of bin2 - min radius of bin1
        float maxDeltaR = std::fabs(rb2 - rb1);

        float deltaPhi = m_config.graph_building_params.min_delta_phi +
                         m_config.graph_building_params.dphi_coeff * maxDeltaR;

        if (maxDeltaR < 60) {
            deltaPhi =
                m_config.graph_building_params.min_delta_phi_low_dr +
                m_config.graph_building_params.dphi_coeff_low_dr * maxDeltaR;
        }
        // splitting large bins into more consistent sizes

        int currBegin_bin1 = begin_bin1;

        int currEnd_bin1 =
            end_bin1 < traccc::device::gbts_consts::node_buffer_length
                ? end_bin1
                : begin_bin1 + traccc::device::gbts_consts::node_buffer_length;

        for (; currEnd_bin1 < end_bin1;
             currEnd_bin1 += traccc::device::gbts_consts::node_buffer_length,
             pairIdx++) {
            unsigned int offset = 4 * pairIdx;

            ctx.h_bin_pair_views[offset] = currBegin_bin1;
            ctx.h_bin_pair_views[1 + offset] = currEnd_bin1;
            ctx.h_bin_pair_views[2 + offset] =
                ctx.h_eta_bin_views[2 * binPair.second];
            ctx.h_bin_pair_views[3 + offset] =
                ctx.h_eta_bin_views[2 * binPair.second + 1];
            ctx.h_bin_pair_dphi[pairIdx] = deltaPhi;

            currBegin_bin1 = currEnd_bin1;
        }
        currEnd_bin1 = end_bin1;

        unsigned int offset = 4 * pairIdx;

        ctx.h_bin_pair_views[offset] = currBegin_bin1;
        ctx.h_bin_pair_views[1 + offset] = currEnd_bin1;
        ctx.h_bin_pair_views[2 + offset] =
            ctx.h_eta_bin_views[2 * binPair.second];
        ctx.h_bin_pair_views[3 + offset] =
            ctx.h_eta_bin_views[2 * binPair.second + 1];
        ctx.h_bin_pair_dphi[pairIdx] = deltaPhi;
        pairIdx++;
    }
    ctx.nUsedBinPairs = pairIdx;
    if (ctx.nUsedBinPairs == 0)
        return {0, m_mr.main};
    ctx.h_eta_bin_views.reset();
    // allocate memory and copy bin pair views and phi cuts to GPU

    unsigned int data_size = ctx.nUsedBinPairs * sizeof(uint4);

    cudaMalloc(&ctx.d_bin_pair_views, data_size);
    cudaMemcpyAsync(&ctx.d_bin_pair_views[0], &ctx.h_bin_pair_views[0],
                    data_size, cudaMemcpyHostToDevice, stream);

    data_size = ctx.nUsedBinPairs * sizeof(float);

    cudaMalloc(&ctx.d_bin_pair_dphi, data_size);
    cudaMemcpyAsync(&ctx.d_bin_pair_dphi[0], &ctx.h_bin_pair_dphi[0], data_size,
                    cudaMemcpyHostToDevice, stream);

    cudaMalloc(&ctx.d_counters, sizeof(unsigned int) * 12);
    cudaMemsetAsync(ctx.d_counters, 0, sizeof(unsigned int) * 12, stream);

    cudaStreamSynchronize(stream);

    // 2. Find edges between spacepoint pairs
    ctx.nMaxEdges = m_config.max_edges_factor * ctx.nNodes;
    cudaMalloc(&ctx.d_edge_params, sizeof(kernels::half4) * ctx.nMaxEdges);
    cudaMalloc(&ctx.d_edge_nodes, sizeof(int2) * ctx.nMaxEdges);

    cudaMalloc(&ctx.d_num_incoming_edges, sizeof(int) * (ctx.nNodes + 1));
    cudaMemsetAsync(ctx.d_num_incoming_edges, 0, sizeof(int) * (ctx.nNodes + 1),
                    stream);

    nBlocks = ctx.nUsedBinPairs;
    nThreads = 128;

    kernels::graphEdgeMakingKernel<<<nBlocks, nThreads, 0, stream>>>(
        ctx.d_bin_pair_views, ctx.d_bin_pair_dphi, ctx.d_node_params,
        ctx.d_graph_building_params, ctx.d_counters, ctx.d_edge_nodes,
        ctx.d_edge_params, ctx.d_num_incoming_edges, ctx.nMaxEdges,
        m_config.n_phi_bins);

    cudaStreamSynchronize(stream);
    cudaFree(ctx.d_node_params);
    cudaFree(ctx.d_bin_pair_views);
    cudaFree(ctx.d_bin_pair_dphi);

    error = cudaGetLastError();

    if (error != cudaSuccess) {
        TRACCC_ERROR("edge making: CUDA error: " << cudaGetErrorString(error));
        return {0, m_mr.main};
    }

    cudaMemcpyAsync(&ctx.nEdges, ctx.d_counters, sizeof(unsigned int),
                    cudaMemcpyDeviceToHost, stream);

    TRACCC_DEBUG("Created " << ctx.nEdges << " edges with a cap of "
                            << ctx.nMaxEdges);

    if (ctx.nEdges > ctx.nMaxEdges) {
        TRACCC_ERROR("Number of edges exceeds the maximum allowed, Removing "
                     << ctx.nEdges - ctx.nMaxEdges << " edges");
        ctx.nEdges = ctx.nMaxEdges;
    } else if (ctx.nEdges == 0) {
        return {0, m_mr.main};
    }

    std::unique_ptr<int[]> cusum = std::make_unique<int[]>(ctx.nNodes + 1);

    data_size = (ctx.nNodes + 1) * sizeof(int);

    cudaMemcpyAsync(&cusum[0], ctx.d_num_incoming_edges, data_size,
                    cudaMemcpyDeviceToHost, stream);

    cudaStreamSynchronize(stream);

    for (unsigned int k = 0; k < ctx.nNodes; k++)
        cusum[k + 1] += cusum[k];

    cudaMemcpyAsync(ctx.d_num_incoming_edges, &cusum[0], data_size,
                    cudaMemcpyHostToDevice, stream);

    cudaStreamSynchronize(stream);

    cusum.reset();

    // 3. link edges and nodes

    data_size = ctx.nEdges * sizeof(int);

    cudaMalloc(&ctx.d_edge_links, data_size);

    nThreads = 256;
    nBlocks = 1 + (ctx.nEdges - 1) / nThreads;

    kernels::graphEdgeLinkingKernel<<<nBlocks, nThreads, 0, stream>>>(
        ctx.d_edge_nodes, ctx.d_edge_links, ctx.d_num_incoming_edges,
        ctx.nEdges);

    cudaStreamSynchronize(stream);

    error = cudaGetLastError();

    if (error != cudaSuccess) {
        TRACCC_ERROR("edge linking: CUDA error: " << cudaGetErrorString(error));
        return {0, m_mr.main};
    }

    // 4. edge matching to create edge-to-edge connections

    data_size = ctx.nEdges * sizeof(unsigned char);

    cudaMalloc(&ctx.d_num_neighbours, data_size);
    cudaMemsetAsync(ctx.d_num_neighbours, 0, data_size, stream);

    data_size = ctx.nEdges * sizeof(int);

    cudaMalloc(&ctx.d_reIndexer, data_size);
    cudaMemsetAsync(ctx.d_reIndexer, 0xFF, data_size, stream);

    data_size = m_config.max_num_neighbours * ctx.nEdges * sizeof(int);
    cudaMalloc(&ctx.d_neighbours, data_size);

    kernels::graphEdgeMatchingKernel<<<nBlocks, nThreads, 0, stream>>>(
        ctx.d_graph_building_params, ctx.d_edge_params, ctx.d_edge_nodes,
        ctx.d_num_incoming_edges, ctx.d_edge_links, ctx.d_num_neighbours,
        ctx.d_neighbours, ctx.d_reIndexer, ctx.d_counters, ctx.nEdges,
        m_config.max_num_neighbours);

    cudaStreamSynchronize(stream);
    cudaFree(ctx.d_num_incoming_edges);
    cudaFree(ctx.d_edge_links);
    cudaFree(ctx.d_edge_params);

    error = cudaGetLastError();

    if (error != cudaSuccess) {
        TRACCC_ERROR(
            "edge matching: CUDA error: " << cudaGetErrorString(error));
        return {0, m_mr.main};
    }

    // 5. Edge re-indexing to keep only edges involved in any connection

    kernels::edgeReIndexingKernel<<<nBlocks, nThreads, 0, stream>>>(
        ctx.d_reIndexer, ctx.d_counters, ctx.nEdges);

    cudaStreamSynchronize(stream);

    error = cudaGetLastError();

    if (error != cudaSuccess) {
        TRACCC_ERROR(
            "edge re-indexing: CUDA error: " << cudaGetErrorString(error));
        return {0, m_mr.main};
    }

    unsigned int nStats[3];

    cudaMemcpyAsync(&nStats[0], ctx.d_counters, 3 * sizeof(unsigned int),
                    cudaMemcpyDeviceToHost, stream);

    ctx.nConnections = nStats[1];
    ctx.nConnectedEdges = nStats[2];

    TRACCC_DEBUG("created " << ctx.nConnections << " edge links, found "
                            << ctx.nConnectedEdges
                            << " connected edges for seed extraction");
    if (ctx.nConnectedEdges == 0)
        return {0, m_mr.main};

    unsigned int nIntsPerEdge = 2 + 1 + m_config.max_num_neighbours;

    data_size = ctx.nConnectedEdges * nIntsPerEdge * sizeof(int);

    cudaMalloc(&ctx.d_output_graph, data_size);

    nThreads = 256;
    unsigned int nEdgesPerBlock = nThreads * 64;

    nBlocks = 1 + (ctx.nEdges - 1) / nEdgesPerBlock;

    kernels::graphCompressionKernel<<<nBlocks, nThreads, 0, stream>>>(
        ctx.d_node_index, ctx.d_edge_nodes, ctx.d_num_neighbours,
        ctx.d_neighbours, ctx.d_reIndexer, ctx.d_output_graph, nEdgesPerBlock,
        ctx.nEdges, m_config.max_num_neighbours);

    cudaStreamSynchronize(stream);

    cudaFree(ctx.d_edge_nodes);
    cudaFree(ctx.d_node_index);
    cudaFree(ctx.d_reIndexer);
    cudaFree(ctx.d_num_neighbours);
    cudaFree(ctx.d_neighbours);

    error = cudaGetLastError();
    if (error != cudaSuccess) {
        TRACCC_ERROR(
            "graph compression: CUDA error: " << cudaGetErrorString(error));
        return {0, m_mr.main};
    }

    // 6. Find longest segments with CCA

    data_size = ctx.nConnectedEdges * sizeof(int);

    cudaMalloc(&ctx.d_active_edges, data_size);

    data_size = 2 * ctx.nConnectedEdges * sizeof(unsigned char);

    // old levels and new levels are kept in opposite halves of the array
    cudaMalloc(&ctx.d_levels, data_size);
    // initalize to 1 so level counts the maxium number of edge segments
    //  for a seed originating at the edge
    cudaMemsetAsync(ctx.d_levels, 0x1, data_size, stream);

    cudaMalloc(&ctx.d_outgoing_paths, ctx.nConnectedEdges * sizeof(short2));

    unsigned int nEdgesLeft = ctx.nConnectedEdges;

    cudaMemcpyAsync(&ctx.d_counters[3], &nEdgesLeft, sizeof(unsigned int),
                    cudaMemcpyHostToDevice, stream);

    if (nEdgesLeft == 0)
        return {0, m_mr.main};

    cudaStreamSynchronize(stream);

    nThreads = 128;
    nBlocks = 1 + (nEdgesLeft - 1) / nThreads;
    for (int iter = 0; iter < traccc::device::gbts_consts::max_cca_iter;
         iter++) {
        kernels::CCA_IterationKernel<<<nBlocks, nThreads, 0, stream>>>(
            ctx.d_output_graph, ctx.d_levels, ctx.d_active_edges,
            ctx.d_outgoing_paths, iter, ctx.nConnectedEdges,
            m_config.max_num_neighbours, m_config.minLevel);

        cudaStreamSynchronize(stream);
    }

    cudaStreamSynchronize(stream);

    cudaFree(ctx.d_active_edges);

    error = cudaGetLastError();

    if (error != cudaSuccess) {
        TRACCC_ERROR(
            "message-passing CCA: CUDA error: " << cudaGetErrorString(error));
        return {0, m_mr.main};
    }

    nThreads = 128;
    nBlocks = 1 + (ctx.nConnectedEdges - 1) / nThreads;

    cudaStreamSynchronize(stream);

    kernels::count_terminus_edges<<<nBlocks, nThreads, 0, stream>>>(
        ctx.d_path_store, ctx.d_outgoing_paths, ctx.d_counters,
        ctx.nConnectedEdges);

    cudaStreamSynchronize(stream);

    // nPaths to terminus, nTerminusEdges
    unsigned int path_sizes[2];
    cudaMemcpyAsync(&path_sizes[0], &ctx.d_counters[6],
                    2 * sizeof(unsigned int), cudaMemcpyDeviceToHost, stream);
    unsigned int pathsPerTerminus = 1 + (path_sizes[0] - 1) / path_sizes[1];

    TRACCC_DEBUG(path_sizes[0] << "size of path store | nTerminusEdges "
                               << path_sizes[1]);

    cudaMalloc(&ctx.d_path_store,
               (path_sizes[0] + path_sizes[1]) * sizeof(int2));
    cudaMalloc(&ctx.d_seed_proposals, path_sizes[0] * sizeof(int2));
    cudaMalloc(&ctx.d_seed_ambiguity, path_sizes[0] * sizeof(char));

    cudaMalloc(&ctx.d_edge_bids,
               ctx.nConnectedEdges * sizeof(unsigned long long int));
    cudaMemsetAsync(ctx.d_edge_bids, 0,
                    ctx.nConnectedEdges * sizeof(unsigned long long int),
                    stream);

    nThreads = 128;
    kernels::add_terminus_to_path_store<<<nBlocks, nThreads, 0, stream>>>(
        ctx.d_path_store, ctx.d_outgoing_paths, ctx.d_counters,
        ctx.nConnectedEdges);

    unsigned int terminusPerBlock = std::min(
        nThreads, 1 + (traccc::device::gbts_consts::live_path_buffer - 1) /
                          pathsPerTerminus);
    nBlocks = 1 + (path_sizes[1] - 1) / terminusPerBlock;

    kernels::fill_path_store<<<nBlocks, nThreads, 0, stream>>>(
        ctx.d_path_store, ctx.d_output_graph, ctx.d_levels, ctx.d_counters,
        path_sizes[1], terminusPerBlock, m_config.max_num_neighbours,
        path_sizes[0] + path_sizes[1]);

    nThreads = 128;
    nBlocks = 1 + (path_sizes[0] + path_sizes[1] - 1) / nThreads;

    kernels::fit_segments<<<nBlocks, nThreads, 0, stream>>>(
        ctx.d_reducedSP, ctx.d_output_graph, ctx.d_path_store,
        ctx.d_seed_proposals, ctx.d_edge_bids, ctx.d_seed_ambiguity,
        ctx.d_levels, ctx.d_counters, path_sizes[1], m_config.minLevel,
        m_config.max_num_neighbours, m_config.seed_extraction_params);

    unsigned int nProps = 0;
    cudaMemcpyAsync(&nProps, &ctx.d_counters[8], sizeof(unsigned int),
                    cudaMemcpyDeviceToHost, stream);

    TRACCC_DEBUG("nProps " << nProps);

    cudaStreamSynchronize(stream);

    error = cudaGetLastError();
    if (error != cudaSuccess) {
        TRACCC_ERROR(
            "seed extraction: CUDA error: " << cudaGetErrorString(error));
        return {0, m_mr.main};
    }

    cudaFree(ctx.d_levels);
    cudaFree(ctx.d_outgoing_paths);

    if (nProps == 0) {
        return {0, m_mr.main};
    }

    nThreads = 128;
    nBlocks = 1 + (nProps - 1) / nThreads;

    kernels::reset_edge_bids<<<nBlocks, nThreads, 0, stream>>>(
        ctx.d_path_store, ctx.d_seed_proposals, ctx.d_edge_bids,
        ctx.d_seed_ambiguity, ctx.d_counters, -1);

    for (int round = 0; round < 5; ++round) {

        cudaMemsetAsync(ctx.d_edge_bids, 0,
                        ctx.nConnectedEdges * sizeof(unsigned long long int),
                        stream);

        kernels::seeds_rebid_for_edges<<<nBlocks, nThreads, 0, stream>>>(
            ctx.d_path_store, ctx.d_seed_proposals, ctx.d_edge_bids,
            ctx.d_seed_ambiguity, nProps);

        kernels::reset_edge_bids<<<nBlocks, nThreads, 0, stream>>>(
            ctx.d_path_store, ctx.d_seed_proposals, ctx.d_edge_bids,
            ctx.d_seed_ambiguity, ctx.d_counters, round);
    }
    unsigned int nRejectedProps = 0;
    cudaMemcpyAsync(&nRejectedProps, &ctx.d_counters[9], sizeof(unsigned int),
                    cudaMemcpyDeviceToHost, stream);
    ctx.nSeeds = nProps - nRejectedProps;

    TRACCC_DEBUG("Rejecetd " << nRejectedProps << " out of " << nProps
                             << " seed proposals");

    cudaFree(ctx.d_edge_bids);
    cudaFree(ctx.d_counters);
    cudaFree(ctx.d_graph_building_params);

    // 8. convert to 3sp seeds and make output buffer
    // allocate extra seed space for hit permutation
    edm::seed_collection::buffer output_seeds(
        2 * ctx.nSeeds, m_mr.main, vecmem::data::buffer_type::resizable);
    m_copy.get().setup(output_seeds)->wait();

    cudaMalloc(&ctx.d_hit_bids, ctx.nSp * sizeof(unsigned long long int));
    cudaMemsetAsync(ctx.d_hit_bids, 0, ctx.nSp * sizeof(unsigned long long int),
                    stream);

    nThreads = 128;
    nBlocks = 1 + (ctx.nSeeds - 1) / nThreads;

    kernels::seeds_bid_for_hits<<<nBlocks, nThreads, 0, stream>>>(
        ctx.d_output_graph, ctx.d_seed_proposals, ctx.d_path_store,
        ctx.d_seed_ambiguity, ctx.d_hit_bids, nProps,
        static_cast<int>(1 + 2 + m_config.max_num_neighbours));

    kernels::gbts_seed_conversion_kernel<<<nBlocks, nThreads, 0, stream>>>(
        ctx.d_seed_proposals, ctx.d_seed_ambiguity, ctx.d_path_store,
        ctx.d_output_graph, ctx.d_reducedSP, output_seeds, ctx.d_hit_bids,
        nProps, m_config.max_num_neighbours,
        m_config.seed_ambi_params.dropout_dcurv_m,
        m_config.seed_ambi_params.force_dropout_max_curv_m,
        m_config.seed_ambi_params.best_hit_frac,
        m_config.seed_ambi_params.tight_bid_cot_threshold,
        m_config.seed_ambi_params.use_dropout);

    cudaStreamSynchronize(stream);

    ctx.nSeeds = m_copy.get().get_size(output_seeds);

    cudaFree(ctx.d_reducedSP);
    cudaFree(ctx.d_output_graph);
    cudaFree(ctx.d_path_store);
    cudaFree(ctx.d_seed_proposals);
    cudaFree(ctx.d_seed_ambiguity);
    cudaFree(ctx.d_hit_bids);

    error = cudaGetLastError();
    if (error != cudaSuccess) {
        TRACCC_ERROR("Seed ambiguity solving: CUDA error: "
                     << cudaGetErrorString(error));
        return {0, m_mr.main};
    }

    TRACCC_DEBUG("GBTS found " << ctx.nSeeds << " seeds");
    return output_seeds;
}

}  // namespace traccc::cuda
