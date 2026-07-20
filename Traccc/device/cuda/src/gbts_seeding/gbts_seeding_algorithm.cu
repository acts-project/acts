/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "../utils/barrier.hpp"
#include "../utils/cuda_error_handling.hpp"
#include "../utils/thread_id.hpp"
#include "../utils/utils.hpp"
#include "traccc/cuda/gbts_seeding/gbts_seeding_algorithm.hpp"

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

// System include(s).
#include <algorithm>
#include <memory_resource>

// Thrust include(s).
#include <thrust/execution_policy.h>
#include <thrust/scan.h>

namespace traccc::cuda {

namespace kernels {

using float4 = traccc::float4;
using uint2 = traccc::uint2;
using int2 = traccc::int2;

// ---------------------------------------------------------------------------
// Stage 1 — nodes-making kernels
// ---------------------------------------------------------------------------

/// CUDA kernel for running @c traccc::device::gbts_count_spacepoints_by_layer
__global__ void gbts_count_spacepoints_by_layer(
    const device::gbts_count_spacepoints_by_layer_payload payload) {
  device::gbts_count_spacepoints_by_layer(details::thread_id1{}, payload);
}

/// CUDA kernel for running @c traccc::device::gbts_bin_spacepoints
__global__ void gbts_bin_spacepoints(
    const device::gbts_bin_spacepoints_payload payload) {
  device::gbts_bin_spacepoints(details::thread_id1{}, payload);
}

/// CUDA kernel for running @c traccc::device::gbts_count_eta_phi_bins
__global__ void gbts_count_eta_phi_bins(
    const device::gbts_count_eta_phi_bins_payload payload) {
  device::gbts_count_eta_phi_bins(details::thread_id1{}, payload);
}

/// CUDA kernel for running @c traccc::device::gbts_prefix_sum_eta_phi_bins
__global__ void gbts_prefix_sum_eta_phi_bins(
    const device::gbts_prefix_sum_eta_phi_bins_payload payload) {
  device::gbts_prefix_sum_eta_phi_bins(details::thread_id1{}, payload);
}

/// CUDA kernel for running @c traccc::device::gbts_sort_nodes
__global__ void gbts_sort_nodes(const device::gbts_sort_nodes_payload payload) {
  device::gbts_sort_nodes(details::thread_id1{}, payload);
}

/// CUDA kernel for running @c traccc::device::gbts_find_minmax_radius
__global__ void gbts_find_minmax_radius(
    const device::gbts_find_minmax_radius_payload payload) {
  device::gbts_find_minmax_radius(details::thread_id1{}, payload);
}

// ---------------------------------------------------------------------------
// Stage 2 — graph-making kernels
// ---------------------------------------------------------------------------

/// CUDA kernel for running @c traccc::device::gbts_make_graph_edges
__global__ void gbts_make_graph_edges(
    const device::gbts_make_graph_edges_payload payload) {
  __shared__ float phi[traccc::device::gbts_consts::node_buffer_length];
  __shared__ float4 node_pack[traccc::device::gbts_consts::node_buffer_length];
  const traccc::cuda::barrier barrier;

  device::gbts_make_graph_edges(
      details::thread_id1{}, barrier, payload,
      {vecmem::data::vector_view<float>(
           traccc::device::gbts_consts::node_buffer_length, phi),
       vecmem::data::vector_view<float4>(
           traccc::device::gbts_consts::node_buffer_length, node_pack)});
}

/// CUDA kernel for running @c traccc::device::gbts_link_graph_edges
__global__ void gbts_link_graph_edges(
    const device::gbts_link_graph_edges_payload payload) {
  device::gbts_link_graph_edges(details::thread_id1{}, payload);
}

/// CUDA kernel for running @c traccc::device::gbts_match_graph_edges
__global__ void gbts_match_graph_edges(
    const device::gbts_match_graph_edges_payload payload) {
  device::gbts_match_graph_edges(details::thread_id1{}, payload);
}

/// CUDA kernel for running @c traccc::device::gbts_reindex_edges
__global__ void gbts_reindex_edges(
    const device::gbts_reindex_edges_payload payload) {
  device::gbts_reindex_edges(details::thread_id1{}, payload);
}

/// CUDA kernel for running @c traccc::device::gbts_compress_graph
__global__ void gbts_compress_graph(
    const device::gbts_compress_graph_payload payload) {
  device::gbts_compress_graph(details::thread_id1{}, payload);
}

// ---------------------------------------------------------------------------
// Stage 3 — graph-processing kernels
// ---------------------------------------------------------------------------

/// CUDA kernel for running @c traccc::device::gbts_run_cca_iteration
__global__ void gbts_run_cca_iteration(
    const device::gbts_run_cca_iteration_payload payload) {
  device::gbts_run_cca_iteration(details::thread_id1{}, payload);
}

/// CUDA kernel for running @c traccc::device::gbts_count_terminus_edges
__global__ void gbts_count_terminus_edges(
    const device::gbts_count_terminus_edges_payload payload) {
  device::gbts_count_terminus_edges(details::thread_id1{}, payload);
}

/// CUDA kernel for running @c traccc::device::gbts_add_terminus_to_path_store
__global__ void gbts_add_terminus_to_path_store(
    const device::gbts_add_terminus_to_path_store_payload payload) {
  device::gbts_add_terminus_to_path_store(details::thread_id1{}, payload);
}

/// CUDA kernel for running @c traccc::device::gbts_fill_path_store
__global__ void gbts_fill_path_store(
    const device::gbts_fill_path_store_payload payload) {
  __shared__ traccc::uint2
      live_paths[traccc::device::gbts_consts::live_path_buffer];
  __shared__ int n_live_paths;
  const traccc::cuda::barrier barrier;

  device::gbts_fill_path_store(
      details::thread_id1{}, barrier, payload,
      {vecmem::data::vector_view<traccc::uint2>(
           traccc::device::gbts_consts::live_path_buffer, live_paths),
       n_live_paths});
}

/// CUDA kernel for running @c traccc::device::gbts_fit_segments
__global__ void gbts_fit_segments(
    const device::gbts_fit_segments_payload payload) {
  device::gbts_fit_segments(details::thread_id1{}, payload);
}

/// CUDA kernel for running @c traccc::device::gbts_reset_edge_bids
__global__ void gbts_reset_edge_bids(
    const device::gbts_reset_edge_bids_payload payload) {
  device::gbts_reset_edge_bids(details::thread_id1{}, payload);
}

/// CUDA kernel for running @c traccc::device::gbts_rebid_seeds_for_edges
__global__ void gbts_rebid_seeds_for_edges(
    const device::gbts_rebid_seeds_for_edges_payload payload) {
  device::gbts_rebid_seeds_for_edges(details::thread_id1{}, payload);
}

/// CUDA kernel for running @c traccc::device::gbts_bid_seeds_for_hits
__global__ void gbts_bid_seeds_for_hits(
    const device::gbts_bid_seeds_for_hits_payload payload) {
  device::gbts_bid_seeds_for_hits(details::thread_id1{}, payload);
}

/// CUDA kernel for running @c traccc::device::gbts_convert_seeds
__global__ void gbts_convert_seeds(
    const device::gbts_convert_seeds_payload payload) {
  device::gbts_convert_seeds(details::thread_id1{}, payload);
}

}  // namespace kernels

// ===========================================================================
// gbts_seeding_algorithm: kernel launchers
// ===========================================================================

gbts_seeding_algorithm::gbts_seeding_algorithm(
    const gbts_seedfinder_config& cfg, const memory_resource& mr,
    const vecmem::copy& copy, const stream_wrapper& str,
    std::unique_ptr<const Logger> logger)
    : device::gbts_seeding_algorithm(cfg, mr, copy, std::move(logger)),
      cuda::algorithm_base{str} {}

void gbts_seeding_algorithm::gbts_count_spacepoints_by_layer_kernel(
    const device::gbts_count_spacepoints_by_layer_payload& payload) const {
  const unsigned int n_threads = 128;
  const unsigned int n_blocks = 1 + (payload.nSp - 1) / n_threads;
  kernels::gbts_count_spacepoints_by_layer<<<n_blocks, n_threads, 0,
                                             details::get_stream(stream())>>>(
      payload);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());  //
  vecmem::device_vector<unsigned int> d_layer_sums(payload.layerCounts);
  thrust::inclusive_scan(
      thrust::cuda::par_nosync(std::pmr::polymorphic_allocator(&(mr().main)))
          .on(details::get_stream(stream())),
      d_layer_sums.begin(), d_layer_sums.end(), d_layer_sums.begin());
}

void gbts_seeding_algorithm::gbts_bin_spacepoints_kernel(
    const device::gbts_bin_spacepoints_payload& payload) const {
  const unsigned int n_threads = 128;
  const unsigned int n_blocks = 1 + (payload.nSp - 1) / n_threads;
  kernels::gbts_bin_spacepoints<<<n_blocks, n_threads, 0,
                                  details::get_stream(stream())>>>(payload);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());  //
}

void gbts_seeding_algorithm::gbts_count_eta_phi_bins_kernel(
    const device::gbts_count_eta_phi_bins_payload& payload) const {
  const unsigned int n_threads = 128;
  const unsigned int n_blocks = 1 + (payload.nEtaBins - 1) / n_threads;
  kernels::gbts_count_eta_phi_bins<<<n_blocks, n_threads, 0,
                                     details::get_stream(stream())>>>(payload);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());  //
  vecmem::device_vector<unsigned int> d_eta_sums(payload.eta_node_counter);
  thrust::inclusive_scan(
      thrust::cuda::par_nosync(std::pmr::polymorphic_allocator(&(mr().main)))
          .on(details::get_stream(stream())),
      d_eta_sums.begin(), d_eta_sums.end(), d_eta_sums.begin());
}

void gbts_seeding_algorithm::gbts_prefix_sum_eta_phi_bins_kernel(
    const device::gbts_prefix_sum_eta_phi_bins_payload& payload) const {
  const unsigned int n_threads = 128;
  const unsigned int n_blocks = 1 + (payload.nEtaBins - 1) / n_threads;
  kernels::gbts_prefix_sum_eta_phi_bins<<<n_blocks, n_threads, 0,
                                          details::get_stream(stream())>>>(
      payload);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());  //
}

void gbts_seeding_algorithm::gbts_sort_nodes_kernel(
    const device::gbts_sort_nodes_payload& payload) const {
  const unsigned int n_threads = 256;
  const unsigned int n_blocks = 1 + (payload.nNodes - 1) / n_threads;
  kernels::gbts_sort_nodes<<<n_blocks, n_threads, 0,
                             details::get_stream(stream())>>>(payload);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());  //
}

void gbts_seeding_algorithm::gbts_find_minmax_radius_kernel(
    const device::gbts_find_minmax_radius_payload& payload) const {
  const unsigned int n_threads = 128;
  const unsigned int n_blocks = 1 + (payload.nEtaBins - 1) / n_threads;
  kernels::gbts_find_minmax_radius<<<n_blocks, n_threads, 0,
                                     details::get_stream(stream())>>>(payload);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());  //
}

void gbts_seeding_algorithm::gbts_make_graph_edges_kernel(
    const device::gbts_make_graph_edges_payload& payload) const {
  const unsigned int n_threads = 128;
  const unsigned int n_blocks = payload.nUsedBinPairs;
  kernels::gbts_make_graph_edges<<<n_blocks, n_threads, 0,
                                   details::get_stream(stream())>>>(payload);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());  //
  vecmem::device_vector<unsigned int> d_num_outgoing_edges(
      payload.num_outgoing_edges);
  thrust::inclusive_scan(
      thrust::cuda::par_nosync(std::pmr::polymorphic_allocator(&(mr().main)))
          .on(details::get_stream(stream())),
      d_num_outgoing_edges.begin(), d_num_outgoing_edges.end(),
      d_num_outgoing_edges.begin());
}

void gbts_seeding_algorithm::gbts_link_graph_edges_kernel(
    const device::gbts_link_graph_edges_payload& payload) const {
  const unsigned int n_threads = 256;
  const unsigned int n_blocks = 1 + (payload.nEdges - 1) / n_threads;
  kernels::gbts_link_graph_edges<<<n_blocks, n_threads, 0,
                                   details::get_stream(stream())>>>(payload);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());  //
}

void gbts_seeding_algorithm::gbts_match_graph_edges_kernel(
    const device::gbts_match_graph_edges_payload& payload) const {
  const unsigned int n_threads = 256;
  const unsigned int n_blocks = 1 + (payload.nEdges - 1) / n_threads;
  kernels::gbts_match_graph_edges<<<n_blocks, n_threads, 0,
                                    details::get_stream(stream())>>>(payload);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());  //
}

void gbts_seeding_algorithm::gbts_reindex_edges_kernel(
    const device::gbts_reindex_edges_payload& payload) const {
  const unsigned int n_threads = 256;
  const unsigned int n_blocks = 1 + (payload.nEdges - 1) / n_threads;
  kernels::gbts_reindex_edges<<<n_blocks, n_threads, 0,
                                details::get_stream(stream())>>>(payload);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());  //
}

void gbts_seeding_algorithm::gbts_compress_graph_kernel(
    const device::gbts_compress_graph_payload& payload) const {
  const unsigned int n_threads = 256;
  const unsigned int n_blocks = 1 + (payload.nEdges - 1) / n_threads;
  kernels::gbts_compress_graph<<<n_blocks, n_threads, 0,
                                 details::get_stream(stream())>>>(payload);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());  //
}

void gbts_seeding_algorithm::gbts_run_cca_iteration_kernel(
    const device::gbts_run_cca_iteration_payload& payload) const {
  const unsigned int n_threads = 128;
  const unsigned int n_blocks = 1 + (payload.nConnectedEdges - 1) / n_threads;

  kernels::gbts_run_cca_iteration<<<n_blocks, n_threads, 0,
                                    details::get_stream(stream())>>>(payload);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());  //
}

void gbts_seeding_algorithm::gbts_count_terminus_edges_kernel(
    const device::gbts_count_terminus_edges_payload& payload) const {
  const unsigned int n_threads = 128;
  const unsigned int n_blocks = 1 + (payload.nConnectedEdges - 1) / n_threads;
  kernels::gbts_count_terminus_edges<<<n_blocks, n_threads, 0,
                                       details::get_stream(stream())>>>(
      payload);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());  //
}

void gbts_seeding_algorithm::gbts_add_terminus_to_path_store_kernel(
    const device::gbts_add_terminus_to_path_store_payload& payload) const {
  const unsigned int n_threads = 128;
  const unsigned int n_blocks = 1 + (payload.nConnectedEdges - 1) / n_threads;
  kernels::gbts_add_terminus_to_path_store<<<n_blocks, n_threads, 0,
                                             details::get_stream(stream())>>>(
      payload);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
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
  kernels::gbts_fill_path_store<<<n_blocks, n_threads, 0,
                                  details::get_stream(stream())>>>(payload);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void gbts_seeding_algorithm::gbts_fit_segments_kernel(
    const device::gbts_fit_segments_payload& payload) const {
  const unsigned int n_threads = 128;
  const unsigned int n_blocks = 1 + (payload.nPaths - 1) / n_threads;
  kernels::gbts_fit_segments<<<n_blocks, n_threads, 0,
                               details::get_stream(stream())>>>(payload);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());  //
}

void gbts_seeding_algorithm::gbts_reset_edge_bids_kernel(
    const device::gbts_reset_edge_bids_payload& payload) const {
  const unsigned int n_threads = 128;
  const unsigned int n_blocks = 1 + (payload.nProps - 1) / n_threads;
  kernels::gbts_reset_edge_bids<<<n_blocks, n_threads, 0,
                                  details::get_stream(stream())>>>(payload);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void gbts_seeding_algorithm::gbts_rebid_seeds_for_edges_kernel(
    const device::gbts_rebid_seeds_for_edges_payload& payload) const {
  const unsigned int n_threads = 128;
  const unsigned int n_blocks = 1 + (payload.nProps - 1) / n_threads;
  kernels::gbts_rebid_seeds_for_edges<<<n_blocks, n_threads, 0,
                                        details::get_stream(stream())>>>(
      payload);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void gbts_seeding_algorithm::gbts_bid_seeds_for_hits_kernel(
    const device::gbts_bid_seeds_for_hits_payload& payload) const {
  const unsigned int n_threads = 128;
  const unsigned int n_blocks = 1 + (payload.nProps - 1) / n_threads;
  kernels::gbts_bid_seeds_for_hits<<<n_blocks, n_threads, 0,
                                     details::get_stream(stream())>>>(payload);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void gbts_seeding_algorithm::gbts_convert_seeds_kernel(
    const device::gbts_convert_seeds_payload& payload) const {
  const unsigned int n_threads = 128;
  const unsigned int n_blocks = 1 + (payload.nProps - 1) / n_threads;
  kernels::gbts_convert_seeds<<<n_blocks, n_threads, 0,
                                details::get_stream(stream())>>>(payload);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());  //
}

}  // namespace traccc::cuda
