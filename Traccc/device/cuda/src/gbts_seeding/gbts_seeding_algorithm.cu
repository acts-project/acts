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

// System include(s).
#include <algorithm>
#include <memory_resource>

// Thrust include(s).
#include <thrust/execution_policy.h>
#include <thrust/scan.h>
#include <thrust/sort.h>

namespace traccc::cuda {

namespace kernels {

using float4 = traccc::float4;
using uint2 = traccc::uint2;
using int2 = traccc::int2;

// ---------------------------------------------------------------------------
// Stage 1 — nodes-making kernels
// ---------------------------------------------------------------------------

/// CUDA kernel for running @c traccc::device::gbts_bin_spacepoints
__global__ void gbts_bin_spacepoints(
    const device::gbts_bin_spacepoints_payload payload) {
  device::gbts_bin_spacepoints(details::thread_id1{}, payload);
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

/// CUDA kernel for running @c traccc::device::gbts_build_edge_work_list
__global__ void gbts_build_edge_work_list(
    const device::gbts_build_edge_work_list_payload payload) {
  __shared__ unsigned int scratch[device::gbts_build_edge_work_list_block_size];
  const traccc::cuda::barrier barrier;
  device::gbts_build_edge_work_list(
      details::thread_id1{}, barrier, payload,
      {vecmem::data::vector_view<unsigned int>(
          device::gbts_build_edge_work_list_block_size, scratch)});
}

/// CUDA kernel for running @c traccc::device::gbts_count_graph_edges
__global__ void gbts_count_graph_edges(
    const device::gbts_count_graph_edges_payload payload) {
  device::gbts_count_graph_edges(details::thread_id1{}, payload);
}

/// CUDA kernel for running @c traccc::device::gbts_fill_graph_edges
__global__ void gbts_fill_graph_edges(
    const device::gbts_fill_graph_edges_payload payload) {
  device::gbts_fill_graph_edges(details::thread_id1{}, payload);
}

/// CUDA kernel for running @c traccc::device::gbts_match_graph_edges
__global__ void gbts_match_graph_edges(
    const device::gbts_match_graph_edges_payload payload) {
  device::gbts_match_graph_edges(details::thread_id1{}, payload);
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

/// CUDA kernel for running @c traccc::device::gbts_finish_cca
__global__ void gbts_finish_cca(const device::gbts_finish_cca_payload payload) {
  device::gbts_finish_cca(details::thread_id1{}, payload);
}

/// CUDA kernel for running @c traccc::device::gbts_count_paths
__global__ void gbts_count_paths(
    const device::gbts_count_paths_payload payload) {
  device::gbts_count_paths(details::thread_id1{}, payload);
}

/// CUDA kernel for running @c traccc::device::gbts_fill_path_store
__global__ void gbts_fill_path_store(
    const device::gbts_fill_path_store_payload payload) {
  device::gbts_fill_path_store(details::thread_id1{}, payload);
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

void gbts_seeding_algorithm::gbts_bin_spacepoints_kernel(
    const device::gbts_bin_spacepoints_payload& payload) const {
  const unsigned int n_threads = 128;
  const unsigned int n_blocks = 1 + (payload.nSp - 1) / n_threads;
  kernels::gbts_bin_spacepoints<<<n_blocks, n_threads, 0,
                                  details::get_stream(stream())>>>(payload);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());

  // Turn the per-bin node counts into the node offsets.
  vecmem::device_vector<unsigned int> d_eta_node_counter(
      payload.eta_node_counter);
  thrust::exclusive_scan(
      thrust::cuda::par_nosync(std::pmr::polymorphic_allocator(&(mr().main)))
          .on(details::get_stream(stream())),
      d_eta_node_counter.begin(), d_eta_node_counter.end(),
      d_eta_node_counter.begin());
}

void gbts_seeding_algorithm::gbts_sort_nodes_kernel(
    const device::gbts_sort_nodes_payload& payload) const {
  // Order the nodes by their (eta bin, phi, spacepoint index bits) keys,
  // carrying the full spacepoint index along as the value.
  vecmem::device_vector<unsigned long long int> d_sort_keys(payload.sort_keys);
  vecmem::device_vector<unsigned int> d_sort_values(payload.sort_values);
  thrust::sort_by_key(
      thrust::cuda::par_nosync(std::pmr::polymorphic_allocator(&(mr().main)))
          .on(details::get_stream(stream())),
      d_sort_keys.begin(),
      d_sort_keys.begin() + static_cast<int>(payload.nNodes),
      d_sort_values.begin());

  const unsigned int n_threads = 256;
  const unsigned int n_blocks = 1 + (payload.nNodes - 1) / n_threads;
  kernels::gbts_sort_nodes<<<n_blocks, n_threads, 0,
                             details::get_stream(stream())>>>(payload);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void gbts_seeding_algorithm::gbts_find_minmax_radius_kernel(
    const device::gbts_find_minmax_radius_payload& payload) const {
  const unsigned int n_threads = 128;
  const unsigned int n_blocks = 1 + (payload.nEtaBins - 1) / n_threads;
  kernels::gbts_find_minmax_radius<<<n_blocks, n_threads, 0,
                                     details::get_stream(stream())>>>(payload);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());  //
}

void gbts_seeding_algorithm::gbts_build_edge_work_list_kernel(
    const device::gbts_build_edge_work_list_payload& payload) const {
  kernels::gbts_build_edge_work_list<<<
      1, device::gbts_build_edge_work_list_block_size, 0,
      details::get_stream(stream())>>>(payload);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void gbts_seeding_algorithm::gbts_count_graph_edges_kernel(
    const device::gbts_count_graph_edges_payload& payload) const {
  // One thread per inner node of a chunk. The blocks stride over the work
  // items.
  const unsigned int n_threads = device::gbts_consts::edge_chunk_size;
  const unsigned int n_blocks =
      std::min(payload.nWorkMax, device::gbts_count_graph_edges_max_blocks);
  kernels::gbts_count_graph_edges<<<n_blocks, n_threads, 0,
                                    details::get_stream(stream())>>>(payload);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());

  // Turn the per-node counts into the edge buckets.
  vecmem::device_vector<unsigned int> d_num_outgoing_edges(
      payload.num_outgoing_edges);
  thrust::inclusive_scan(
      thrust::cuda::par_nosync(std::pmr::polymorphic_allocator(&(mr().main)))
          .on(details::get_stream(stream())),
      d_num_outgoing_edges.begin(), d_num_outgoing_edges.end(),
      d_num_outgoing_edges.begin());
}

void gbts_seeding_algorithm::gbts_fill_graph_edges_kernel(
    const device::gbts_fill_graph_edges_payload& payload) const {
  const unsigned int n_threads = device::gbts_consts::edge_chunk_size;
  const unsigned int n_blocks =
      std::min(payload.nWorkMax, device::gbts_fill_graph_edges_max_blocks);
  kernels::gbts_fill_graph_edges<<<n_blocks, n_threads, 0,
                                   details::get_stream(stream())>>>(payload);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void gbts_seeding_algorithm::gbts_match_graph_edges_kernel(
    const device::gbts_match_graph_edges_payload& payload) const {
  const unsigned int n_threads = 256;
  const unsigned int n_blocks = 1u + (payload.nEdgesMax - 1u) / n_threads;
  kernels::gbts_match_graph_edges<<<n_blocks, n_threads, 0,
                                    details::get_stream(stream())>>>(payload);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());

  // Compact the kept edges with a prefix sum over their 0/1 flags.
  thrust::inclusive_scan(
      thrust::cuda::par_nosync(std::pmr::polymorphic_allocator(&(mr().main)))
          .on(details::get_stream(stream())),
      payload.kept.ptr(), payload.kept.ptr() + payload.nEdgesMax,
      payload.reIndexer.ptr());
}

void gbts_seeding_algorithm::gbts_compress_graph_kernel(
    const device::gbts_compress_graph_payload& payload) const {
  const unsigned int n_threads = 256;
  const unsigned int n_blocks = 1u + (payload.nEdgesMax - 1u) / n_threads;
  kernels::gbts_compress_graph<<<n_blocks, n_threads, 0,
                                 details::get_stream(stream())>>>(payload);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void gbts_seeding_algorithm::gbts_run_cca_iteration_kernel(
    const device::gbts_run_cca_iteration_payload& payload) const {
  const unsigned int n_threads = 256;
  const unsigned int n_blocks = 1 + (payload.nConnectedEdges - 1) / n_threads;

  kernels::gbts_run_cca_iteration<<<n_blocks, n_threads, 0,
                                    details::get_stream(stream())>>>(payload);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void gbts_seeding_algorithm::gbts_finish_cca_kernel(
    const device::gbts_finish_cca_payload& payload) const {
  const unsigned int n_threads = 256;
  const unsigned int n_blocks = 1 + (payload.nConnectedEdges - 1) / n_threads;
  kernels::gbts_finish_cca<<<n_blocks, n_threads, 0,
                             details::get_stream(stream())>>>(payload);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void gbts_seeding_algorithm::gbts_count_paths_kernel(
    const device::gbts_count_paths_payload& payload) const {
  const unsigned int n_threads = 128;
  const unsigned int n_blocks = 1 + (payload.nConnectedEdges - 1) / n_threads;
  kernels::gbts_count_paths<<<n_blocks, n_threads, 0,
                              details::get_stream(stream())>>>(payload);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
  // Path offsets of the path store.
  vecmem::device_vector<unsigned int> d_path_counts(payload.path_counts);
  thrust::inclusive_scan(
      thrust::cuda::par_nosync(std::pmr::polymorphic_allocator(&(mr().main)))
          .on(details::get_stream(stream())),
      d_path_counts.begin(), d_path_counts.begin() + payload.nConnectedEdges,
      d_path_counts.begin());
}

void gbts_seeding_algorithm::gbts_fill_path_store_kernel(
    const device::gbts_fill_path_store_payload& payload) const {
  const unsigned int n_threads = 128;
  const unsigned int n_blocks = 1 + (payload.nPathsGrid - 1) / n_threads;
  kernels::gbts_fill_path_store<<<n_blocks, n_threads, 0,
                                  details::get_stream(stream())>>>(payload);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void gbts_seeding_algorithm::gbts_bid_seeds_for_hits_kernel(
    const device::gbts_bid_seeds_for_hits_payload& payload) const {
  const unsigned int n_threads = 128;
  const unsigned int n_blocks = 1 + (payload.nPathsGrid - 1) / n_threads;
  kernels::gbts_bid_seeds_for_hits<<<n_blocks, n_threads, 0,
                                     details::get_stream(stream())>>>(payload);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void gbts_seeding_algorithm::gbts_convert_seeds_kernel(
    const device::gbts_convert_seeds_payload& payload) const {
  const unsigned int n_threads = 128;
  const unsigned int n_blocks = 1 + (payload.nPathsGrid - 1) / n_threads;
  kernels::gbts_convert_seeds<<<n_blocks, n_threads, 0,
                                details::get_stream(stream())>>>(payload);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void gbts_seeding_algorithm::synchronize() const {
  stream().synchronize();
}

}  // namespace traccc::cuda
