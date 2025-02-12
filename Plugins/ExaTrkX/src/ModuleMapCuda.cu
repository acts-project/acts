// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/ModuleMapCuda.hpp"
#include "Acts/Plugins/ExaTrkX/detail/CudaUtils.cuh"
#include "Acts/Plugins/ExaTrkX/detail/ModuleMapUtils.cuh"
#include "Acts/Plugins/ExaTrkX/detail/TensorVectorConversion.hpp"

#include <CUDA_graph_creator>
#include <CUDA_module_map_doublet>
#include <CUDA_module_map_triplet>
#include <TTree_hits>
#include <chrono>

#include <c10/cuda/CUDAGuard.h>
#include <thrust/execution_policy.h>
#include <thrust/functional.h>
#include <thrust/scan.h>
#include <thrust/sort.h>
#include <thrust/transform_scan.h>

using Clock = std::chrono::high_resolution_clock;

using namespace torch::indexing;

namespace Acts {

ModuleMapCuda::ModuleMapCuda(const Config &cfg,
                             std::unique_ptr<const Acts::Logger> logger_)
    : m_logger(std::move(logger_)), m_cfg(cfg) {
  module_map_triplet<float> moduleMapCpu;
  moduleMapCpu.read_TTree(cfg.moduleMapPath.c_str());
  if (!moduleMapCpu) {
    throw std::runtime_error("Cannot retrieve ModuleMap from " +
                             cfg.moduleMapPath);
  }

  ACTS_DEBUG("ModuleMap GPU block dim: " << m_cfg.gpuBlocks);

  m_cudaModuleMapDoublet =
      std::make_unique<CUDA_module_map_doublet<float>>(moduleMapCpu);
  m_cudaModuleMapDoublet->HostToDevice();
  m_cudaModuleMapTriplet =
      std::make_unique<CUDA_module_map_triplet<float>>(moduleMapCpu);
  m_cudaModuleMapTriplet->HostToDevice();

  ACTS_DEBUG("# of modules = " << moduleMapCpu.module_map().size());

  // check if we actually have a module map
  std::map<std::uint64_t, int> test;

  std::vector<std::uint64_t> keys;
  keys.reserve(m_cudaModuleMapDoublet->module_map().size());
  std::vector<int> vals;
  vals.reserve(m_cudaModuleMapDoublet->module_map().size());

  for (auto [key, value] : m_cudaModuleMapDoublet->module_map()) {
    auto [it, success] = test.insert({key, value});
    if (!success) {
      throw std::runtime_error("Duplicate key in module map");
    }
    keys.push_back(key);
    vals.push_back(value);
  }

  // copy module map to device
  m_cudaModuleMapSize = m_cudaModuleMapDoublet->module_map().size();
  cudaMalloc(&m_cudaModuleMapKeys, m_cudaModuleMapSize * sizeof(std::uint64_t));
  cudaMalloc(&m_cudaModuleMapVals, m_cudaModuleMapSize * sizeof(int));

  cudaMemcpy(m_cudaModuleMapKeys, keys.data(),
             m_cudaModuleMapSize * sizeof(std::uint64_t),
             cudaMemcpyHostToDevice);
  cudaMemcpy(m_cudaModuleMapVals, vals.data(),
             m_cudaModuleMapSize * sizeof(int), cudaMemcpyHostToDevice);
}

ModuleMapCuda::~ModuleMapCuda() {
  cudaFree(m_cudaModuleMapKeys);
  cudaFree(m_cudaModuleMapVals);
}

namespace {}  // namespace

std::tuple<std::any, std::any, std::any> ModuleMapCuda::operator()(
    std::vector<float> &inputValues, std::size_t numNodes,
    const std::vector<std::uint64_t> &moduleIds,
    const ExecutionContext &execContext) {
  const auto nHits = moduleIds.size();
  const auto nFeatures = inputValues.size() / moduleIds.size();
  auto &features = inputValues;

  const dim3 blockDim = m_cfg.gpuBlocks;
  const dim3 gridDimHits = (nHits + blockDim.x - 1) / blockDim.x;
  ACTS_VERBOSE("gridDimHits: " << gridDimHits.x
                               << ", blockDim: " << blockDim.x);

  // Get stream if available, otherwise use default stream
  cudaStream_t stream;
#if 0
  ACTS_CUDA_CHECK(cudaStreamCreate(&stream));
#else
  assert(execContext.stream);
  std::optional<c10::cuda::CUDAStreamGuard> streamGuard;
  if (execContext.stream) {
    ACTS_VERBOSE("Got stream " << *execContext.stream);
    stream = execContext.stream->stream();
    streamGuard.emplace(*execContext.stream);
  }
#endif

  /////////////////////////
  // Prepare input data
  ////////////////////////

  // Full node features to device
  float *cudaNodeFeatures{};
  ACTS_CUDA_CHECK(cudaMallocAsync(&cudaNodeFeatures,
                                  features.size() * sizeof(float), stream));
  ACTS_CUDA_CHECK(cudaMemcpyAsync(cudaNodeFeatures, features.data(),
                                  features.size() * sizeof(float),
                                  cudaMemcpyHostToDevice, stream));

  // Module IDs to device
  std::uint64_t *cudaModuleIds;
  ACTS_CUDA_CHECK(
      cudaMallocAsync(&cudaModuleIds, nHits * sizeof(std::uint64_t), stream));
  ACTS_CUDA_CHECK(cudaMemcpyAsync(cudaModuleIds, moduleIds.data(),
                                  nHits * sizeof(std::uint64_t),
                                  cudaMemcpyHostToDevice, stream));

  // Node features for module map graph
  constexpr std::size_t rOffset = 0;
  constexpr std::size_t phiOffset = 1;
  constexpr std::size_t zOffset = 2;
  constexpr std::size_t etaOffset = 3;

  const auto srcStride = sizeof(float) * nFeatures;
  const auto dstStride = sizeof(float);  // contiguous in destination
  const auto width = sizeof(float);      // only copy 1 column
  const auto height = nHits;

  Acts::detail::CUDA_hit_data<float> inputData;
  inputData.m_size = nHits;

  ACTS_CUDA_CHECK(
      cudaMallocAsync(&inputData.m_cuda_x, nHits * sizeof(float), stream));

  ACTS_CUDA_CHECK(
      cudaMallocAsync(&inputData.m_cuda_y, nHits * sizeof(float), stream));

  ACTS_CUDA_CHECK(cudaMallocAsync(&inputData.m_cuda_hit_id,
                                  nHits * sizeof(std::uint64_t), stream));

  ACTS_CUDA_CHECK(
      cudaMallocAsync(&inputData.m_cuda_R, nHits * sizeof(float), stream));
  ACTS_CUDA_CHECK(cudaMemcpy2DAsync(
      inputData.m_cuda_R, dstStride, cudaNodeFeatures + rOffset, srcStride,
      width, height, cudaMemcpyDeviceToDevice, stream));

  ACTS_CUDA_CHECK(
      cudaMallocAsync(&inputData.m_cuda_phi, nHits * sizeof(float), stream));
  ACTS_CUDA_CHECK(cudaMemcpy2DAsync(
      inputData.m_cuda_phi, dstStride, cudaNodeFeatures + phiOffset, srcStride,
      width, height, cudaMemcpyDeviceToDevice, stream));

  ACTS_CUDA_CHECK(
      cudaMallocAsync(&inputData.m_cuda_z, nHits * sizeof(float), stream));
  ACTS_CUDA_CHECK(cudaMemcpy2DAsync(
      inputData.m_cuda_z, dstStride, cudaNodeFeatures + zOffset, srcStride,
      width, height, cudaMemcpyDeviceToDevice, stream));

  ACTS_CUDA_CHECK(
      cudaMallocAsync(&inputData.m_cuda_eta, nHits * sizeof(float), stream));
  ACTS_CUDA_CHECK(cudaMemcpy2DAsync(
      inputData.m_cuda_eta, dstStride, cudaNodeFeatures + etaOffset, srcStride,
      width, height, cudaMemcpyDeviceToDevice, stream));

  // Allocate helper nb hits memory
  int *cudaNbHits;
  ACTS_CUDA_CHECK(cudaMallocAsync(
      &cudaNbHits, (m_cudaModuleMapSize + 1) * sizeof(int), stream));
  ACTS_CUDA_CHECK(cudaMemsetAsync(
      cudaNbHits, 0, (m_cudaModuleMapSize + 1) * sizeof(int), stream));

  // Synchronize to ensure memory is there (TODO is this not garantueed?)
  ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));

  // Preprocess features
  detail::rescaleFeature<<<gridDimHits, blockDim, 0, stream>>>(
      nHits, inputData.m_cuda_z, 1000.f);
  ACTS_CUDA_CHECK(cudaGetLastError());
  detail::rescaleFeature<<<gridDimHits, blockDim, 0, stream>>>(
      nHits, inputData.m_cuda_R, 1000.f);
  ACTS_CUDA_CHECK(cudaGetLastError());
  detail::rescaleFeature<<<gridDimHits, blockDim, 0, stream>>>(
      nHits, inputData.m_cuda_phi, 3.14159f);
  ACTS_CUDA_CHECK(cudaGetLastError());

  detail::computeXandY<<<gridDimHits, blockDim, 0, stream>>>(
      nHits, inputData.m_cuda_x, inputData.m_cuda_y, inputData.m_cuda_R,
      inputData.m_cuda_phi);
  ACTS_CUDA_CHECK(cudaGetLastError());

  detail::setHitId<<<gridDimHits, blockDim, 0, stream>>>(
      nHits, inputData.m_cuda_hit_id);
  ACTS_CUDA_CHECK(cudaGetLastError());

  detail::mapModuleIdsToNbHits<<<gridDimHits, blockDim, 0, stream>>>(
      cudaNbHits, nHits, cudaModuleIds, m_cudaModuleMapSize,
      m_cudaModuleMapKeys, m_cudaModuleMapVals);
  ACTS_CUDA_CHECK(cudaGetLastError());

  thrust::exclusive_scan(thrust::device.on(stream), cudaNbHits,
                         cudaNbHits + m_cudaModuleMapSize + 1, cudaNbHits);
  ACTS_CUDA_CHECK(cudaGetLastError());
  int *cudaHitIndice = cudaNbHits;

  ///////////////////////////////////
  // Perform module map inference
  ////////////////////////////////////
  const auto edgeData = makeEdges(inputData, cudaHitIndice, stream);
  ACTS_CUDA_CHECK(cudaGetLastError());

  ACTS_CUDA_CHECK(cudaFreeAsync(inputData.cuda_R(), stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(inputData.cuda_phi(), stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(inputData.cuda_z(), stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(inputData.cuda_eta(), stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(inputData.cuda_x(), stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(inputData.cuda_y(), stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(inputData.cuda_hit_id(), stream));

  dim3 gridDimEdges = (edgeData.nEdges + blockDim.x - 1) / blockDim.x;

  // Make edge features
  float *edgeFeaturePtr{};
  ACTS_CUDA_CHECK(cudaMallocAsync(&edgeFeaturePtr,
                                  6 * edgeData.nEdges * sizeof(float), stream));

  detail::makeEdgeFeatures<<<gridDimEdges, blockDim, 0, stream>>>(
      edgeData.nEdges, edgeData.cudaEdgePtr,
      edgeData.cudaEdgePtr + edgeData.nEdges, nFeatures, cudaNodeFeatures,
      edgeFeaturePtr);
  ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));
  ACTS_CUDA_CHECK(cudaGetLastError());

  //  Make torch tensors
  auto edgeFeatures = torch::from_blob(
      edgeFeaturePtr, 6 * static_cast<long>(edgeData.nEdges),
      [stream](void *ptr) { ACTS_CUDA_CHECK(cudaFreeAsync(ptr, stream)); },
      at::TensorOptions().device(at::kCUDA).dtype(at::kFloat));
  edgeFeatures = edgeFeatures.reshape({static_cast<long>(edgeData.nEdges), 6});
  ACTS_CUDA_CHECK(cudaGetLastError());

  ACTS_VERBOSE("edge features:\n"
               << edgeFeatures.index({Slice(None, 5), Slice()}));
  auto edgeFeaturesNew = edgeFeatures.clone();

  auto edgeIndex =
      torch::from_blob(edgeData.cudaEdgePtr,
                       {2, static_cast<long>(edgeData.nEdges)},
                       at::TensorOptions().device(at::kCUDA).dtype(at::kInt))
          .to(torch::kLong);
  ACTS_CUDA_CHECK(cudaFreeAsync(edgeData.cudaEdgePtr, stream));
  ACTS_VERBOSE("edge index:\n" << edgeIndex.index({Slice(), Slice(0, 10)}));

  auto nodeFeatures = torch::from_blob(
      cudaNodeFeatures, inputValues.size(),
      [stream](void *ptr) { ACTS_CUDA_CHECK(cudaFreeAsync(ptr, stream)); },
      at::TensorOptions().device(at::kCUDA).dtype(at::kFloat));
  nodeFeatures = nodeFeatures.reshape({(long)nHits, (long)nFeatures});
  ACTS_CUDA_CHECK(cudaGetLastError());

  return {nodeFeatures, edgeIndex, edgeFeatures};
}

struct ArgsortFun {
  int *cudaPtr = nullptr;
  bool __device__ operator()(int l, int r) { return cudaPtr[l] < cudaPtr[r]; }
};

struct CastBoolToInt {
  int __device__ operator()(bool b) { return static_cast<int>(b); }
};

detail::CUDA_edge_data<float> ModuleMapCuda::makeEdges(
    detail::CUDA_hit_data<float> cuda_TThits, int *cuda_hit_indice,
    cudaStream_t &stream) const {
  const dim3 block_dim = m_cfg.gpuBlocks;
  // ----------------------------------
  // memory allocation for hits + edges
  // ----------------------------------
  const int nb_doublets = m_cudaModuleMapDoublet->size();
  ACTS_DEBUG("nb doublets " << nb_doublets);
  dim3 grid_dim = ((nb_doublets + block_dim.x - 1) / block_dim.x);

  // hit buffer pointers need to be visible outside the scope
  int *cuda_reduced_M1_hits, *cuda_reduced_M2_hits;

  // count the edges per doublet
  int *cuda_nb_edges_per_doublet;
  ACTS_CUDA_CHECK(cudaMallocAsync(&cuda_nb_edges_per_doublet,
                                  (nb_doublets + 1) * sizeof(int), stream));

  int *cuda_edge_sum;
  ACTS_CUDA_CHECK(
      cudaMallocAsync(&cuda_edge_sum, (nb_doublets + 1) * sizeof(int), stream));

  int nb_doublet_edges{};
  // ---------------------------------------------
  // A: Allocate fixed memory for edges
  // ---------------------------------------------
  if (m_cfg.maxEdgesAllocate > 0) {
    throw std::runtime_error("Not implemented");
  }
  // ---------------------------------------------
  // B: Allocate memory for edges dynamically
  // ---------------------------------------------
  else {
    ACTS_DEBUG("Allocate memory for edges dynamically");
    // Allocate one integer on device and set it to 0
    int *cuda_nb_doublet_edges;
    ACTS_CUDA_CHECK(
        cudaMallocAsync(&cuda_nb_doublet_edges, sizeof(int), stream));
    ACTS_CUDA_CHECK(
        cudaMemsetAsync(cuda_nb_doublet_edges, 0, sizeof(int), stream));

    detail::count_doublet_edges<float><<<grid_dim, block_dim, 0, stream>>>(
        nb_doublets, m_cudaModuleMapDoublet->cuda_module1(),
        m_cudaModuleMapDoublet->cuda_module2(), cuda_TThits.cuda_R(),
        cuda_TThits.cuda_z(), cuda_TThits.cuda_eta(), cuda_TThits.cuda_phi(),
        m_cudaModuleMapDoublet->cuda_z0_min(),
        m_cudaModuleMapDoublet->cuda_z0_max(),
        m_cudaModuleMapDoublet->cuda_deta_min(),
        m_cudaModuleMapDoublet->cuda_deta_max(),
        m_cudaModuleMapDoublet->cuda_phi_slope_min(),
        m_cudaModuleMapDoublet->cuda_phi_slope_max(),
        m_cudaModuleMapDoublet->cuda_dphi_min(),
        m_cudaModuleMapDoublet->cuda_dphi_max(), cuda_hit_indice, TMath::Pi(),
        cuda_nb_doublet_edges, cuda_edge_sum);
    ACTS_CUDA_CHECK(cudaGetLastError());

    // Copy the number of edges to the host, synchronize and allocate
    ACTS_CUDA_CHECK(cudaMemcpyAsync(&nb_doublet_edges, cuda_nb_doublet_edges,
                                    sizeof(int), cudaMemcpyDeviceToHost,
                                    stream));
    ACTS_CUDA_CHECK(cudaFreeAsync(cuda_nb_doublet_edges, stream));
    ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));
    ACTS_DEBUG("nb_doublet_edges: " << nb_doublet_edges);
    ACTS_DEBUG("Allocate " << (2ul * nb_doublet_edges * sizeof(int)) * 1.0e-6
                           << " MB for edges");
    ACTS_CUDA_CHECK(cudaMallocAsync(
        &cuda_reduced_M1_hits,
        static_cast<std::size_t>(nb_doublet_edges) * sizeof(int), stream));
    ACTS_CUDA_CHECK(cudaMallocAsync(
        &cuda_reduced_M2_hits,
        static_cast<std::size_t>(nb_doublet_edges) * sizeof(int), stream));

    // Prefix sum to get the edge offset for each doublet
    thrust::exclusive_scan(thrust::device.on(stream), cuda_edge_sum,
                           cuda_edge_sum + (nb_doublets + 1), cuda_edge_sum);

    detail::doublet_cuts_new<float><<<grid_dim, block_dim, 0, stream>>>(
        nb_doublets, m_cudaModuleMapDoublet->cuda_module1(),
        m_cudaModuleMapDoublet->cuda_module2(), cuda_TThits.cuda_R(),
        cuda_TThits.cuda_z(), cuda_TThits.cuda_eta(), cuda_TThits.cuda_phi(),
        m_cudaModuleMapDoublet->cuda_z0_min(),
        m_cudaModuleMapDoublet->cuda_z0_max(),
        m_cudaModuleMapDoublet->cuda_deta_min(),
        m_cudaModuleMapDoublet->cuda_deta_max(),
        m_cudaModuleMapDoublet->cuda_phi_slope_min(),
        m_cudaModuleMapDoublet->cuda_phi_slope_max(),
        m_cudaModuleMapDoublet->cuda_dphi_min(),
        m_cudaModuleMapDoublet->cuda_dphi_max(), cuda_hit_indice, TMath::Pi(),
        cuda_reduced_M1_hits, cuda_reduced_M2_hits, cuda_edge_sum);
    ACTS_CUDA_CHECK(cudaGetLastError());
  }

  //---------------------------------------------
  // reduce nb of hits and sort by hid id M2 hits
  //---------------------------------------------

  int *cuda_sorted_M2_hits;
  ACTS_CUDA_CHECK(cudaMallocAsync(&cuda_sorted_M2_hits,
                                  nb_doublet_edges * sizeof(int), stream));
  ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));
  grid_dim = ((nb_doublet_edges + block_dim.x - 1) / block_dim.x);
  init_vector<<<grid_dim, block_dim, 0, stream>>>(cuda_sorted_M2_hits,
                                                  nb_doublet_edges);
  ACTS_CUDA_CHECK(cudaGetLastError());

  grid_dim = ((nb_doublets + block_dim.x - 1) / block_dim.x);
  partial_quick_sort<<<grid_dim, block_dim, 0, stream>>>(
      cuda_sorted_M2_hits, cuda_reduced_M2_hits, cuda_edge_sum, nb_doublets);
  ACTS_CUDA_CHECK(cudaGetLastError());

  // -----------------------------
  // build doublets geometric cuts
  // -----------------------------

  float *cuda_z0, *cuda_phi_slope, *cuda_deta, *cuda_dphi;
  ACTS_CUDA_CHECK(
      cudaMallocAsync(&cuda_z0, nb_doublet_edges * sizeof(float), stream));
  ACTS_CUDA_CHECK(cudaMallocAsync(&cuda_phi_slope,
                                  nb_doublet_edges * sizeof(float), stream));
  ACTS_CUDA_CHECK(
      cudaMallocAsync(&cuda_deta, nb_doublet_edges * sizeof(float), stream));
  ACTS_CUDA_CHECK(
      cudaMallocAsync(&cuda_dphi, nb_doublet_edges * sizeof(float), stream));
  ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));
  grid_dim = ((nb_doublet_edges + block_dim.x - 1) / block_dim.x);
  hits_geometric_cuts<<<grid_dim, block_dim, 0, stream>>>(
      cuda_z0, cuda_phi_slope, cuda_deta, cuda_dphi, cuda_reduced_M1_hits,
      cuda_reduced_M2_hits, cuda_TThits.cuda_R(), cuda_TThits.cuda_z(),
      cuda_TThits.cuda_eta(), cuda_TThits.cuda_phi(), TMath::Pi(),
      nb_doublet_edges);
  ACTS_CUDA_CHECK(cudaGetLastError());
  // free mem at the end

  ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));

  // CUDA_mask cuda_edge_tag(nb_doublet_edges);
  bool *cuda_mask;
  ACTS_CUDA_CHECK(cudaMallocAsync(
      &cuda_mask, (nb_doublet_edges + 1) * sizeof(bool), stream));
  ACTS_CUDA_CHECK(
      cudaMemset(cuda_mask, 0, (nb_doublet_edges + 1) * sizeof(bool)));
  ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));

  // -------------------------
  // loop over module triplets
  // -------------------------
  int nb_triplets = m_cudaModuleMapTriplet->size();
  grid_dim = ((nb_triplets + block_dim.x - 1) / block_dim.x);

  bool *cuda_vertices;
  ACTS_CUDA_CHECK(cudaMallocAsync(&cuda_vertices,
                                  cuda_TThits.size() * sizeof(bool), stream));
  ACTS_CUDA_CHECK(cudaMemsetAsync(cuda_vertices, 0,
                                  cuda_TThits.size() * sizeof(bool), stream));
  ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));

  Acts::detail::triplet_cuts_new<float><<<grid_dim, block_dim, 0, stream>>>(
      nb_triplets, m_cudaModuleMapTriplet->cuda_module12_map(),
      m_cudaModuleMapTriplet->cuda_module23_map(), cuda_TThits.cuda_x(),
      cuda_TThits.cuda_y(), cuda_TThits.cuda_z(), cuda_TThits.cuda_R(), cuda_z0,
      cuda_phi_slope, cuda_deta, cuda_dphi,
      m_cudaModuleMapTriplet->module12().cuda_z0_min(),
      m_cudaModuleMapTriplet->module12().cuda_z0_max(),
      m_cudaModuleMapTriplet->module12().cuda_deta_min(),
      m_cudaModuleMapTriplet->module12().cuda_deta_max(),
      m_cudaModuleMapTriplet->module12().cuda_phi_slope_min(),
      m_cudaModuleMapTriplet->module12().cuda_phi_slope_max(),
      m_cudaModuleMapTriplet->module12().cuda_dphi_min(),
      m_cudaModuleMapTriplet->module12().cuda_dphi_max(),
      m_cudaModuleMapTriplet->module23().cuda_z0_min(),
      m_cudaModuleMapTriplet->module23().cuda_z0_max(),
      m_cudaModuleMapTriplet->module23().cuda_deta_min(),
      m_cudaModuleMapTriplet->module23().cuda_deta_max(),
      m_cudaModuleMapTriplet->module23().cuda_phi_slope_min(),
      m_cudaModuleMapTriplet->module23().cuda_phi_slope_max(),
      m_cudaModuleMapTriplet->module23().cuda_dphi_min(),
      m_cudaModuleMapTriplet->module23().cuda_dphi_max(),
      m_cudaModuleMapTriplet->cuda_diff_dydx_min(),
      m_cudaModuleMapTriplet->cuda_diff_dydx_max(),
      m_cudaModuleMapTriplet->cuda_diff_dzdr_min(),
      m_cudaModuleMapTriplet->cuda_diff_dzdr_max(), TMath::Pi(),
      cuda_reduced_M1_hits, cuda_reduced_M2_hits, cuda_sorted_M2_hits,
      cuda_edge_sum, cuda_vertices, cuda_mask);
  ACTS_CUDA_CHECK(cudaGetLastError());

  //----------------
  // edges reduction
  //----------------
  int *cuda_mask_sum;
  ACTS_CUDA_CHECK(cudaMallocAsync(
      &cuda_mask_sum, (nb_doublet_edges + 1) * sizeof(int), stream));
  ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));
  thrust::transform_exclusive_scan(
      thrust::device.on(stream), cuda_mask, cuda_mask + (nb_doublet_edges + 1),
      cuda_mask_sum, CastBoolToInt{}, 0, thrust::plus<int>());
  int nb_graph_edges;
  ACTS_CUDA_CHECK(cudaMemcpyAsync(&nb_graph_edges,
                                  &cuda_mask_sum[nb_doublet_edges], sizeof(int),
                                  cudaMemcpyDeviceToHost, stream));
  ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));

  ACTS_VERBOSE("nb_graph_edges: " << nb_graph_edges);

  int *cuda_graph_edge_ptr{};
  cudaMallocAsync(&cuda_graph_edge_ptr, 2 * nb_graph_edges * sizeof(int),
                  stream);
  int *cuda_graph_M1_hits = cuda_graph_edge_ptr;
  int *cuda_graph_M2_hits = cuda_graph_edge_ptr + nb_graph_edges;

  grid_dim = ((nb_doublet_edges + block_dim.x - 1) / block_dim.x);
  ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));
  compact_stream<<<grid_dim, block_dim, 0, stream>>>(
      cuda_graph_M1_hits, cuda_reduced_M1_hits, cuda_mask, cuda_mask_sum,
      nb_doublet_edges);
  ACTS_CUDA_CHECK(cudaGetLastError());
  compact_stream<<<grid_dim, block_dim, 0, stream>>>(
      cuda_graph_M2_hits, cuda_reduced_M2_hits, cuda_mask, cuda_mask_sum,
      nb_doublet_edges);
  ACTS_CUDA_CHECK(cudaGetLastError());

  ACTS_CUDA_CHECK(cudaFreeAsync(cuda_mask, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(cuda_mask_sum, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(cuda_dphi, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(cuda_deta, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(cuda_phi_slope, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(cuda_z0, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(cuda_sorted_M2_hits, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(cuda_reduced_M2_hits, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(cuda_reduced_M1_hits, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(cuda_edge_sum, stream));
  ACTS_CUDA_CHECK(cudaFreeAsync(cuda_nb_edges_per_doublet, stream));

  detail::CUDA_edge_data<float> edge_data;
  edge_data.nEdges = nb_graph_edges;
  edge_data.cudaEdgePtr = cuda_graph_edge_ptr;

  return edge_data;
}

}  // namespace Acts
