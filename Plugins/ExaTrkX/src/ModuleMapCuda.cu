// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/ModuleMapCuda.hpp"
#include "Acts/Plugins/ExaTrkX/detail/CudaUtils.cuh"
#include "Acts/Plugins/ExaTrkX/detail/CudaUtils.hpp"
#include "Acts/Plugins/ExaTrkX/detail/ModuleMapUtils.cuh"
#include "Acts/Plugins/ExaTrkX/detail/TensorVectorConversion.hpp"

#include <CUDA_graph_creator>
#include <CUDA_module_map_doublet>
#include <CUDA_module_map_triplet>
#include <TTree_hits>
#include <chrono>

#include <c10/cuda/CUDAGuard.h>
#include <cub/block/block_merge_sort.cuh>
#include <thrust/execution_policy.h>
#include <thrust/functional.h>
#include <thrust/scan.h>
#include <thrust/sort.h>
#include <thrust/transform_scan.h>

using Clock = std::chrono::high_resolution_clock;

using namespace torch::indexing;

namespace {

template <typename T>
class ScopedCudaPtr {
  cudaStream_t *m_stream = nullptr;
  T *m_ptr = nullptr;

 public:
  ScopedCudaPtr(std::size_t n, cudaStream_t &stream) : m_stream(&stream) {
    ACTS_CUDA_CHECK(cudaMallocAsync(&m_ptr, n * sizeof(T), stream));
  }

  ScopedCudaPtr(const ScopedCudaPtr &) = delete;
  ScopedCudaPtr(ScopedCudaPtr &&) = delete;

  ~ScopedCudaPtr() { ACTS_CUDA_CHECK(cudaFreeAsync(m_ptr, *m_stream)); }

  T *data() { return m_ptr; }
  const T *data() const { return m_ptr; }
};

}  // namespace

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
  ACTS_CUDA_CHECK(cudaMalloc(&m_cudaModuleMapKeys,
                             m_cudaModuleMapSize * sizeof(std::uint64_t)));
  ACTS_CUDA_CHECK(
      cudaMalloc(&m_cudaModuleMapVals, m_cudaModuleMapSize * sizeof(int)));

  ACTS_CUDA_CHECK(cudaMemcpy(m_cudaModuleMapKeys, keys.data(),
                             m_cudaModuleMapSize * sizeof(std::uint64_t),
                             cudaMemcpyHostToDevice));
  ACTS_CUDA_CHECK(cudaMemcpy(m_cudaModuleMapVals, vals.data(),
                             m_cudaModuleMapSize * sizeof(int),
                             cudaMemcpyHostToDevice));
}

ModuleMapCuda::~ModuleMapCuda() {
  ACTS_CUDA_CHECK(cudaFree(m_cudaModuleMapKeys));
  ACTS_CUDA_CHECK(cudaFree(m_cudaModuleMapVals));
}

namespace {}  // namespace

std::tuple<std::any, std::any, std::any> ModuleMapCuda::operator()(
    std::vector<float> &inputValues, std::size_t numNodes,
    const std::vector<std::uint64_t> &moduleIds,
    const ExecutionContext &execContext) {
  auto t0 = std::chrono::high_resolution_clock::now();

  if (moduleIds.empty()) {
    throw NoEdgesError{};
  }

  const auto nHits = moduleIds.size();
  assert(inputValues.size() % moduleIds.size() == 0);
  const auto nFeatures = inputValues.size() / moduleIds.size();
  auto &features = inputValues;

  const dim3 blockDim = m_cfg.gpuBlocks;
  const dim3 gridDimHits = (nHits + blockDim.x - 1) / blockDim.x;
  ACTS_VERBOSE("gridDimHits: " << gridDimHits.x
                               << ", blockDim: " << blockDim.x);

  // Get stream if available, otherwise use default stream
  cudaStream_t stream = cudaStreamLegacy;
  if (execContext.stream) {
    ACTS_DEBUG("Got stream " << *execContext.stream);
    stream = execContext.stream.value();
  }

  /////////////////////////
  // Prepare input data
  ////////////////////////

  // Full node features to device
  auto nodeFeatures = torch::empty(
      features.size(),
      torch::TensorOptions().device(torch::kCUDA).dtype(torch::kFloat32));
  nodeFeatures = nodeFeatures.reshape(
      {static_cast<long>(nHits), static_cast<long>(nFeatures)});
  float *cudaNodeFeaturePtr = nodeFeatures.data_ptr<float>();
  ACTS_CUDA_CHECK(cudaMemcpyAsync(cudaNodeFeaturePtr, features.data(),
                                  features.size() * sizeof(float),
                                  cudaMemcpyHostToDevice, stream));
  ACTS_VERBOSE("Slice of node features[0:9,0:9]:\n"
               << nodeFeatures.index({Slice(0, std::min(9ul, nHits)),
                                      Slice(0, std::min(9ul, nFeatures))}));

  // Module IDs to device
  ScopedCudaPtr<std::uint64_t> cudaModuleIds(nHits, stream);
  ACTS_CUDA_CHECK(cudaMemcpyAsync(cudaModuleIds.data(), moduleIds.data(),
                                  nHits * sizeof(std::uint64_t),
                                  cudaMemcpyHostToDevice, stream));

  // Allocate memory for transposed node features that are needed for the
  // module map kernels in one block
  ScopedCudaPtr<float> cudaNodeFeaturesTransposed(6 * nHits, stream);

  Acts::detail::CUDA_hit_data<float> inputData{};
  inputData.m_size = nHits;
  inputData.m_cuda_R = cudaNodeFeaturesTransposed.data() + 0 * nHits;
  inputData.m_cuda_phi = cudaNodeFeaturesTransposed.data() + 1 * nHits;
  inputData.m_cuda_z = cudaNodeFeaturesTransposed.data() + 2 * nHits;
  inputData.m_cuda_x = cudaNodeFeaturesTransposed.data() + 3 * nHits;
  inputData.m_cuda_y = cudaNodeFeaturesTransposed.data() + 4 * nHits;
  inputData.m_cuda_eta = cudaNodeFeaturesTransposed.data() + 5 * nHits;

  // Node features for module map graph
  constexpr std::size_t rOffset = 0;
  constexpr std::size_t phiOffset = 1;
  constexpr std::size_t zOffset = 2;
  constexpr std::size_t etaOffset = 3;

  const auto srcStride = sizeof(float) * nFeatures;
  const auto dstStride = sizeof(float);  // contiguous in destination
  const auto width = sizeof(float);      // only copy 1 column
  const auto height = nHits;

  ACTS_CUDA_CHECK(cudaMemcpy2DAsync(
      inputData.cuda_R(), dstStride, cudaNodeFeaturePtr + rOffset, srcStride,
      width, height, cudaMemcpyDeviceToDevice, stream));
  ACTS_CUDA_CHECK(cudaMemcpy2DAsync(
      inputData.cuda_phi(), dstStride, cudaNodeFeaturePtr + phiOffset,
      srcStride, width, height, cudaMemcpyDeviceToDevice, stream));
  ACTS_CUDA_CHECK(cudaMemcpy2DAsync(
      inputData.cuda_z(), dstStride, cudaNodeFeaturePtr + zOffset, srcStride,
      width, height, cudaMemcpyDeviceToDevice, stream));
  ACTS_CUDA_CHECK(cudaMemcpy2DAsync(
      inputData.cuda_eta(), dstStride, cudaNodeFeaturePtr + etaOffset,
      srcStride, width, height, cudaMemcpyDeviceToDevice, stream));

  // Allocate helper nb hits memory
  ScopedCudaPtr<int> cudaNbHits(m_cudaModuleMapSize + 1, stream);
  ACTS_CUDA_CHECK(cudaMemsetAsync(
      cudaNbHits.data(), 0, (m_cudaModuleMapSize + 1) * sizeof(int), stream));

  // Preprocess features
  detail::rescaleFeature<<<gridDimHits, blockDim, 0, stream>>>(
      nHits, inputData.cuda_z(), m_cfg.zScale);
  ACTS_CUDA_CHECK(cudaGetLastError());
  detail::rescaleFeature<<<gridDimHits, blockDim, 0, stream>>>(
      nHits, inputData.cuda_R(), m_cfg.rScale);
  ACTS_CUDA_CHECK(cudaGetLastError());
  detail::rescaleFeature<<<gridDimHits, blockDim, 0, stream>>>(
      nHits, inputData.cuda_phi(), m_cfg.phiScale);
  ACTS_CUDA_CHECK(cudaGetLastError());

  detail::computeXandY<<<gridDimHits, blockDim, 0, stream>>>(
      nHits, inputData.cuda_x(), inputData.cuda_y(), inputData.cuda_R(),
      inputData.cuda_phi());
  ACTS_CUDA_CHECK(cudaGetLastError());

  ScopedCudaPtr<std::uint64_t> cudaHitId(nHits, stream);
  detail::iota<<<gridDimHits, blockDim, 0, stream>>>(nHits, cudaHitId.data());
  ACTS_CUDA_CHECK(cudaGetLastError());

  detail::mapModuleIdsToNbHits<<<gridDimHits, blockDim, 0, stream>>>(
      cudaNbHits.data(), nHits, cudaModuleIds.data(), m_cudaModuleMapSize,
      m_cudaModuleMapKeys, m_cudaModuleMapVals);
  ACTS_CUDA_CHECK(cudaGetLastError());

  thrust::exclusive_scan(thrust::device.on(stream), cudaNbHits.data(),
                         cudaNbHits.data() + m_cudaModuleMapSize + 1,
                         cudaNbHits.data());
  ACTS_CUDA_CHECK(cudaGetLastError());
  int *cudaHitIndice = cudaNbHits.data();

  ///////////////////////////////////
  // Perform module map inference
  ////////////////////////////////////

  ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));
  auto t1 = std::chrono::high_resolution_clock::now();

  // TODO refactor this to avoid that inputData type in this form

  const auto edgeData = makeEdges(inputData, cudaHitIndice, stream);
  ACTS_CUDA_CHECK(cudaGetLastError());
  ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));

  auto t2 = std::chrono::high_resolution_clock::now();

  if (edgeData.nEdges == 0) {
    throw NoEdgesError{};
  }

  dim3 gridDimEdges = (edgeData.nEdges + blockDim.x - 1) / blockDim.x;
  ACTS_DEBUG("gridDimEdges: " << gridDimEdges.x
                              << ", blockDim: " << blockDim.x);

  // Make edge features
  auto edgeFeatures = torch::empty(
      6 * edgeData.nEdges,
      torch::TensorOptions().device(torch::kCUDA).dtype(torch::kFloat32));
  edgeFeatures = edgeFeatures.reshape({static_cast<long>(edgeData.nEdges), 6});
  ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));

  detail::makeEdgeFeatures<<<gridDimEdges, blockDim, 0, stream>>>(
      edgeData.nEdges, edgeData.cudaEdgePtr,
      edgeData.cudaEdgePtr + edgeData.nEdges, nFeatures, cudaNodeFeaturePtr,
      edgeFeatures.data_ptr<float>());
  ACTS_CUDA_CHECK(cudaGetLastError());
  ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));

  //  Make torch tensors
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

  ACTS_CUDA_CHECK(cudaGetLastError());
  ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));
  auto t3 = std::chrono::high_resolution_clock::now();

  auto ms = [](auto a, auto b) {
    return std::chrono::duration<double, std::milli>(b - a).count();
  };
  ACTS_DEBUG("Preparation: " << ms(t0, t1));
  ACTS_DEBUG("Inference: " << ms(t1, t2));
  ACTS_DEBUG("Postprocessing: " << ms(t2, t3));

  return {nodeFeatures, edgeIndex, edgeFeatures};
}

struct ArgsortFun {
  int *cudaPtr = nullptr;
  bool __device__ operator()(int l, int r) { return cudaPtr[l] < cudaPtr[r]; }
};

struct CastBoolToInt {
  int __device__ operator()(bool b) { return static_cast<int>(b); }
};

std::string debugPrintEdges(std::size_t nbEdges, const int *cudaSrc,
                            const int *cudaDst) {
  std::stringstream ss;
  if (nbEdges == 0) {
    return "zero edges remained";
  }
  nbEdges = std::min(10ul, nbEdges);
  std::vector<int> src(nbEdges), dst(nbEdges);
  ACTS_CUDA_CHECK(cudaDeviceSynchronize());
  ACTS_CUDA_CHECK(cudaMemcpy(src.data(), cudaSrc, nbEdges * sizeof(int),
                             cudaMemcpyDeviceToHost));
  ACTS_CUDA_CHECK(cudaMemcpy(dst.data(), cudaDst, nbEdges * sizeof(int),
                             cudaMemcpyDeviceToHost));
  for (std::size_t i = 0; i < nbEdges; ++i) {
    ss << src.at(i) << " ";
  }
  ss << "\n";
  for (std::size_t i = 0; i < nbEdges; ++i) {
    ss << dst.at(i) << " ";
  }
  return ss.str();
}

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
  // Should be definitively filled after the if block
  std::optional<ScopedCudaPtr<int>> cuda_reduced_M1_hits, cuda_reduced_M2_hits;

  ScopedCudaPtr<int> cuda_edge_sum(nb_doublets + 1, stream);

  int nb_doublet_edges{};
  // ---------------------------------------------
  // A: New method to exploit more parallelism
  // ---------------------------------------------
  if (m_cfg.moreParallel) {
    // Algorithm to build edges parallel for each hit+doublet combination
    // ==================================================================
    //
    // The motivation for this is the work imbalance for different doublets
    // By essentially pulling out the outer loop over the hits on a module
    // we have a better change that the work is more evenly distributed
    //
    // a) Assume hit ids
    // 0 1 2 3 4 5
    //
    // b) Assume module ids for hits
    // 0 0 1 1 2 3
    //
    // c) Assume doublet module map
    // 0: 0 -> 1
    // 1: 0 -> 2
    // 2: 1 -> 2
    // 3: 2 -> 3
    //
    // d) Count start hits per doublet
    // 2 2 2 1
    //
    // e) Prefix sum
    // 0 2 4 6 7

    ScopedCudaPtr<int> cuda_nb_src_hits_per_doublet(nb_doublets + 1, stream);

    detail::count_src_hits_per_doublet<<<grid_dim, block_dim, 0, stream>>>(
        nb_doublets, m_cudaModuleMapDoublet->cuda_module1(), cuda_hit_indice,
        cuda_nb_src_hits_per_doublet.data());
    ACTS_CUDA_CHECK(cudaGetLastError());

    thrust::exclusive_scan(
        thrust::device.on(stream), cuda_nb_src_hits_per_doublet.data(),
        cuda_nb_src_hits_per_doublet.data() + nb_doublets + 1,
        cuda_nb_src_hits_per_doublet.data());

    int sum_nb_src_hits_per_doublet{};
    ACTS_CUDA_CHECK(
        cudaMemcpyAsync(&sum_nb_src_hits_per_doublet,
                        &cuda_nb_src_hits_per_doublet.data()[nb_doublets],
                        sizeof(int), cudaMemcpyDeviceToHost, stream));
    ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));
    ACTS_DEBUG("sum_nb_hits_per_doublet: " << sum_nb_src_hits_per_doublet);

    if (sum_nb_src_hits_per_doublet == 0) {
      throw NoEdgesError{};
    }

    ScopedCudaPtr<int> cuda_edge_sum_per_src_hit(
        sum_nb_src_hits_per_doublet + 1, stream);

    dim3 grid_dim_shpd =
        ((sum_nb_src_hits_per_doublet + block_dim.x - 1) / block_dim.x);
    detail::count_doublet_edges_new<float>
        <<<grid_dim_shpd, block_dim, 0, stream>>>(
            sum_nb_src_hits_per_doublet, nb_doublets,
            cuda_nb_src_hits_per_doublet.data(),
            m_cudaModuleMapDoublet->cuda_module1(),
            m_cudaModuleMapDoublet->cuda_module2(), cuda_TThits.cuda_R(),
            cuda_TThits.cuda_z(), cuda_TThits.cuda_eta(),
            cuda_TThits.cuda_phi(), m_cudaModuleMapDoublet->cuda_z0_min(),
            m_cudaModuleMapDoublet->cuda_z0_max(),
            m_cudaModuleMapDoublet->cuda_deta_min(),
            m_cudaModuleMapDoublet->cuda_deta_max(),
            m_cudaModuleMapDoublet->cuda_phi_slope_min(),
            m_cudaModuleMapDoublet->cuda_phi_slope_max(),
            m_cudaModuleMapDoublet->cuda_dphi_min(),
            m_cudaModuleMapDoublet->cuda_dphi_max(), cuda_hit_indice,
            cuda_edge_sum_per_src_hit.data(), m_cfg.epsilon);
    ACTS_CUDA_CHECK(cudaGetLastError());

    thrust::exclusive_scan(
        thrust::device.on(stream), cuda_edge_sum_per_src_hit.data(),
        cuda_edge_sum_per_src_hit.data() + sum_nb_src_hits_per_doublet + 1,
        cuda_edge_sum_per_src_hit.data());

    ACTS_CUDA_CHECK(cudaMemcpyAsync(
        &nb_doublet_edges,
        &cuda_edge_sum_per_src_hit.data()[sum_nb_src_hits_per_doublet],
        sizeof(int), cudaMemcpyDeviceToHost, stream));
    ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));
    ACTS_DEBUG("nb_doublet_edges: " << nb_doublet_edges);
    ACTS_DEBUG("Allocate " << (2ul * nb_doublet_edges * sizeof(int)) * 1.0e-6
                           << " MB for edges");
    cuda_reduced_M1_hits.emplace(nb_doublet_edges, stream);
    cuda_reduced_M2_hits.emplace(nb_doublet_edges, stream);

    detail::build_doublet_edges_new<float>
        <<<grid_dim_shpd, block_dim, 0, stream>>>(
            sum_nb_src_hits_per_doublet, nb_doublets,
            cuda_nb_src_hits_per_doublet.data(),
            m_cudaModuleMapDoublet->cuda_module1(),
            m_cudaModuleMapDoublet->cuda_module2(), cuda_TThits.cuda_R(),
            cuda_TThits.cuda_z(), cuda_TThits.cuda_eta(),
            cuda_TThits.cuda_phi(), m_cudaModuleMapDoublet->cuda_z0_min(),
            m_cudaModuleMapDoublet->cuda_z0_max(),
            m_cudaModuleMapDoublet->cuda_deta_min(),
            m_cudaModuleMapDoublet->cuda_deta_max(),
            m_cudaModuleMapDoublet->cuda_phi_slope_min(),
            m_cudaModuleMapDoublet->cuda_phi_slope_max(),
            m_cudaModuleMapDoublet->cuda_dphi_min(),
            m_cudaModuleMapDoublet->cuda_dphi_max(), cuda_hit_indice,
            cuda_reduced_M1_hits->data(), cuda_reduced_M2_hits->data(),
            cuda_edge_sum_per_src_hit.data(), m_cfg.epsilon);
    ACTS_CUDA_CHECK(cudaGetLastError());

    detail::computeDoubletEdgeSum<<<grid_dim, block_dim, 0, stream>>>(
        nb_doublets, cuda_nb_src_hits_per_doublet.data(),
        cuda_edge_sum_per_src_hit.data(), cuda_edge_sum.data());
    ACTS_CUDA_CHECK(cudaGetLastError());
  } else {
    // Allocate one integer on device and set it to 0
    ScopedCudaPtr<int> cuda_nb_doublet_edges(1, stream);
    ACTS_CUDA_CHECK(
        cudaMemsetAsync(cuda_nb_doublet_edges.data(), 0, sizeof(int), stream));

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
        cuda_nb_doublet_edges.data(), cuda_edge_sum.data(), m_cfg.epsilon);
    ACTS_CUDA_CHECK(cudaGetLastError());

    // Copy the number of edges to the host, synchronize and allocate
    ACTS_CUDA_CHECK(cudaMemcpyAsync(&nb_doublet_edges,
                                    cuda_nb_doublet_edges.data(), sizeof(int),
                                    cudaMemcpyDeviceToHost, stream));
    ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));
    ACTS_DEBUG("nb_doublet_edges: " << nb_doublet_edges);

    ACTS_DEBUG("Allocate " << (2ul * nb_doublet_edges * sizeof(int)) * 1.0e-6
                           << " MB for edges");
    cuda_reduced_M1_hits.emplace(nb_doublet_edges, stream);
    cuda_reduced_M2_hits.emplace(nb_doublet_edges, stream);

    // Prefix sum to get the edge offset for each doublet
    thrust::exclusive_scan(thrust::device.on(stream), cuda_edge_sum.data(),
                           cuda_edge_sum.data() + (nb_doublets + 1),
                           cuda_edge_sum.data());

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
        cuda_reduced_M1_hits->data(), cuda_reduced_M2_hits->data(),
        cuda_edge_sum.data(), m_cfg.epsilon);
    ACTS_CUDA_CHECK(cudaGetLastError());
  }

  ACTS_VERBOSE("First 10 doublet edges:\n"
               << debugPrintEdges(nb_doublet_edges,
                                  cuda_reduced_M1_hits->data(),
                                  cuda_reduced_M2_hits->data()));
  if (nb_doublet_edges == 0) {
    throw NoEdgesError{};
  }
  //---------------------------------------------
  // reduce nb of hits and sort by hid id M2 hits
  //---------------------------------------------

  ScopedCudaPtr<int> cuda_sorted_M2_hits(nb_doublet_edges, stream);

  grid_dim = ((nb_doublet_edges + block_dim.x - 1) / block_dim.x);
  init_vector<<<grid_dim, block_dim, 0, stream>>>(cuda_sorted_M2_hits.data(),
                                                  nb_doublet_edges);
  ACTS_CUDA_CHECK(cudaGetLastError());

  if (m_cfg.moreParallel) {
    dim3 grid_dim_sort = nb_doublets;
    dim3 block_dim_even_odd = 64;
    ACTS_DEBUG("Using block_odd_even_sort, grid_dim.x = "
               << nb_doublets << ", block_dim.x = " << block_dim_even_odd.x);
    detail::block_odd_even_sort<<<nb_doublets, block_dim_even_odd, 0, stream>>>(
        cuda_reduced_M2_hits->data(), cuda_edge_sum.data(),
        cuda_sorted_M2_hits.data());
  } else {
    grid_dim = ((nb_doublets + block_dim.x - 1) / block_dim.x);
    partial_quick_sort<<<grid_dim, block_dim, 0, stream>>>(
        cuda_sorted_M2_hits.data(), cuda_reduced_M2_hits->data(),
        cuda_edge_sum.data(), nb_doublets);
    ACTS_CUDA_CHECK(cudaGetLastError());
  }

  // -----------------------------
  // build doublets geometric cuts
  // -----------------------------

  ScopedCudaPtr<float> cuda_z0(nb_doublet_edges, stream),
      cuda_phi_slope(nb_doublet_edges, stream),
      cuda_deta(nb_doublet_edges, stream), cuda_dphi(nb_doublet_edges, stream);
  grid_dim = ((nb_doublet_edges + block_dim.x - 1) / block_dim.x);
  hits_geometric_cuts<<<grid_dim, block_dim, 0, stream>>>(
      cuda_z0.data(), cuda_phi_slope.data(), cuda_deta.data(), cuda_dphi.data(),
      cuda_reduced_M1_hits->data(), cuda_reduced_M2_hits->data(),
      cuda_TThits.cuda_R(), cuda_TThits.cuda_z(), cuda_TThits.cuda_eta(),
      cuda_TThits.cuda_phi(), TMath::Pi(), nb_doublet_edges, m_cfg.epsilon);
  ACTS_CUDA_CHECK(cudaGetLastError());
  ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));

  ScopedCudaPtr<bool> cuda_mask(nb_doublet_edges + 1, stream);
  ACTS_CUDA_CHECK(
      cudaMemset(cuda_mask.data(), 0, (nb_doublet_edges + 1) * sizeof(bool)));
  ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));

  // -------------------------
  // loop over module triplets
  // -------------------------
  int nb_triplets = m_cudaModuleMapTriplet->size();
  grid_dim = ((nb_triplets + block_dim.x - 1) / block_dim.x);

  ScopedCudaPtr<bool> cuda_vertices(cuda_TThits.size(), stream);
  ACTS_CUDA_CHECK(cudaMemsetAsync(cuda_vertices.data(), 0,
                                  cuda_TThits.size() * sizeof(bool), stream));
  ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));

  if (m_cfg.moreParallel) {
    // Allocate memory for the number of hits per triplet
    ScopedCudaPtr<int> cuda_src_hits_per_triplet(nb_triplets + 1, stream);

    detail::count_triplet_hits<<<grid_dim, block_dim, 0, stream>>>(
        nb_triplets, m_cudaModuleMapTriplet->cuda_module12_map(),
        m_cudaModuleMapTriplet->cuda_module23_map(), cuda_edge_sum.data(),
        cuda_src_hits_per_triplet.data());
    ACTS_CUDA_CHECK(cudaGetLastError());

    // Perform prefix sum to get the offset for each triplet
    thrust::exclusive_scan(thrust::device.on(stream),
                           cuda_src_hits_per_triplet.data(),
                           cuda_src_hits_per_triplet.data() + nb_triplets + 1,
                           cuda_src_hits_per_triplet.data());

    int nb_src_hits_per_triplet_sum{};
    ACTS_CUDA_CHECK(
        cudaMemcpyAsync(&nb_src_hits_per_triplet_sum,
                        &cuda_src_hits_per_triplet.data()[nb_triplets],
                        sizeof(int), cudaMemcpyDeviceToHost, stream));
    ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));
    ACTS_DEBUG("nb_src_hits_per_triplet_sum: " << nb_src_hits_per_triplet_sum);
    if (nb_src_hits_per_triplet_sum == 0) {
      throw NoEdgesError{};
    }

    dim3 grid_dim_shpt =
        ((nb_src_hits_per_triplet_sum + block_dim.x - 1) / block_dim.x);

    detail::triplet_cuts_new2<float><<<grid_dim_shpt, block_dim, 0, stream>>>(
        nb_src_hits_per_triplet_sum, nb_triplets,
        cuda_src_hits_per_triplet.data(),
        m_cudaModuleMapTriplet->cuda_module12_map(),
        m_cudaModuleMapTriplet->cuda_module23_map(), cuda_TThits.cuda_x(),
        cuda_TThits.cuda_y(), cuda_TThits.cuda_z(), cuda_TThits.cuda_R(),
        cuda_z0.data(), cuda_phi_slope.data(), cuda_deta.data(),
        cuda_dphi.data(), m_cudaModuleMapTriplet->module12().cuda_z0_min(),
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
        cuda_reduced_M1_hits->data(), cuda_reduced_M2_hits->data(),
        cuda_sorted_M2_hits.data(), cuda_edge_sum.data(), cuda_vertices.data(),
        cuda_mask.data(), m_cfg.epsilon);
    ACTS_CUDA_CHECK(cudaGetLastError());
  } else {
    Acts::detail::triplet_cuts_new<float><<<grid_dim, block_dim, 0, stream>>>(
        nb_triplets, m_cudaModuleMapTriplet->cuda_module12_map(),
        m_cudaModuleMapTriplet->cuda_module23_map(), cuda_TThits.cuda_x(),
        cuda_TThits.cuda_y(), cuda_TThits.cuda_z(), cuda_TThits.cuda_R(),
        cuda_z0.data(), cuda_phi_slope.data(), cuda_deta.data(),
        cuda_dphi.data(), m_cudaModuleMapTriplet->module12().cuda_z0_min(),
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
        cuda_reduced_M1_hits->data(), cuda_reduced_M2_hits->data(),
        cuda_sorted_M2_hits.data(), cuda_edge_sum.data(), cuda_vertices.data(),
        cuda_mask.data(), m_cfg.epsilon);
    ACTS_CUDA_CHECK(cudaGetLastError());
  }

  //----------------
  // edges reduction
  //----------------
  ScopedCudaPtr<int> cuda_mask_sum(nb_doublet_edges + 1, stream);
  ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));
  thrust::transform_exclusive_scan(thrust::device.on(stream), cuda_mask.data(),
                                   cuda_mask.data() + (nb_doublet_edges + 1),
                                   cuda_mask_sum.data(), CastBoolToInt{}, 0,
                                   thrust::plus<int>());

  int nb_graph_edges{};
  ACTS_CUDA_CHECK(cudaMemcpyAsync(&nb_graph_edges,
                                  &cuda_mask_sum.data()[nb_doublet_edges],
                                  sizeof(int), cudaMemcpyDeviceToHost, stream));
  ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));

  ACTS_DEBUG("nb_graph_edges: " << nb_graph_edges);
  if (nb_graph_edges == 0) {
    throw NoEdgesError{};
  }

  // Leave this as a bare pointer for now, since there is only very simple
  // control flow after here and we can keep the interface clean form the
  // ScopedCudaPtr type
  int *cuda_graph_edge_ptr{};
  ACTS_CUDA_CHECK(cudaMallocAsync(&cuda_graph_edge_ptr,
                                  2 * nb_graph_edges * sizeof(int), stream));
  int *cuda_graph_M1_hits = cuda_graph_edge_ptr;
  int *cuda_graph_M2_hits = cuda_graph_edge_ptr + nb_graph_edges;

  grid_dim = ((nb_doublet_edges + block_dim.x - 1) / block_dim.x);
  ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));
  compact_stream<<<grid_dim, block_dim, 0, stream>>>(
      cuda_graph_M1_hits, cuda_reduced_M1_hits->data(), cuda_mask.data(),
      cuda_mask_sum.data(), nb_doublet_edges);
  ACTS_CUDA_CHECK(cudaGetLastError());
  compact_stream<<<grid_dim, block_dim, 0, stream>>>(
      cuda_graph_M2_hits, cuda_reduced_M2_hits->data(), cuda_mask.data(),
      cuda_mask_sum.data(), nb_doublet_edges);
  ACTS_CUDA_CHECK(cudaGetLastError());

  ACTS_VERBOSE("First 10 doublet edges:\n"
               << debugPrintEdges(nb_graph_edges, cuda_graph_M1_hits,
                                  cuda_graph_M2_hits));

  detail::CUDA_edge_data<float> edge_data{};
  edge_data.nEdges = nb_graph_edges;
  edge_data.cudaEdgePtr = cuda_graph_edge_ptr;

  return edge_data;
}

}  // namespace Acts
