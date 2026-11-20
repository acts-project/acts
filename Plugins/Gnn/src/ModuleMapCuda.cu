// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Gnn/ModuleMapCuda.hpp"
#include "ActsPlugins/Gnn/detail/CudaUtils.cuh"
#include "ActsPlugins/Gnn/detail/CudaUtils.hpp"
#include "ActsPlugins/Gnn/detail/ModuleMapKernels.cuh"
#include "ActsPlugins/Gnn/detail/ModuleMapUtils.cuh"

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <MMG/CUDA_module_map_doublet>
#include <MMG/CUDA_module_map_triplet>
#include <cub/device/device_scan.cuh>
#include <thrust/iterator/transform_iterator.h>

// The double-precision reference arithmetic in ModuleMapKernels.cuh matches
// the module map's data_type2
static_assert(std::is_same_v<data_type2, double>);

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

template <typename T>
struct CUDA_hit_data {
  T *m_cuda_x;
  T *m_cuda_y;
  T *m_cuda_R;
  T *m_cuda_phi;
  T *m_cuda_z;
  T *m_cuda_eta;

  T *cuda_x() { return m_cuda_x; }
  T *cuda_y() { return m_cuda_y; }
  T *cuda_z() { return m_cuda_z; }
  T *cuda_R() { return m_cuda_R; }
  T *cuda_phi() { return m_cuda_phi; }
  T *cuda_eta() { return m_cuda_eta; }
};

struct CastBoolToInt {
  int __host__ __device__ operator()(bool b) const {
    return static_cast<int>(b);
  }
};

/// Exclusive prefix sum with CUB. Unlike thrust with the default policy, this
/// takes its scratch memory from the stream-ordered allocator and does not
/// synchronize the stream. In-place operation (in == out) is supported.
template <typename InputIt, typename OutputIt>
void exclusiveSum(InputIt in, OutputIt out, std::size_t n,
                  cudaStream_t &stream) {
  std::size_t tmpBytes = 0;
  ACTS_CUDA_CHECK(
      cub::DeviceScan::ExclusiveSum(nullptr, tmpBytes, in, out, n, stream));
  ScopedCudaPtr<std::byte> tmp(tmpBytes, stream);
  ACTS_CUDA_CHECK(
      cub::DeviceScan::ExclusiveSum(tmp.data(), tmpBytes, in, out, n, stream));
}

template <typename T>
std::string debugPrintEdges(std::size_t nbEdges, const T *cudaSrc,
                            const T *cudaDst) {
  std::stringstream ss;
  if (nbEdges == 0) {
    return "zero edges remained";
  }
  nbEdges = std::min(10ul, nbEdges);
  std::vector<T> src(nbEdges), dst(nbEdges);
  ACTS_CUDA_CHECK(cudaDeviceSynchronize());
  ACTS_CUDA_CHECK(cudaMemcpy(src.data(), cudaSrc, nbEdges * sizeof(T),
                             cudaMemcpyDeviceToHost));
  ACTS_CUDA_CHECK(cudaMemcpy(dst.data(), cudaDst, nbEdges * sizeof(T),
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

}  // namespace

using namespace Acts;

namespace ActsPlugins {

class ModuleMapCuda::Impl {
 public:
  std::unique_ptr<CUDA_module_map_doublet<float>> cudaModuleMapDoublet;
  std::unique_ptr<CUDA_module_map_triplet<float>> cudaModuleMapTriplet;

  std::uint64_t *cudaModuleMapKeys{};
  int *cudaModuleMapVals{};
  std::size_t cudaModuleMapSize{};

  /// Returns the [2, nEdges] edge index and the [nEdges, 6] edge features
  std::pair<Tensor<std::int64_t>, Tensor<float>> makeEdges(
      CUDA_hit_data<float> cuda_TThits, int *cuda_hit_indice,
      const float *cudaNodeFeatures, std::size_t nNodeFeatures,
      cudaStream_t &stream, const ExecutionContext &execContext,
      const ModuleMapCuda::Config &cfg, const Logger &logger) const;
};

ModuleMapCuda::ModuleMapCuda(const Config &cfg,
                             std::unique_ptr<const Logger> logger_)
    : m_impl(std::make_unique<Impl>()),
      m_cfg(cfg),
      m_logger(std::move(logger_)) {
  if (m_cfg.gpuBlocks <= 0 || m_cfg.gpuBlocks > detail::kMaxBlockSize) {
    throw std::invalid_argument(
        "ModuleMapCuda: gpuBlocks (threads per block) must be in [1, " +
        std::to_string(detail::kMaxBlockSize) + "], got " +
        std::to_string(m_cfg.gpuBlocks));
  }

  module_map_triplet<float> moduleMapCpu;
  moduleMapCpu.read_TTree(cfg.moduleMapPath.c_str());
  if (!moduleMapCpu) {
    throw std::runtime_error("Cannot retrieve ModuleMap from " +
                             cfg.moduleMapPath);
  }

  ACTS_DEBUG("ModuleMap GPU block dim: " << m_cfg.gpuBlocks);

  if (m_cfg.retainMemPool) {
    // By default the stream-ordered memory pool returns all free memory to the
    // system at every synchronization, so each cudaMallocAsync after a sync
    // becomes a real allocation. Keep the memory in the pool instead.
    cudaMemPool_t memPool{};
    ACTS_CUDA_CHECK(cudaDeviceGetDefaultMemPool(&memPool, m_cfg.gpuDevice));
    std::uint64_t threshold = std::numeric_limits<std::uint64_t>::max();
    ACTS_CUDA_CHECK(cudaMemPoolSetAttribute(
        memPool, cudaMemPoolAttrReleaseThreshold, &threshold));
  }

  m_impl->cudaModuleMapDoublet =
      std::make_unique<CUDA_module_map_doublet<float>>(moduleMapCpu);
  m_impl->cudaModuleMapDoublet->HostToDevice();
  m_impl->cudaModuleMapTriplet =
      std::make_unique<CUDA_module_map_triplet<float>>(moduleMapCpu);
  m_impl->cudaModuleMapTriplet->HostToDevice();

  ACTS_DEBUG("# of modules = " << moduleMapCpu.module_map().size());

  // check if we actually have a module map
  std::map<std::uint64_t, int> test;

  std::vector<std::uint64_t> keys;
  keys.reserve(m_impl->cudaModuleMapDoublet->module_map().size());
  std::vector<int> vals;
  vals.reserve(m_impl->cudaModuleMapDoublet->module_map().size());

  for (auto [key, value] : m_impl->cudaModuleMapDoublet->module_map()) {
    auto [it, success] = test.insert({key, value});
    if (!success) {
      throw std::runtime_error("Duplicate key in module map");
    }
    keys.push_back(key);
    vals.push_back(value);
  }

  // copy module map to device
  m_impl->cudaModuleMapSize = m_impl->cudaModuleMapDoublet->module_map().size();
  ACTS_CUDA_CHECK(
      cudaMalloc(&m_impl->cudaModuleMapKeys,
                 m_impl->cudaModuleMapSize * sizeof(std::uint64_t)));
  ACTS_CUDA_CHECK(cudaMalloc(&m_impl->cudaModuleMapVals,
                             m_impl->cudaModuleMapSize * sizeof(int)));

  ACTS_CUDA_CHECK(cudaMemcpy(m_impl->cudaModuleMapKeys, keys.data(),
                             m_impl->cudaModuleMapSize * sizeof(std::uint64_t),
                             cudaMemcpyHostToDevice));
  ACTS_CUDA_CHECK(cudaMemcpy(m_impl->cudaModuleMapVals, vals.data(),
                             m_impl->cudaModuleMapSize * sizeof(int),
                             cudaMemcpyHostToDevice));
}

ModuleMapCuda::~ModuleMapCuda() {
  ACTS_CUDA_CHECK(cudaFree(m_impl->cudaModuleMapKeys));
  ACTS_CUDA_CHECK(cudaFree(m_impl->cudaModuleMapVals));
}

PipelineTensors ModuleMapCuda::operator()(
    std::vector<float> &inputValues, std::size_t numNodes,
    const std::vector<std::uint64_t> &moduleIds,
    const ExecutionContext &execContext) {
  auto t0 = std::chrono::high_resolution_clock::now();
  assert(execContext.device.isCuda());

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

  auto nodeFeatures = Tensor<float>::Create({nHits, nFeatures}, execContext);
  float *cudaNodeFeaturePtr = nodeFeatures.data();
  ACTS_CUDA_CHECK(cudaMemcpyAsync(cudaNodeFeaturePtr, features.data(),
                                  features.size() * sizeof(float),
                                  cudaMemcpyHostToDevice, stream));

  // Module IDs to device
  ScopedCudaPtr<std::uint64_t> cudaModuleIds(nHits, stream);
  ACTS_CUDA_CHECK(cudaMemcpyAsync(cudaModuleIds.data(), moduleIds.data(),
                                  nHits * sizeof(std::uint64_t),
                                  cudaMemcpyHostToDevice, stream));

  // Allocate memory for transposed node features that are needed for the
  // module map kernels in one block
  ScopedCudaPtr<float> cudaNodeFeaturesTransposed(6 * nHits, stream);

  CUDA_hit_data<float> inputData{};
  inputData.m_cuda_R = cudaNodeFeaturesTransposed.data() + 0 * nHits;
  inputData.m_cuda_phi = cudaNodeFeaturesTransposed.data() + 1 * nHits;
  inputData.m_cuda_z = cudaNodeFeaturesTransposed.data() + 2 * nHits;
  inputData.m_cuda_x = cudaNodeFeaturesTransposed.data() + 3 * nHits;
  inputData.m_cuda_y = cudaNodeFeaturesTransposed.data() + 4 * nHits;
  inputData.m_cuda_eta = cudaNodeFeaturesTransposed.data() + 5 * nHits;

  // Allocate helper nb hits memory
  ScopedCudaPtr<int> cudaNbHits(m_impl->cudaModuleMapSize + 1, stream);
  ACTS_CUDA_CHECK(cudaMemsetAsync(cudaNbHits.data(), 0,
                                  (m_impl->cudaModuleMapSize + 1) * sizeof(int),
                                  stream));

  detail::preprocessHitFeatures<<<gridDimHits, blockDim, 0, stream>>>(
      nHits, nFeatures, cudaNodeFeaturePtr, inputData.cuda_R(),
      inputData.cuda_phi(), inputData.cuda_z(), inputData.cuda_eta(),
      inputData.cuda_x(), inputData.cuda_y(), m_cfg.rScale, m_cfg.phiScale,
      m_cfg.zScale);
  ACTS_CUDA_CHECK(cudaGetLastError());

  detail::mapModuleIdsToNbHits<<<gridDimHits, blockDim, 0, stream>>>(
      cudaNbHits.data(), nHits, cudaModuleIds.data(), m_impl->cudaModuleMapSize,
      m_impl->cudaModuleMapKeys, m_impl->cudaModuleMapVals);
  ACTS_CUDA_CHECK(cudaGetLastError());

  exclusiveSum(cudaNbHits.data(), cudaNbHits.data(),
               m_impl->cudaModuleMapSize + 1, stream);
  int *cudaHitIndice = cudaNbHits.data();

  ///////////////////////////////////
  // Perform module map inference
  ////////////////////////////////////

  if (m_cfg.debugSynchronize) {
    ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));
  }
  auto t1 = std::chrono::high_resolution_clock::now();

  // Builds the edges, and in its final pass writes the edge index and edge
  // features directly
  auto [edgeIndex, edgeFeatures] =
      m_impl->makeEdges(inputData, cudaHitIndice, cudaNodeFeaturePtr, nFeatures,
                        stream, execContext, m_cfg, logger());
  ACTS_CUDA_CHECK(cudaGetLastError());
  if (m_cfg.debugSynchronize) {
    ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));
  }

  auto t2 = std::chrono::high_resolution_clock::now();

  auto ms = [](auto a, auto b) {
    return std::chrono::duration<double, std::milli>(b - a).count();
  };
  ACTS_DEBUG("Preparation: " << ms(t0, t1));
  ACTS_DEBUG("Inference + postprocessing: " << ms(t1, t2));

  return {std::move(nodeFeatures),
          std::move(edgeIndex),
          std::move(edgeFeatures),
          {}};
}

std::pair<Tensor<std::int64_t>, Tensor<float>> ModuleMapCuda::Impl::makeEdges(
    CUDA_hit_data<float> cuda_TThits, int *cuda_hit_indice,
    const float *cudaNodeFeatures, std::size_t nNodeFeatures,
    cudaStream_t &stream, const ExecutionContext &execContext,
    const ModuleMapCuda::Config &cfg, const Logger &logger) const {
  const dim3 block_dim = cfg.gpuBlocks;
  auto gridFor = [&](std::size_t nThreads) {
    return dim3(
        static_cast<unsigned int>((nThreads + block_dim.x - 1) / block_dim.x));
  };

  // ---------------------------------------------
  // count source hits per module doublet (sync 1)
  // ---------------------------------------------
  const int nb_doublets = cudaModuleMapDoublet->size();
  ACTS_DEBUG("nb doublets " << nb_doublets);

  ScopedCudaPtr<int> cuda_nb_src_hits_per_doublet(nb_doublets + 1, stream);

  detail::
      countSourceHitsPerDoublet<<<gridFor(nb_doublets), block_dim, 0, stream>>>(
          cuda_nb_src_hits_per_doublet.data(),
          cudaModuleMapDoublet->cuda_module1(), cuda_hit_indice, nb_doublets);
  ACTS_CUDA_CHECK(cudaGetLastError());

  exclusiveSum(cuda_nb_src_hits_per_doublet.data(),
               cuda_nb_src_hits_per_doublet.data(), nb_doublets + 1, stream);

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

  // ------------------------------------------------------------
  // count doublet edges and triplet work items together (sync 2)
  // ------------------------------------------------------------
  ScopedCudaPtr<int> cuda_src_work_to_doublet(sum_nb_src_hits_per_doublet,
                                              stream);
  constexpr int kThreadsPerDoublet = 32;
  detail::buildWorkToItem<kThreadsPerDoublet>
      <<<gridFor(static_cast<std::size_t>(nb_doublets) * kThreadsPerDoublet),
         block_dim, 0, stream>>>(cuda_src_work_to_doublet.data(),
                                 cuda_nb_src_hits_per_doublet.data(),
                                 nb_doublets);
  ACTS_CUDA_CHECK(cudaGetLastError());

  ScopedCudaPtr<int> cuda_edge_sum_per_src_hit(sum_nb_src_hits_per_doublet + 1,
                                               stream);
  // Accepted hit pairs per source hit, filled by the count pass and reused by
  // the build pass
  ScopedCudaPtr<std::uint64_t> cuda_pair_masks(sum_nb_src_hits_per_doublet,
                                               stream);
  detail::countDoubletEdges<float>
      <<<gridFor(sum_nb_src_hits_per_doublet), block_dim, 0, stream>>>(
          cuda_edge_sum_per_src_hit.data(), cuda_pair_masks.data(),
          cuda_src_work_to_doublet.data(), cuda_nb_src_hits_per_doublet.data(),
          cudaModuleMapDoublet->cuda_module1(),
          cudaModuleMapDoublet->cuda_module2(), cuda_TThits.cuda_R(),
          cuda_TThits.cuda_z(), cuda_TThits.cuda_eta(), cuda_TThits.cuda_phi(),
          cudaModuleMapDoublet->cuda_z0_min(),
          cudaModuleMapDoublet->cuda_deta_min(),
          cudaModuleMapDoublet->cuda_phi_slope_min(),
          cudaModuleMapDoublet->cuda_dphi_min(),
          cudaModuleMapDoublet->cuda_z0_max(),
          cudaModuleMapDoublet->cuda_deta_max(),
          cudaModuleMapDoublet->cuda_phi_slope_max(),
          cudaModuleMapDoublet->cuda_dphi_max(), cuda_hit_indice, detail::g_pi,
          cfg.epsilon, sum_nb_src_hits_per_doublet);
  ACTS_CUDA_CHECK(cudaGetLastError());

  exclusiveSum(cuda_edge_sum_per_src_hit.data(),
               cuda_edge_sum_per_src_hit.data(),
               sum_nb_src_hits_per_doublet + 1, stream);

  ScopedCudaPtr<int> cuda_active_src_flags(sum_nb_src_hits_per_doublet + 1,
                                           stream);
  ScopedCudaPtr<int> cuda_active_src_offsets(sum_nb_src_hits_per_doublet + 1,
                                             stream);
  detail::markActiveSrcWork<<<gridFor(sum_nb_src_hits_per_doublet), block_dim,
                              0, stream>>>(cuda_active_src_flags.data(),
                                           cuda_edge_sum_per_src_hit.data(),
                                           sum_nb_src_hits_per_doublet);
  ACTS_CUDA_CHECK(cudaGetLastError());
  exclusiveSum(cuda_active_src_flags.data(), cuda_active_src_offsets.data(),
               sum_nb_src_hits_per_doublet + 1, stream);

  // Edge offsets per module doublet. These only depend on the edge counts, so
  // the triplet work can be counted before the doublet edges are built.
  ScopedCudaPtr<int> cuda_edge_sum(nb_doublets + 1, stream);
  detail::doubletEdgeSum<<<gridFor(nb_doublets + 1), block_dim, 0, stream>>>(
      cuda_edge_sum.data(), cuda_nb_src_hits_per_doublet.data(),
      cuda_edge_sum_per_src_hit.data(), nb_doublets);
  ACTS_CUDA_CHECK(cudaGetLastError());

  const int nb_triplets = cudaModuleMapTriplet->size();
  ScopedCudaPtr<int> cuda_src_hits_per_triplet(nb_triplets + 1, stream);
  detail::countTripletHits<<<gridFor(nb_triplets), block_dim, 0, stream>>>(
      cuda_src_hits_per_triplet.data(),
      cudaModuleMapTriplet->cuda_module12_map(),
      cudaModuleMapTriplet->cuda_module23_map(), cuda_edge_sum.data(),
      nb_triplets);
  ACTS_CUDA_CHECK(cudaGetLastError());

  exclusiveSum(cuda_src_hits_per_triplet.data(),
               cuda_src_hits_per_triplet.data(), nb_triplets + 1, stream);

  int nb_doublet_edges{};
  int nb_active_src{};
  int nb_src_hits_per_triplet_sum{};
  ACTS_CUDA_CHECK(cudaMemcpyAsync(
      &nb_doublet_edges,
      &cuda_edge_sum_per_src_hit.data()[sum_nb_src_hits_per_doublet],
      sizeof(int), cudaMemcpyDeviceToHost, stream));
  ACTS_CUDA_CHECK(cudaMemcpyAsync(
      &nb_active_src,
      &cuda_active_src_offsets.data()[sum_nb_src_hits_per_doublet], sizeof(int),
      cudaMemcpyDeviceToHost, stream));
  ACTS_CUDA_CHECK(
      cudaMemcpyAsync(&nb_src_hits_per_triplet_sum,
                      &cuda_src_hits_per_triplet.data()[nb_triplets],
                      sizeof(int), cudaMemcpyDeviceToHost, stream));
  ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));
  ACTS_DEBUG("nb_doublet_edges: " << nb_doublet_edges);
  ACTS_DEBUG("nb_active_src: " << nb_active_src);
  ACTS_DEBUG("nb_src_hits_per_triplet_sum: " << nb_src_hits_per_triplet_sum);

  if (nb_doublet_edges == 0 || nb_src_hits_per_triplet_sum == 0) {
    throw NoEdgesError{};
  }

  // ------------------
  // build doublet edges
  // ------------------
  ACTS_DEBUG("Allocate " << (2ul * nb_doublet_edges * sizeof(int)) * 1.0e-6
                         << " MB for edges");
  ScopedCudaPtr<int> cuda_reduced_M1_hits(nb_doublet_edges, stream);
  ScopedCudaPtr<int> cuda_reduced_M2_hits(nb_doublet_edges, stream);

  // nb_doublet_edges > 0 implies nb_active_src > 0
  ScopedCudaPtr<int> cuda_active_src_work(nb_active_src, stream);
  detail::compactActiveDoublets<<<gridFor(sum_nb_src_hits_per_doublet),
                                  block_dim, 0, stream>>>(
      cuda_active_src_work.data(), cuda_active_src_flags.data(),
      cuda_active_src_offsets.data(), sum_nb_src_hits_per_doublet);
  ACTS_CUDA_CHECK(cudaGetLastError());

  detail::buildDoubletEdgesActive<float>
      <<<gridFor(nb_active_src), block_dim, 0, stream>>>(
          cuda_reduced_M1_hits.data(), cuda_reduced_M2_hits.data(),
          nb_active_src, cuda_pair_masks.data(), cuda_active_src_work.data(),
          cuda_src_work_to_doublet.data(), cuda_nb_src_hits_per_doublet.data(),
          cuda_hit_indice, cuda_edge_sum_per_src_hit.data(),
          cudaModuleMapDoublet->cuda_module1(),
          cudaModuleMapDoublet->cuda_module2(), cuda_TThits.cuda_R(),
          cuda_TThits.cuda_z(), cuda_TThits.cuda_eta(), cuda_TThits.cuda_phi(),
          cudaModuleMapDoublet->cuda_z0_min(),
          cudaModuleMapDoublet->cuda_deta_min(),
          cudaModuleMapDoublet->cuda_phi_slope_min(),
          cudaModuleMapDoublet->cuda_dphi_min(),
          cudaModuleMapDoublet->cuda_z0_max(),
          cudaModuleMapDoublet->cuda_deta_max(),
          cudaModuleMapDoublet->cuda_phi_slope_max(),
          cudaModuleMapDoublet->cuda_dphi_max(), detail::g_pi, cfg.epsilon);
  ACTS_CUDA_CHECK(cudaGetLastError());

  ACTS_VERBOSE("First 10 doublet edges:\n"
               << debugPrintEdges(nb_doublet_edges, cuda_reduced_M1_hits.data(),
                                  cuda_reduced_M2_hits.data()));

  // -----------------------------
  // build doublets geometric cuts
  // -----------------------------
  ScopedCudaPtr<float4> cuda_geo(nb_doublet_edges, stream);
  ScopedCudaPtr<float4> cuda_edge_slope(nb_doublet_edges, stream);

  detail::hitsGeometricCutsPacked<<<gridFor(nb_doublet_edges), block_dim, 0,
                                    stream>>>(
      cuda_geo.data(), cuda_edge_slope.data(), cuda_reduced_M1_hits.data(),
      cuda_reduced_M2_hits.data(), cuda_TThits.cuda_R(), cuda_TThits.cuda_z(),
      cuda_TThits.cuda_x(), cuda_TThits.cuda_y(), cuda_TThits.cuda_eta(),
      cuda_TThits.cuda_phi(), detail::g_pi, cfg.epsilon, nb_doublet_edges);
  ACTS_CUDA_CHECK(cudaGetLastError());
  if (cfg.debugSynchronize) {
    ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));
  }

  ScopedCudaPtr<bool> cuda_mask(nb_doublet_edges + 1, stream);
  ACTS_CUDA_CHECK(cudaMemsetAsync(
      cuda_mask.data(), 0, (nb_doublet_edges + 1) * sizeof(bool), stream));

  // -------------------------
  // loop over module triplets
  // -------------------------
  ScopedCudaPtr<int> cuda_work_to_triplet(nb_src_hits_per_triplet_sum, stream);
  constexpr int kThreadsPerTriplet = 8;
  detail::buildWorkToItem<kThreadsPerTriplet>
      <<<gridFor(static_cast<std::size_t>(nb_triplets) * kThreadsPerTriplet),
         block_dim, 0, stream>>>(cuda_work_to_triplet.data(),
                                 cuda_src_hits_per_triplet.data(), nb_triplets);
  ACTS_CUDA_CHECK(cudaGetLastError());

  // Pairs whose cuts cannot be decided from the float geometry are queued by
  // the triplet kernel and resolved by tripletPairCutsFallback
  ScopedCudaPtr<detail::TripletFallbackPair> cuda_fallback_pairs(
      detail::kTripletFallbackCapacity, stream);
  ScopedCudaPtr<int> cuda_fallback_count(1, stream);
  ACTS_CUDA_CHECK(
      cudaMemsetAsync(cuda_fallback_count.data(), 0, sizeof(int), stream));

  auto launch_triplet_cuts = [&](auto defer_fallback) {
    constexpr bool kDefer = decltype(defer_fallback)::value;
    detail::tripletPairCutsFusedM23<float, kDefer>
        <<<gridFor(nb_src_hits_per_triplet_sum), block_dim, 0, stream>>>(
            cuda_mask.data(), nb_src_hits_per_triplet_sum,
            cuda_src_hits_per_triplet.data(), cuda_work_to_triplet.data(),
            cudaModuleMapTriplet->cuda_module12_map(),
            cudaModuleMapTriplet->cuda_module23_map(), cuda_geo.data(),
            cuda_edge_slope.data(),
            cudaModuleMapTriplet->module12().cuda_z0_min(),
            cudaModuleMapTriplet->module12().cuda_phi_slope_min(),
            cudaModuleMapTriplet->module12().cuda_deta_min(),
            cudaModuleMapTriplet->module12().cuda_dphi_min(),
            cudaModuleMapTriplet->module12().cuda_z0_max(),
            cudaModuleMapTriplet->module12().cuda_phi_slope_max(),
            cudaModuleMapTriplet->module12().cuda_deta_max(),
            cudaModuleMapTriplet->module12().cuda_dphi_max(),
            cudaModuleMapTriplet->module23().cuda_z0_min(),
            cudaModuleMapTriplet->module23().cuda_phi_slope_min(),
            cudaModuleMapTriplet->module23().cuda_deta_min(),
            cudaModuleMapTriplet->module23().cuda_dphi_min(),
            cudaModuleMapTriplet->module23().cuda_z0_max(),
            cudaModuleMapTriplet->module23().cuda_phi_slope_max(),
            cudaModuleMapTriplet->module23().cuda_deta_max(),
            cudaModuleMapTriplet->module23().cuda_dphi_max(),
            cudaModuleMapTriplet->cuda_diff_dydx_min(),
            cudaModuleMapTriplet->cuda_diff_dydx_max(),
            cudaModuleMapTriplet->cuda_diff_dzdr_min(),
            cudaModuleMapTriplet->cuda_diff_dzdr_max(),
            cuda_reduced_M1_hits.data(), cuda_reduced_M2_hits.data(),
            cuda_edge_sum.data(), cuda_nb_src_hits_per_doublet.data(),
            cudaModuleMapDoublet->cuda_module1(), cuda_hit_indice,
            cuda_edge_sum_per_src_hit.data(), cuda_TThits.cuda_R(),
            cuda_TThits.cuda_z(), cuda_TThits.cuda_x(), cuda_TThits.cuda_y(),
            cuda_TThits.cuda_phi(), detail::g_pi, cfg.epsilon,
            cuda_fallback_pairs.data(), cuda_fallback_count.data(),
            detail::kTripletFallbackCapacity);
    ACTS_CUDA_CHECK(cudaGetLastError());
  };

  launch_triplet_cuts(std::true_type{});

  constexpr int kFallbackBlocks = 64;
  constexpr int kFallbackThreads = 256;
  detail::tripletPairCutsFallback<float>
      <<<kFallbackBlocks, kFallbackThreads, 0, stream>>>(
          cuda_mask.data(), cuda_fallback_pairs.data(),
          cuda_fallback_count.data(), detail::kTripletFallbackCapacity,
          cuda_geo.data(), cudaModuleMapTriplet->module12().cuda_z0_min(),
          cudaModuleMapTriplet->module12().cuda_phi_slope_min(),
          cudaModuleMapTriplet->module12().cuda_deta_min(),
          cudaModuleMapTriplet->module12().cuda_dphi_min(),
          cudaModuleMapTriplet->module12().cuda_z0_max(),
          cudaModuleMapTriplet->module12().cuda_phi_slope_max(),
          cudaModuleMapTriplet->module12().cuda_deta_max(),
          cudaModuleMapTriplet->module12().cuda_dphi_max(),
          cudaModuleMapTriplet->module23().cuda_z0_min(),
          cudaModuleMapTriplet->module23().cuda_phi_slope_min(),
          cudaModuleMapTriplet->module23().cuda_deta_min(),
          cudaModuleMapTriplet->module23().cuda_dphi_min(),
          cudaModuleMapTriplet->module23().cuda_z0_max(),
          cudaModuleMapTriplet->module23().cuda_phi_slope_max(),
          cudaModuleMapTriplet->module23().cuda_deta_max(),
          cudaModuleMapTriplet->module23().cuda_dphi_max(),
          cudaModuleMapTriplet->cuda_diff_dydx_min(),
          cudaModuleMapTriplet->cuda_diff_dydx_max(),
          cudaModuleMapTriplet->cuda_diff_dzdr_min(),
          cudaModuleMapTriplet->cuda_diff_dzdr_max(),
          cuda_reduced_M1_hits.data(), cuda_reduced_M2_hits.data(),
          cuda_TThits.cuda_R(), cuda_TThits.cuda_z(), cuda_TThits.cuda_x(),
          cuda_TThits.cuda_y(), cuda_TThits.cuda_phi(), detail::g_pi,
          cfg.epsilon);
  ACTS_CUDA_CHECK(cudaGetLastError());

  //------------------------
  // edges reduction (sync 3)
  //------------------------
  ScopedCudaPtr<int> cuda_mask_sum(nb_doublet_edges + 1, stream);
  auto scan_mask = [&]() {
    exclusiveSum(
        thrust::make_transform_iterator(cuda_mask.data(), CastBoolToInt{}),
        cuda_mask_sum.data(), nb_doublet_edges + 1, stream);
  };
  scan_mask();

  int nb_graph_edges{};
  int nb_fallback_pairs{};
  ACTS_CUDA_CHECK(cudaMemcpyAsync(&nb_graph_edges,
                                  &cuda_mask_sum.data()[nb_doublet_edges],
                                  sizeof(int), cudaMemcpyDeviceToHost, stream));
  ACTS_CUDA_CHECK(cudaMemcpyAsync(&nb_fallback_pairs,
                                  cuda_fallback_count.data(), sizeof(int),
                                  cudaMemcpyDeviceToHost, stream));
  ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));
  ACTS_DEBUG("nb_fallback_pairs: " << nb_fallback_pairs);

  if (nb_fallback_pairs > detail::kTripletFallbackCapacity) {
    // Some queued pairs were dropped. All tags set so far are correct (they
    // are only set for pairs that pass), so rerunning the triplet cuts with
    // inline fallback completes the mask.
    ACTS_WARNING("Triplet fallback queue overflow ("
                 << nb_fallback_pairs << " > "
                 << detail::kTripletFallbackCapacity
                 << "), rerunning triplet cuts with inline fallback");
    launch_triplet_cuts(std::false_type{});
    scan_mask();
    ACTS_CUDA_CHECK(cudaMemcpyAsync(
        &nb_graph_edges, &cuda_mask_sum.data()[nb_doublet_edges], sizeof(int),
        cudaMemcpyDeviceToHost, stream));
    ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));
  }

  ACTS_DEBUG("nb_graph_edges: " << nb_graph_edges);
  if (nb_graph_edges == 0) {
    throw NoEdgesError{};
  }

  // Compact the surviving edges straight into the output tensors and compute
  // the edge features in the same pass
  const auto nEdges = static_cast<std::size_t>(nb_graph_edges);
  auto edgeIndex = Tensor<std::int64_t>::Create({2, nEdges}, execContext);
  auto edgeFeatures = Tensor<float>::Create(
      {nEdges, static_cast<std::size_t>(detail::g_nEdgeFeatures)}, execContext);

  detail::compactEdgesAndMakeFeatures<<<gridFor(nb_doublet_edges), block_dim, 0,
                                        stream>>>(
      nb_doublet_edges, cuda_mask.data(), cuda_mask_sum.data(),
      cuda_reduced_M1_hits.data(), cuda_reduced_M2_hits.data(), nEdges,
      nNodeFeatures, cudaNodeFeatures, edgeIndex.data(), edgeFeatures.data());
  ACTS_CUDA_CHECK(cudaGetLastError());

  ACTS_VERBOSE(
      "First 10 graph edges:\n"
      << debugPrintEdges(nEdges, edgeIndex.data(), edgeIndex.data() + nEdges));

  return {std::move(edgeIndex), std::move(edgeFeatures)};
}

}  // namespace ActsPlugins

// clang-format off
// EdgeLayerConnector implementation included here due to ODR violations in ModuleMapGraph
#include "EdgeLayerConnector.cu"
// clang-format on
