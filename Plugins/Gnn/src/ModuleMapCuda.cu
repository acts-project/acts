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
#include "ActsPlugins/Gnn/detail/DeviceMemory.cuh"
#include "ActsPlugins/Gnn/detail/ModuleMapUtils.cuh"

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <map>
#include <numbers>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <MMG/CUDA_graph_creator>
#include <MMG/CUDA_module_map_doublet>
#include <MMG/CUDA_module_map_triplet>
#include <thrust/execution_policy.h>
#include <thrust/functional.h>
#include <thrust/scan.h>
#include <thrust/sort.h>
#include <thrust/transform_scan.h>
#include <vecmem/containers/data/vector_buffer.hpp>
#include <vecmem/containers/data/vector_view.hpp>
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/utils/cuda/copy.hpp>

using Clock = std::chrono::high_resolution_clock;

namespace {

template <typename T>
struct CUDA_hit_data {
  std::size_t m_size;
  T *m_cuda_x;
  T *m_cuda_y;
  T *m_cuda_R;
  T *m_cuda_phi;
  T *m_cuda_z;
  T *m_cuda_eta;

  std::size_t size() { return m_size; }
  T *cuda_x() { return m_cuda_x; }
  T *cuda_y() { return m_cuda_y; }
  T *cuda_z() { return m_cuda_z; }
  T *cuda_R() { return m_cuda_R; }
  T *cuda_phi() { return m_cuda_phi; }
  T *cuda_eta() { return m_cuda_eta; }
};

/// Source nodes followed by target nodes
struct CUDA_edge_data {
  std::size_t nEdges = 0;
  vecmem::data::vector_buffer<int> edges;
};

struct CastBoolToInt {
  int __device__ operator()(bool b) { return static_cast<int>(b); }
};

std::string debugPrintEdges(std::size_t nbEdges, const int *cudaSrc,
                            const int *cudaDst,
                            ActsPlugins::detail::DeviceMemory &mem) {
  std::stringstream ss;
  if (nbEdges == 0) {
    return "zero edges remained";
  }
  nbEdges = std::min(10ul, nbEdges);
  std::vector<int> src(nbEdges), dst(nbEdges);
  const auto n = static_cast<unsigned int>(nbEdges);
  mem.copy()(vecmem::data::vector_view<const int>(n, cudaSrc),
             vecmem::get_data(src), vecmem::copy::type::device_to_host)
      ->wait();
  mem.copy()(vecmem::data::vector_view<const int>(n, cudaDst),
             vecmem::get_data(dst), vecmem::copy::type::device_to_host)
      ->wait();
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

  /// Memory of the module map lookup table
  vecmem::cuda::device_memory_resource moduleMapMemory;
  vecmem::data::vector_buffer<std::uint64_t> cudaModuleMapKeys;
  vecmem::data::vector_buffer<int> cudaModuleMapVals;
  std::size_t cudaModuleMapSize{};

  CUDA_edge_data makeEdges(CUDA_hit_data<float> cuda_TThits,
                           int *cuda_hit_indice, detail::DeviceMemory &mem,
                           const ModuleMapCuda::Config &cfg,
                           const Logger &logger) const;
};

ModuleMapCuda::ModuleMapCuda(const Config &cfg,
                             std::unique_ptr<const Logger> logger_)
    : m_impl(std::make_unique<Impl>()),
      m_cfg(cfg),
      m_logger(std::move(logger_)) {
  module_map_triplet<float> moduleMapCpu;
  moduleMapCpu.read_TTree(cfg.moduleMapPath.c_str());
  if (!moduleMapCpu) {
    throw std::runtime_error("Cannot retrieve ModuleMap from " +
                             cfg.moduleMapPath);
  }

  ACTS_DEBUG("ModuleMap GPU block dim: " << m_cfg.gpuBlocks);

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
  const auto size = static_cast<vecmem::data::vector_buffer<int>::size_type>(
      m_impl->cudaModuleMapSize);
  m_impl->cudaModuleMapKeys =
      vecmem::data::vector_buffer<std::uint64_t>(size, m_impl->moduleMapMemory);
  m_impl->cudaModuleMapVals =
      vecmem::data::vector_buffer<int>(size, m_impl->moduleMapMemory);
  vecmem::cuda::copy copy;
  copy.setup(m_impl->cudaModuleMapKeys)->ignore();
  copy.setup(m_impl->cudaModuleMapVals)->ignore();
  copy(vecmem::get_data(keys), m_impl->cudaModuleMapKeys,
       vecmem::copy::type::host_to_device)
      ->wait();
  copy(vecmem::get_data(vals), m_impl->cudaModuleMapVals,
       vecmem::copy::type::host_to_device)
      ->wait();
}

ModuleMapCuda::~ModuleMapCuda() = default;

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

  detail::DeviceMemory mem(execContext);
  const cudaStream_t stream = mem.stream();
  ACTS_DEBUG("Got stream " << stream);

  /////////////////////////
  // Prepare input data
  ////////////////////////

  // Full node features to device

  auto nodeFeatures = Tensor<float>::Create({nHits, nFeatures}, execContext);
  float *cudaNodeFeaturePtr = nodeFeatures.data();
  mem.copy()(
         vecmem::get_data(features),
         vecmem::data::vector_view<float>(
             static_cast<unsigned int>(features.size()), cudaNodeFeaturePtr),
         vecmem::copy::type::host_to_device)
      ->ignore();

  // Module IDs to device
  vecmem::data::vector_buffer<std::uint64_t> cudaModuleIds(
      static_cast<unsigned int>(nHits), mem.mr());
  mem.copy().setup(cudaModuleIds)->ignore();
  mem.copy()(vecmem::get_data(moduleIds), cudaModuleIds,
             vecmem::copy::type::host_to_device)
      ->ignore();

  // Allocate memory for transposed node features that are needed for the
  // module map kernels in one block
  vecmem::data::vector_buffer<float> cudaNodeFeaturesTransposed(
      static_cast<unsigned int>(6 * nHits), mem.mr());
  mem.copy().setup(cudaNodeFeaturesTransposed)->ignore();

  CUDA_hit_data<float> inputData{};
  inputData.m_size = nHits;
  inputData.m_cuda_R = cudaNodeFeaturesTransposed.ptr() + 0 * nHits;
  inputData.m_cuda_phi = cudaNodeFeaturesTransposed.ptr() + 1 * nHits;
  inputData.m_cuda_z = cudaNodeFeaturesTransposed.ptr() + 2 * nHits;
  inputData.m_cuda_x = cudaNodeFeaturesTransposed.ptr() + 3 * nHits;
  inputData.m_cuda_y = cudaNodeFeaturesTransposed.ptr() + 4 * nHits;
  inputData.m_cuda_eta = cudaNodeFeaturesTransposed.ptr() + 5 * nHits;

  // Node features for module map graph
  constexpr std::size_t rOffset = 0;
  constexpr std::size_t phiOffset = 1;
  constexpr std::size_t zOffset = 2;
  constexpr std::size_t etaOffset = 3;

  // Column copies, which vecmem does not provide
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
  vecmem::data::vector_buffer<int> cudaNbHits(
      static_cast<unsigned int>(m_impl->cudaModuleMapSize + 1), mem.mr());
  mem.copy().setup(cudaNbHits)->ignore();
  mem.copy().memset(cudaNbHits, 0)->ignore();

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

  detail::mapModuleIdsToNbHits<<<gridDimHits, blockDim, 0, stream>>>(
      cudaNbHits.ptr(), nHits, cudaModuleIds.ptr(), m_impl->cudaModuleMapSize,
      m_impl->cudaModuleMapKeys.ptr(), m_impl->cudaModuleMapVals.ptr());
  ACTS_CUDA_CHECK(cudaGetLastError());

  thrust::exclusive_scan(mem.policy(), cudaNbHits.ptr(),
                         cudaNbHits.ptr() + m_impl->cudaModuleMapSize + 1,
                         cudaNbHits.ptr());
  ACTS_CUDA_CHECK(cudaGetLastError());
  int *cudaHitIndice = cudaNbHits.ptr();

  ///////////////////////////////////
  // Perform module map inference
  ////////////////////////////////////

  mem.synchronize();
  auto t1 = std::chrono::high_resolution_clock::now();

  const auto edgeData =
      m_impl->makeEdges(inputData, cudaHitIndice, mem, m_cfg, logger());
  ACTS_CUDA_CHECK(cudaGetLastError());
  mem.synchronize();

  auto t2 = std::chrono::high_resolution_clock::now();

  if (edgeData.nEdges == 0) {
    throw NoEdgesError{};
  }

  dim3 gridDimEdges = (edgeData.nEdges + blockDim.x - 1) / blockDim.x;
  ACTS_DEBUG("gridDimEdges: " << gridDimEdges.x
                              << ", blockDim: " << blockDim.x);

  // Make edge features
  auto edgeFeatures = Tensor<float>::Create({edgeData.nEdges, 6}, execContext);
  mem.synchronize();

  detail::makeEdgeFeatures<<<gridDimEdges, blockDim, 0, stream>>>(
      edgeData.nEdges, edgeData.edges.ptr(),
      edgeData.edges.ptr() + edgeData.nEdges, nFeatures, cudaNodeFeaturePtr,
      edgeFeatures.data());
  ACTS_CUDA_CHECK(cudaGetLastError());
  mem.synchronize();

  auto edgeIndex =
      Tensor<std::int64_t>::Create({2, edgeData.nEdges}, execContext);
  thrust::transform(mem.policy(), edgeData.edges.ptr(),
                    edgeData.edges.ptr() + 2 * edgeData.nEdges,
                    edgeIndex.data(),
                    [] __device__(int i) -> std::int64_t { return i; });

  ACTS_CUDA_CHECK(cudaGetLastError());
  mem.synchronize();
  auto t3 = std::chrono::high_resolution_clock::now();

  auto ms = [](auto a, auto b) {
    return std::chrono::duration<double, std::milli>(b - a).count();
  };
  ACTS_DEBUG("Preparation: " << ms(t0, t1));
  ACTS_DEBUG("Inference: " << ms(t1, t2));
  ACTS_DEBUG("Postprocessing: " << ms(t2, t3));

  return {std::move(nodeFeatures),
          std::move(edgeIndex),
          std::move(edgeFeatures),
          {}};
}

CUDA_edge_data ModuleMapCuda::Impl::makeEdges(CUDA_hit_data<float> cuda_TThits,
                                              int *cuda_hit_indice,
                                              detail::DeviceMemory &mem,
                                              const ModuleMapCuda::Config &cfg,
                                              const Logger &logger) const {
  const cudaStream_t stream = mem.stream();
  const dim3 block_dim = cfg.gpuBlocks;
  // ----------------------------------
  // memory allocation for hits + edges
  // ----------------------------------
  const int nb_doublets = cudaModuleMapDoublet->size();
  ACTS_DEBUG("nb doublets " << nb_doublets);
  dim3 grid_dim = ((nb_doublets + block_dim.x - 1) / block_dim.x);

  // hit buffer pointers need to be visible outside the scope
  // Should be definitively filled after the if block
  vecmem::data::vector_buffer<int> cuda_reduced_M1_hits, cuda_reduced_M2_hits;

  vecmem::data::vector_buffer<int> cuda_edge_sum(
      static_cast<unsigned int>(nb_doublets + 1), mem.mr());
  mem.copy().setup(cuda_edge_sum)->ignore();

  int nb_doublet_edges{};

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

  vecmem::data::vector_buffer<int> cuda_nb_src_hits_per_doublet(
      static_cast<unsigned int>(nb_doublets + 1), mem.mr());
  mem.copy().setup(cuda_nb_src_hits_per_doublet)->ignore();

  count_source_hits_per_doublet<<<grid_dim, block_dim, 0, stream>>>(
      cuda_nb_src_hits_per_doublet.ptr(), cudaModuleMapDoublet->cuda_module1(),
      cuda_hit_indice, nb_doublets);
  ACTS_CUDA_CHECK(cudaGetLastError());

  thrust::exclusive_scan(mem.policy(), cuda_nb_src_hits_per_doublet.ptr(),
                         cuda_nb_src_hits_per_doublet.ptr() + nb_doublets + 1,
                         cuda_nb_src_hits_per_doublet.ptr());

  int sum_nb_src_hits_per_doublet{};
  mem.copy()(vecmem::data::vector_view<const int>(
                 1u, cuda_nb_src_hits_per_doublet.ptr() + nb_doublets),
             vecmem::data::vector_view<int>(1u, &sum_nb_src_hits_per_doublet),
             vecmem::copy::type::device_to_host)
      ->wait();
  ACTS_DEBUG("sum_nb_hits_per_doublet: " << sum_nb_src_hits_per_doublet);

  if (sum_nb_src_hits_per_doublet == 0) {
    throw NoEdgesError{};
  }

  vecmem::data::vector_buffer<int> cuda_edge_sum_per_src_hit(
      static_cast<unsigned int>(sum_nb_src_hits_per_doublet + 1), mem.mr());
  mem.copy().setup(cuda_edge_sum_per_src_hit)->ignore();

  dim3 grid_dim_shpd =
      ((sum_nb_src_hits_per_doublet + block_dim.x - 1) / block_dim.x);
  count_doublet_edges<float><<<grid_dim_shpd, block_dim, 0, stream>>>(
      cuda_edge_sum_per_src_hit.ptr(), nb_doublets,
      cuda_nb_src_hits_per_doublet.ptr(), cudaModuleMapDoublet->cuda_module1(),
      cudaModuleMapDoublet->cuda_module2(), cuda_TThits.cuda_R(),
      cuda_TThits.cuda_z(), cuda_TThits.cuda_eta(), cuda_TThits.cuda_phi(),
      cudaModuleMapDoublet->cuda_z0_min(),
      cudaModuleMapDoublet->cuda_deta_min(),
      cudaModuleMapDoublet->cuda_phi_slope_min(),
      cudaModuleMapDoublet->cuda_dphi_min(),
      cudaModuleMapDoublet->cuda_z0_max(),
      cudaModuleMapDoublet->cuda_deta_max(),
      cudaModuleMapDoublet->cuda_phi_slope_max(),
      cudaModuleMapDoublet->cuda_dphi_max(), cuda_hit_indice,
      std::numbers::pi_v<float>, cfg.epsilon, sum_nb_src_hits_per_doublet);
  ACTS_CUDA_CHECK(cudaGetLastError());

  thrust::exclusive_scan(
      mem.policy(), cuda_edge_sum_per_src_hit.ptr(),
      cuda_edge_sum_per_src_hit.ptr() + sum_nb_src_hits_per_doublet + 1,
      cuda_edge_sum_per_src_hit.ptr());

  mem.copy()(
         vecmem::data::vector_view<const int>(
             1u, cuda_edge_sum_per_src_hit.ptr() + sum_nb_src_hits_per_doublet),
         vecmem::data::vector_view<int>(1u, &nb_doublet_edges),
         vecmem::copy::type::device_to_host)
      ->wait();
  ACTS_DEBUG("nb_doublet_edges: " << nb_doublet_edges);
  ACTS_DEBUG("Allocate " << (2ul * nb_doublet_edges * sizeof(int)) * 1.0e-6
                         << " MB for edges");
  cuda_reduced_M1_hits = vecmem::data::vector_buffer<int>(
      static_cast<unsigned int>(nb_doublet_edges), mem.mr());
  mem.copy().setup(cuda_reduced_M1_hits)->ignore();
  cuda_reduced_M2_hits = vecmem::data::vector_buffer<int>(
      static_cast<unsigned int>(nb_doublet_edges), mem.mr());
  mem.copy().setup(cuda_reduced_M2_hits)->ignore();

  build_doublet_edges<float><<<grid_dim_shpd, block_dim, 0, stream>>>(
      cuda_reduced_M1_hits.ptr(), cuda_reduced_M2_hits.ptr(), nb_doublets,
      cuda_nb_src_hits_per_doublet.ptr(), cuda_hit_indice,
      cuda_edge_sum_per_src_hit.ptr(), cudaModuleMapDoublet->cuda_module1(),
      cudaModuleMapDoublet->cuda_module2(), cuda_TThits.cuda_R(),
      cuda_TThits.cuda_z(), cuda_TThits.cuda_eta(), cuda_TThits.cuda_phi(),
      cudaModuleMapDoublet->cuda_z0_min(),
      cudaModuleMapDoublet->cuda_deta_min(),
      cudaModuleMapDoublet->cuda_phi_slope_min(),
      cudaModuleMapDoublet->cuda_dphi_min(),
      cudaModuleMapDoublet->cuda_z0_max(),
      cudaModuleMapDoublet->cuda_deta_max(),
      cudaModuleMapDoublet->cuda_phi_slope_max(),
      cudaModuleMapDoublet->cuda_dphi_max(), std::numbers::pi_v<float>,
      cfg.epsilon, sum_nb_src_hits_per_doublet);
  ACTS_CUDA_CHECK(cudaGetLastError());

  doublet_edge_sum<<<grid_dim, block_dim, 0, stream>>>(
      cuda_edge_sum.ptr(), cuda_nb_src_hits_per_doublet.ptr(),
      cuda_edge_sum_per_src_hit.ptr(), nb_doublets);
  ACTS_CUDA_CHECK(cudaGetLastError());

  ACTS_VERBOSE("First 10 doublet edges:\n"
               << debugPrintEdges(nb_doublet_edges, cuda_reduced_M1_hits.ptr(),
                                  cuda_reduced_M2_hits.ptr(), mem));
  if (nb_doublet_edges == 0) {
    throw NoEdgesError{};
  }

  // -----------------------------
  // build doublets geometric cuts
  // -----------------------------

  vecmem::data::vector_buffer<float> cuda_z0(
      static_cast<unsigned int>(nb_doublet_edges), mem.mr());
  mem.copy().setup(cuda_z0)->ignore();
  vecmem::data::vector_buffer<float> cuda_phi_slope(
      static_cast<unsigned int>(nb_doublet_edges), mem.mr());
  mem.copy().setup(cuda_phi_slope)->ignore();
  vecmem::data::vector_buffer<float> cuda_deta(
      static_cast<unsigned int>(nb_doublet_edges), mem.mr());
  mem.copy().setup(cuda_deta)->ignore();
  vecmem::data::vector_buffer<float> cuda_dphi(
      static_cast<unsigned int>(nb_doublet_edges), mem.mr());
  mem.copy().setup(cuda_dphi)->ignore();
  grid_dim = ((nb_doublet_edges + block_dim.x - 1) / block_dim.x);
  hits_geometric_cuts<<<grid_dim, block_dim, 0, stream>>>(
      cuda_z0.ptr(), cuda_phi_slope.ptr(), cuda_deta.ptr(), cuda_dphi.ptr(),
      cuda_reduced_M1_hits.ptr(), cuda_reduced_M2_hits.ptr(),
      cuda_TThits.cuda_R(), cuda_TThits.cuda_z(), cuda_TThits.cuda_eta(),
      cuda_TThits.cuda_phi(), std::numbers::pi_v<float>, cfg.epsilon,
      nb_doublet_edges);
  ACTS_CUDA_CHECK(cudaGetLastError());
  mem.synchronize();

  // vecmem has no bool containers, the flags are bytes for the MMG kernels
  vecmem::data::vector_buffer<unsigned char> cuda_mask_bytes(
      static_cast<unsigned int>(nb_doublet_edges + 1), mem.mr());
  mem.copy().setup(cuda_mask_bytes)->ignore();
  mem.copy().memset(cuda_mask_bytes, 0)->ignore();
  bool *cuda_mask = reinterpret_cast<bool *>(cuda_mask_bytes.ptr());
  mem.synchronize();

  // -------------------------
  // loop over module triplets
  // -------------------------
  int nb_triplets = cudaModuleMapTriplet->size();
  grid_dim = ((nb_triplets + block_dim.x - 1) / block_dim.x);

  vecmem::data::vector_buffer<unsigned char> cuda_vertices_bytes(
      static_cast<unsigned int>(cuda_TThits.size()), mem.mr());
  mem.copy().setup(cuda_vertices_bytes)->ignore();
  mem.copy().memset(cuda_vertices_bytes, 0)->ignore();
  bool *cuda_vertices = reinterpret_cast<bool *>(cuda_vertices_bytes.ptr());
  mem.synchronize();

  // Allocate memory for the number of hits per triplet
  vecmem::data::vector_buffer<int> cuda_src_hits_per_triplet(
      static_cast<unsigned int>(nb_triplets + 1), mem.mr());
  mem.copy().setup(cuda_src_hits_per_triplet)->ignore();

  count_triplet_hits<<<grid_dim, block_dim, 0, stream>>>(
      cuda_src_hits_per_triplet.ptr(),
      cudaModuleMapTriplet->cuda_module12_map(),
      cudaModuleMapTriplet->cuda_module23_map(), cuda_edge_sum.ptr(),
      nb_triplets);
  ACTS_CUDA_CHECK(cudaGetLastError());

  // Perform prefix sum to get the offset for each triplet
  thrust::exclusive_scan(mem.policy(), cuda_src_hits_per_triplet.ptr(),
                         cuda_src_hits_per_triplet.ptr() + nb_triplets + 1,
                         cuda_src_hits_per_triplet.ptr());

  int nb_src_hits_per_triplet_sum{};
  mem.copy()(vecmem::data::vector_view<const int>(
                 1u, cuda_src_hits_per_triplet.ptr() + nb_triplets),
             vecmem::data::vector_view<int>(1u, &nb_src_hits_per_triplet_sum),
             vecmem::copy::type::device_to_host)
      ->wait();
  ACTS_DEBUG("nb_src_hits_per_triplet_sum: " << nb_src_hits_per_triplet_sum);
  if (nb_src_hits_per_triplet_sum == 0) {
    throw NoEdgesError{};
  }

  dim3 grid_dim_shpt =
      ((nb_src_hits_per_triplet_sum + block_dim.x - 1) / block_dim.x);

  triplet_cuts<float><<<grid_dim_shpt, block_dim, 0, stream>>>(
      cuda_vertices, cuda_mask, nb_triplets, cuda_src_hits_per_triplet.ptr(),
      cudaModuleMapTriplet->cuda_module12_map(),
      cudaModuleMapTriplet->cuda_module23_map(), cuda_TThits.cuda_x(),
      cuda_TThits.cuda_y(), cuda_TThits.cuda_z(), cuda_TThits.cuda_R(),
      cuda_z0.ptr(), cuda_phi_slope.ptr(), cuda_deta.ptr(), cuda_dphi.ptr(),
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
      cudaModuleMapTriplet->cuda_diff_dzdr_max(), cuda_reduced_M1_hits.ptr(),
      cuda_reduced_M2_hits.ptr(), cuda_edge_sum.ptr(),
      std::numbers::pi_v<float>, cfg.epsilon, nb_src_hits_per_triplet_sum);
  ACTS_CUDA_CHECK(cudaGetLastError());

  //----------------
  // edges reduction
  //----------------
  vecmem::data::vector_buffer<int> cuda_mask_sum(
      static_cast<unsigned int>(nb_doublet_edges + 1), mem.mr());
  mem.copy().setup(cuda_mask_sum)->ignore();
  mem.synchronize();
  thrust::transform_exclusive_scan(
      mem.policy(), cuda_mask, cuda_mask + (nb_doublet_edges + 1),
      cuda_mask_sum.ptr(), CastBoolToInt{}, 0, thrust::plus<int>());

  int nb_graph_edges{};
  mem.copy()(vecmem::data::vector_view<const int>(
                 1u, cuda_mask_sum.ptr() + nb_doublet_edges),
             vecmem::data::vector_view<int>(1u, &nb_graph_edges),
             vecmem::copy::type::device_to_host)
      ->wait();
  ACTS_DEBUG("nb_graph_edges: " << nb_graph_edges);
  if (nb_graph_edges == 0) {
    throw NoEdgesError{};
  }

  CUDA_edge_data edge_data;
  edge_data.nEdges = nb_graph_edges;
  edge_data.edges = vecmem::data::vector_buffer<int>(
      static_cast<unsigned int>(2 * nb_graph_edges), mem.mr());
  mem.copy().setup(edge_data.edges)->ignore();
  int *cuda_graph_M1_hits = edge_data.edges.ptr();
  int *cuda_graph_M2_hits = edge_data.edges.ptr() + nb_graph_edges;

  grid_dim = ((nb_doublet_edges + block_dim.x - 1) / block_dim.x);
  mem.synchronize();
  compact_stream<<<grid_dim, block_dim, 0, stream>>>(
      cuda_graph_M1_hits, cuda_reduced_M1_hits.ptr(), cuda_mask,
      cuda_mask_sum.ptr(), nb_doublet_edges);
  ACTS_CUDA_CHECK(cudaGetLastError());
  compact_stream<<<grid_dim, block_dim, 0, stream>>>(
      cuda_graph_M2_hits, cuda_reduced_M2_hits.ptr(), cuda_mask,
      cuda_mask_sum.ptr(), nb_doublet_edges);
  ACTS_CUDA_CHECK(cudaGetLastError());

  ACTS_VERBOSE("First 10 doublet edges:\n"
               << debugPrintEdges(nb_graph_edges, cuda_graph_M1_hits,
                                  cuda_graph_M2_hits, mem));

  return edge_data;
}

}  // namespace ActsPlugins

// clang-format off
// EdgeLayerConnector implementation included here due to ODR violations in ModuleMapGraph
#include "EdgeLayerConnector.cu"
// clang-format on
