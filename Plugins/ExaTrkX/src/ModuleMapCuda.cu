// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/ModuleMapCuda.hpp"
#include "Acts/Plugins/ExaTrkX/detail/CudaUtils.cuh"
#include "Acts/Plugins/ExaTrkX/detail/TensorVectorConversion.hpp"

#include <CUDA_graph_creator>
#include <CUDA_module_map_doublet>
#include <CUDA_module_map_triplet>
#include <TTree_hits>

using namespace torch::indexing;

namespace Acts {

ModuleMapCuda::ModuleMapCuda(const Config &cfg,
                             std::unique_ptr<const Acts::Logger> logger_)
    : m_logger(std::move(logger_)) {
  module_map_triplet<float> moduleMapCpu;
  moduleMapCpu.read_TTree(cfg.moduleMapPath.c_str());
  if (!moduleMapCpu) {
    throw std::runtime_error("Cannot retrieve ModuleMap from " +
                             cfg.moduleMapPath);
  }

  m_cudaModuleMapDoublet =
      std::make_unique<CUDA_module_map_doublet<float>>(moduleMapCpu);
  m_cudaModuleMapDoublet->HostToDevice();
  m_cudaModuleMapTriplet =
      std::make_unique<CUDA_module_map_triplet<float>>(moduleMapCpu);
  m_cudaModuleMapTriplet->HostToDevice();

  ACTS_DEBUG("# of modules = " << moduleMapCpu.module_map().size());
}

ModuleMapCuda::~ModuleMapCuda() {}

namespace {
std::vector<int> makeHitIndexPrefixSum(
    const std::vector<std::uint64_t> &moduleIds,
    const std::multimap<std::uint64_t, int> module_map) {
  std::vector<int> hit_indice;
  std::vector<int> nb_hits(module_map.size(), 0);
  std::vector<int> hits_bool_mask(moduleIds.size(), 0);

  for (auto i = 0ul; i < moduleIds.size(); ++i) {
    auto it = module_map.find(moduleIds.at(i));
    if (it != module_map.end() && hits_bool_mask.at(i) == 0) {
      hits_bool_mask.at(i) = 1;
      nb_hits[it->second] += 1;
    }
  }

  hit_indice.push_back(0);
  for (std::size_t i = 0; i < nb_hits.size(); i++) {
    hit_indice.push_back(hit_indice[i] + nb_hits[i]);
  }

  assert(!std::all_of(hit_indice.begin(), hit_indice.end(),
                      [](auto v) { return v == 0; }));
  return hit_indice;
}
}  // namespace

std::tuple<std::any, std::any, std::any> ModuleMapCuda::operator()(
    std::vector<float> &inputValues, std::size_t numNodes,
    const std::vector<std::uint64_t> &moduleIds, torch::Device device) {
  const auto nHits = moduleIds.size();
  const auto nFeatures = inputValues.size() / moduleIds.size();
  auto &features = inputValues;

  const dim3 blockDim = m_cfg.gpuBlocks;
  const dim3 gridDimHits = (nHits + blockDim.x - 1) / blockDim.x;
  ACTS_VERBOSE("gridDimHits: " << gridDimHits.x
                               << ", blockDim: " << blockDim.x);

  /////////////////////////
  // Prepare input data
  ////////////////////////
  constexpr std::size_t rOffset = 0;
  constexpr std::size_t phiOffset = 1;
  constexpr std::size_t zOffset = 2;
  constexpr std::size_t etaOffset = 3;

  // const auto srcStride = sizeof(float) * nFeatures;
  // const auto dstStride = sizeof(float);  // contiguous in destination
  // const auto width = sizeof(float);      // only copy 1 column
  // const auto height = nHits;

  std::vector<float> hostR(nHits), hostPhi(nHits), hostZ(nHits), hostEta(nHits);
  for (auto i = 0ul; i < nHits; ++i) {
    hostR.at(i) = features.at(i * nFeatures + rOffset);
    hostPhi.at(i) = features.at(i * nFeatures + phiOffset);
    hostZ.at(i) = features.at(i * nFeatures + zOffset);
    hostEta.at(i) = features.at(i * nFeatures + etaOffset);
  }

  Acts::detail::CUDA_hit_data<float> inputData;
  inputData.m_size = nHits;

  CUDA_CHECK(cudaMalloc(&inputData.m_cuda_R, nHits * sizeof(float)));
  CUDA_CHECK(cudaMemcpy(inputData.m_cuda_R, hostR.data(), nHits * sizeof(float),
                        cudaMemcpyHostToDevice));
  // CUDA_CHECK(cudaMemcpy2D(inputData.cuda_R, dstStride,
  //                         features.data() + rOffset, srcStride, width,
  //                         height, cudaMemcpyHostToDevice));

  CUDA_CHECK(cudaMalloc(&inputData.m_cuda_phi, nHits * sizeof(float)));
  CUDA_CHECK(cudaMemcpy(inputData.m_cuda_phi, hostPhi.data(),
                        nHits * sizeof(float), cudaMemcpyHostToDevice));
  // CUDA_CHECK(cudaMemcpy2D(inputData.cuda_phi, dstStride,
  //                         features.data() + phiOffset, srcStride, width,
  //                         height, cudaMemcpyHostToDevice));

  CUDA_CHECK(cudaMalloc(&inputData.m_cuda_z, nHits * sizeof(float)));
  CUDA_CHECK(cudaMemcpy(inputData.m_cuda_z, hostZ.data(), nHits * sizeof(float),
                        cudaMemcpyHostToDevice));
  // CUDA_CHECK(cudaMemcpy2D(inputData.cuda_z, dstStride,
  //                         features.data() + zOffset, srcStride, width,
  //                         height, cudaMemcpyHostToDevice));

  CUDA_CHECK(cudaMalloc(&inputData.m_cuda_eta, nHits * sizeof(float)));
  CUDA_CHECK(cudaMemcpy(inputData.m_cuda_eta, hostEta.data(),
                        nHits * sizeof(float), cudaMemcpyHostToDevice));
  // CUDA_CHECK(cudaMemcpy2D(inputData.cuda_eta, dstStride,
  //                         features.data() + etaOffset, srcStride, width,
  //                         height, cudaMemcpyHostToDevice));

  CUDA_CHECK(cudaDeviceSynchronize());
  Acts::detail::rescaleFeature<<<gridDimHits, blockDim>>>(
      nHits, inputData.m_cuda_z, 1000.f);
  CUDA_CHECK(cudaGetLastError());
  Acts::detail::rescaleFeature<<<gridDimHits, blockDim>>>(
      nHits, inputData.m_cuda_R, 1000.f);
  CUDA_CHECK(cudaGetLastError());
  Acts::detail::rescaleFeature<<<gridDimHits, blockDim>>>(
      nHits, inputData.m_cuda_phi, 3.14159f);
  CUDA_CHECK(cudaGetLastError());

  CUDA_CHECK(cudaMalloc(&inputData.m_cuda_x, nHits * sizeof(float)));
  CUDA_CHECK(cudaMalloc(&inputData.m_cuda_y, nHits * sizeof(float)));

  Acts::detail::computeXandY<<<gridDimHits, blockDim>>>(
      nHits, inputData.m_cuda_x, inputData.m_cuda_y, inputData.m_cuda_R,
      inputData.m_cuda_phi);
  CUDA_CHECK(cudaGetLastError());

  CUDA_CHECK(
      cudaMalloc(&inputData.m_cuda_hit_id, nHits * sizeof(std::uint64_t)));
  Acts::detail::setHitId<<<gridDimHits, blockDim>>>(nHits,
                                                    inputData.m_cuda_hit_id);
  CUDA_CHECK(cudaGetLastError());

  const auto hit_indice =
      makeHitIndexPrefixSum(moduleIds, m_cudaModuleMapDoublet->module_map());
  int *cuda_hit_indice = nullptr;
  CUDA_CHECK(cudaMalloc(&cuda_hit_indice, hit_indice.size() * sizeof(int)));
  CUDA_CHECK(cudaMemcpy(cuda_hit_indice, hit_indice.data(),
                        hit_indice.size() * sizeof(int),
                        cudaMemcpyHostToDevice));

  ///////////////////////////////////
  // Perform module map inference
  ////////////////////////////////////
  cudaDeviceSynchronize();
  const auto [edgeData, hitData] = makeEdges(inputData, cuda_hit_indice);
  CUDA_CHECK(cudaGetLastError());
  CUDA_CHECK(cudaDeviceSynchronize());

  dim3 gridDimEdges = (edgeData.size + blockDim.x - 1) / blockDim.x;

  Acts::detail::remapEdges<<<gridDimEdges, blockDim>>>(
      edgeData.size, edgeData.cuda_graph_M1_hits, edgeData.cuda_graph_M2_hits,
      hitData.m_cuda_hit_id, moduleIds.size(), hitData.m_size);
  CUDA_CHECK(cudaGetLastError());
  CUDA_CHECK(cudaDeviceSynchronize());

  auto M1 =
      torch::from_blob(edgeData.cuda_graph_M1_hits, edgeData.size,
                       at::TensorOptions().device(at::kCUDA).dtype(at::kInt));
  auto M2 =
      torch::from_blob(edgeData.cuda_graph_M2_hits, edgeData.size,
                       at::TensorOptions().device(at::kCUDA).dtype(at::kInt));
  auto edgeIndex = torch::stack({M1, M2}, 1).transpose(1, 0).to(torch::kLong);
  ACTS_VERBOSE("edge index reshaped:\n"
               << edgeIndex.index({Slice(), Slice(0, 10)}));

  // Make node features
  float *cudaNodeFeatures{};
  CUDA_CHECK(cudaMalloc(&cudaNodeFeatures, features.size() * sizeof(float)));
  CUDA_CHECK(cudaMemcpy(cudaNodeFeatures, features.data(),
                        features.size() * sizeof(float),
                        cudaMemcpyHostToDevice));

  auto nodeFeatures =
      torch::from_blob(cudaNodeFeatures, inputValues.size(),
                       //[](void *ptr) { CUDA_CHECK(cudaFree(ptr)); },
                       at::TensorOptions().device(at::kCUDA).dtype(at::kFloat));
  nodeFeatures = nodeFeatures.reshape({(long)nHits, (long)nFeatures});

  // Make edge features
  float *edgeFeaturePtr{};
  CUDA_CHECK(cudaMalloc(&edgeFeaturePtr, 6 * edgeData.size * sizeof(float)));

  Acts::detail::makeEdgeFeatures<<<gridDimEdges, blockDim>>>(
      edgeData.size, edgeData.cuda_graph_M1_hits, edgeData.cuda_graph_M2_hits,
      nFeatures, cudaNodeFeatures, edgeFeaturePtr);
  CUDA_CHECK(cudaDeviceSynchronize());
  CUDA_CHECK(cudaGetLastError());

  auto edgeFeatures =
      torch::from_blob(edgeFeaturePtr, 6 * static_cast<long>(edgeData.size),
                       //[](void *ptr) { CUDA_CHECK(cudaFree(ptr)); },
                       at::TensorOptions().device(at::kCUDA).dtype(at::kFloat));
  edgeFeatures = edgeFeatures.reshape({static_cast<long>(edgeData.size), 6});

  ACTS_VERBOSE("edge featuers reshaped:\n"
               << edgeFeatures.index({Slice(None, 5), Slice()}));
  auto edgeFeaturesNew = edgeFeatures.clone();
  // copyFromDeviceAndPrint(edgeData.cuda_graph_dR, edgeData.size,
  // "cuda_graph_dR");
  //std::cout << "edgeIndex:\n" << edgeIndex << std::endl;
  CUDA_CHECK(cudaDeviceSynchronize());

  return {nodeFeatures, edgeIndex, edgeFeatures};
}

std::pair<detail::CUDA_edge_data<float>, detail::CUDA_hit_data<float>>
ModuleMapCuda::makeEdges(detail::CUDA_hit_data<float> cuda_TThits,
                         int *cuda_hit_indice) const {
  const dim3 block_dim = m_cfg.gpuBlocks;
  // ----------------------------------
  // memory allocation for hits + edges
  // ----------------------------------
  int nb_doublets = m_cudaModuleMapDoublet->size();
  int *cuda_M1_hits, *cuda_M2_hits;
  CUDA_CHECK(cudaMalloc(&cuda_M1_hits,
                        nb_doublets * m_cfg.maxEdgesAllocate * sizeof(int)));
  CUDA_CHECK(cudaMalloc(&cuda_M2_hits,
                        nb_doublets * m_cfg.maxEdgesAllocate * sizeof(int)));

  int *cuda_nb_edges;
  CUDA_CHECK(cudaMalloc(&cuda_nb_edges, nb_doublets * sizeof(int)));

  // ---------------------------------------------
  // loop over module doublets to kill connections
  // ---------------------------------------------

  CUDA_PRINTV(cuda_TThits.cuda_R(), 10);
  CUDA_PRINTV(cuda_TThits.cuda_R() + cuda_TThits.size() - 10, 10);
  CUDA_PRINTV(cuda_TThits.cuda_phi(), 10);
  CUDA_PRINTV(cuda_TThits.cuda_phi() + cuda_TThits.size() - 10, 10);
  CUDA_PRINTV(cuda_TThits.cuda_eta(), 10);
  CUDA_PRINTV(cuda_TThits.cuda_eta() + cuda_TThits.size() - 10, 10);
  CUDA_PRINTV(cuda_TThits.cuda_z(), 10);
  CUDA_PRINTV(cuda_TThits.cuda_z() + cuda_TThits.size() - 10, 10);
  CUDA_PRINTV(cuda_hit_indice, 10);
  CUDA_PRINTV(
      cuda_hit_indice + m_cudaModuleMapDoublet->module_map().size() - 10, 10);

  dim3 grid_dim = ((nb_doublets + block_dim.x - 1) / block_dim.x);
  doublet_cuts<float><<<grid_dim, block_dim>>>(
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
      std::numeric_limits<float>::max(), cuda_M1_hits, cuda_M2_hits,
      cuda_nb_edges, m_cfg.maxEdgesAllocate);
  CUDA_CHECK(cudaGetLastError());
  cudaDeviceSynchronize();

  CUDA_PRINTV(cuda_nb_edges, 10)
  CUDA_PRINTV(cuda_M1_hits, 10)
  //----------------
  // sum nb of edges
  //----------------

  int *cuda_edge_sum;
  CUDA_CHECK(cudaMalloc(&cuda_edge_sum, (nb_doublets + 1) * sizeof(int)));
  scan(cuda_edge_sum, cuda_nb_edges, m_cfg.gpuBlocks, nb_doublets);

  //---------------------------------------------
  // reduce nb of hits and sort by hid id M2 hits
  //---------------------------------------------

  int nb_doublet_edges;
  CUDA_CHECK(cudaMemcpy(&nb_doublet_edges, &(cuda_edge_sum[nb_doublets]),
                        sizeof(int), cudaMemcpyDeviceToHost));

  int *cuda_reduced_M1_hits, *cuda_reduced_M2_hits;
  CUDA_CHECK(cudaMalloc(&cuda_reduced_M1_hits, nb_doublet_edges * sizeof(int)));
  CUDA_CHECK(cudaMalloc(&cuda_reduced_M2_hits, nb_doublet_edges * sizeof(int)));
  compact_stream<<<grid_dim, block_dim>>>(cuda_reduced_M1_hits, cuda_M1_hits,
                                          cuda_edge_sum, m_cfg.maxEdgesAllocate,
                                          nb_doublets);
  CUDA_CHECK(cudaGetLastError());
  compact_stream<<<grid_dim, block_dim>>>(cuda_reduced_M2_hits, cuda_M2_hits,
                                          cuda_edge_sum, m_cfg.maxEdgesAllocate,
                                          nb_doublets);
  CUDA_CHECK(cudaGetLastError());

  int *cuda_sorted_M2_hits;
  CUDA_CHECK(cudaMalloc(&cuda_sorted_M2_hits, nb_doublet_edges * sizeof(int)));

  grid_dim = ((nb_doublet_edges + block_dim.x - 1) / block_dim.x);
  init_vector<<<grid_dim, block_dim>>>(cuda_sorted_M2_hits, nb_doublet_edges);
  CUDA_CHECK(cudaGetLastError());

  grid_dim = ((nb_doublets + block_dim.x - 1) / block_dim.x);
  partial_quick_sort<<<grid_dim, block_dim>>>(
      cuda_sorted_M2_hits, cuda_reduced_M2_hits, cuda_edge_sum, nb_doublets);
  CUDA_CHECK(cudaGetLastError());

  // -----------------------------
  // build doublets geometric cuts
  // -----------------------------

  float *cuda_z0, *cuda_phi_slope, *cuda_deta, *cuda_dphi;
  CUDA_CHECK(cudaMalloc(&cuda_z0, nb_doublet_edges * sizeof(float)));
  CUDA_CHECK(cudaMalloc(&cuda_phi_slope, nb_doublet_edges * sizeof(float)));
  CUDA_CHECK(cudaMalloc(&cuda_deta, nb_doublet_edges * sizeof(float)));
  CUDA_CHECK(cudaMalloc(&cuda_dphi, nb_doublet_edges * sizeof(float)));
  grid_dim = ((nb_doublet_edges + block_dim.x - 1) / block_dim.x);
  hits_geometric_cuts<<<grid_dim, block_dim>>>(
      cuda_z0, cuda_phi_slope, cuda_deta, cuda_dphi, cuda_reduced_M1_hits,
      cuda_reduced_M2_hits, cuda_TThits.cuda_R(), cuda_TThits.cuda_z(),
      cuda_TThits.cuda_eta(), cuda_TThits.cuda_phi(), TMath::Pi(),
      std::numeric_limits<float>::max(), nb_doublet_edges);
  // free mem at the end

  CUDA_mask cuda_edge_tag(nb_doublet_edges);
  CUDA_CHECK(cudaMemset(cuda_edge_tag.cuda_mask(), 0,
                        nb_doublet_edges * sizeof(bool)));

  cudaDeviceSynchronize();

  // -------------------------
  // loop over module triplets
  // -------------------------
  CUDA_PRINTV(cuda_sorted_M2_hits, 10)
  int nb_triplets = m_cudaModuleMapTriplet->size();
  grid_dim = ((nb_triplets + block_dim.x - 1) / block_dim.x);

  int *cuda_vertices;
  CUDA_CHECK(cudaMalloc(&cuda_vertices, cuda_TThits.size() * sizeof(int)));
  CUDA_CHECK(cudaMemset(cuda_vertices, 0, cuda_TThits.size() * sizeof(int)));

  triplet_cuts<float><<<grid_dim, block_dim>>>(
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
      std::numeric_limits<float>::max(), cuda_reduced_M1_hits,
      cuda_reduced_M2_hits, cuda_sorted_M2_hits, cuda_edge_sum, cuda_vertices,
      cuda_edge_tag.cuda_mask());
  CUDA_CHECK(cudaGetLastError());
  CUDA_PRINTV(cuda_edge_tag.cuda_mask(), 20)
  cudaDeviceSynchronize();

  //----------------
  // edges reduction
  //----------------

  const int nb_graph_edges = cuda_edge_tag.scan(m_cfg.gpuBlocks);

  int *cuda_graph_M1_hits, *cuda_graph_M2_hits;
  float *cuda_graph_dR, *cuda_graph_dz, *cuda_graph_deta, *cuda_graph_dphi,
      *cuda_graph_phi_slope, *cuda_graph_r_phi_slope;
  cudaMalloc(&cuda_graph_M1_hits, nb_graph_edges * sizeof(int));
  cudaMalloc(&cuda_graph_M2_hits, nb_graph_edges * sizeof(int));
  cudaMalloc(&cuda_graph_deta, nb_graph_edges * sizeof(float));
  cudaMalloc(&cuda_graph_dphi, nb_graph_edges * sizeof(float));
  cudaMalloc(&cuda_graph_phi_slope, nb_graph_edges * sizeof(float));

  grid_dim = ((nb_doublet_edges + block_dim.x - 1) / block_dim.x);
  compact_stream<<<grid_dim, block_dim>>>(
      cuda_graph_M1_hits, cuda_reduced_M1_hits, cuda_edge_tag.cuda_mask(),
      cuda_edge_tag.cuda_mask_sum(), nb_doublet_edges);
  CUDA_CHECK(cudaGetLastError());
  compact_stream<<<grid_dim, block_dim>>>(
      cuda_graph_M2_hits, cuda_reduced_M2_hits, cuda_edge_tag.cuda_mask(),
      cuda_edge_tag.cuda_mask_sum(), nb_doublet_edges);
  CUDA_CHECK(cudaGetLastError());
  compact_stream<<<grid_dim, block_dim>>>(
      cuda_graph_deta, cuda_deta, cuda_edge_tag.cuda_mask(),
      cuda_edge_tag.cuda_mask_sum(), nb_doublet_edges);
  CUDA_CHECK(cudaGetLastError());
  compact_stream<<<grid_dim, block_dim>>>(
      cuda_graph_dphi, cuda_dphi, cuda_edge_tag.cuda_mask(),
      cuda_edge_tag.cuda_mask_sum(), nb_doublet_edges);
  CUDA_CHECK(cudaGetLastError());
  compact_stream<<<grid_dim, block_dim>>>(
      cuda_graph_phi_slope, cuda_phi_slope, cuda_edge_tag.cuda_mask(),
      cuda_edge_tag.cuda_mask_sum(), nb_doublet_edges);
  CUDA_CHECK(cudaGetLastError());

  ACTS_VERBOSE("nb_graph_edges: " << nb_graph_edges);
  CUDA_PRINTV(cuda_graph_M1_hits, 10);
  CUDA_PRINTV(cuda_graph_M2_hits, 10);
  CUDA_PRINTV(cuda_graph_M1_hits + nb_graph_edges - 10, 10);
  CUDA_PRINTV(cuda_graph_M2_hits + nb_graph_edges - 10, 10);

  //----------------
  // nodes reduction
  //----------------

  int nb_hits = cuda_TThits.size();
  int *cuda_graph_vertices_sum;
  cudaMalloc(&cuda_graph_vertices_sum, (nb_hits + 1) * sizeof(int));
  scan(cuda_graph_vertices_sum, cuda_vertices, m_cfg.gpuBlocks, nb_hits);
  int nb_graph_hits;
  cudaMemcpy(&nb_graph_hits, &(cuda_graph_vertices_sum[nb_hits]), sizeof(int),
             cudaMemcpyDeviceToHost);

  std::uint64_t *cuda_graph_hitsID;
  float *cuda_graph_R, *cuda_graph_z, *cuda_graph_eta, *cuda_graph_phi;
  cudaMalloc(&cuda_graph_hitsID, nb_graph_hits * sizeof(std::uint64_t));
  cudaMalloc(&cuda_graph_R, nb_graph_hits * sizeof(float));
  cudaMalloc(&cuda_graph_z, nb_graph_hits * sizeof(float));
  cudaMalloc(&cuda_graph_eta, nb_graph_hits * sizeof(float));
  cudaMalloc(&cuda_graph_phi, nb_graph_hits * sizeof(float));

  grid_dim = ((nb_hits + block_dim.x - 1) / block_dim.x);
  compact_stream<<<grid_dim, block_dim>>>(
      cuda_graph_hitsID, cuda_TThits.cuda_hit_id(), cuda_vertices,
      cuda_graph_vertices_sum, nb_hits);
  CUDA_CHECK(cudaGetLastError());
  compact_stream<<<grid_dim, block_dim>>>(cuda_graph_R, cuda_TThits.cuda_R(),
                                          cuda_vertices,
                                          cuda_graph_vertices_sum, nb_hits);
  CUDA_CHECK(cudaGetLastError());
  compact_stream<<<grid_dim, block_dim>>>(cuda_graph_z, cuda_TThits.cuda_z(),
                                          cuda_vertices,
                                          cuda_graph_vertices_sum, nb_hits);
  CUDA_CHECK(cudaGetLastError());
  compact_stream<<<grid_dim, block_dim>>>(cuda_graph_phi,
                                          cuda_TThits.cuda_phi(), cuda_vertices,
                                          cuda_graph_vertices_sum, nb_hits);
  CUDA_CHECK(cudaGetLastError());

  CUDA_PRINTV(cuda_graph_hitsID, 20);

  cudaMalloc(&cuda_graph_dR, nb_graph_edges * sizeof(float));
  cudaMalloc(&cuda_graph_dz, nb_graph_edges * sizeof(float));
  cudaMalloc(&cuda_graph_r_phi_slope, nb_graph_edges * sizeof(float));
  grid_dim = ((nb_graph_edges + block_dim.x - 1) / block_dim.x);
  edge_features<<<grid_dim, block_dim>>>(
      cuda_graph_dR, cuda_graph_dz, cuda_graph_r_phi_slope, cuda_graph_M1_hits,
      cuda_graph_M2_hits, cuda_graph_vertices_sum, cuda_TThits.cuda_R(),
      cuda_TThits.cuda_z(), cuda_graph_phi_slope, TMath::Pi(),
      std::numeric_limits<float>::max(), nb_graph_edges);
  CUDA_CHECK(cudaGetLastError());

  cudaDeviceSynchronize();

  cudaFree(cuda_dphi);
  cudaFree(cuda_deta);
  cudaFree(cuda_phi_slope);
  cudaFree(cuda_z0);
  cudaFree(cuda_sorted_M2_hits);
  cudaFree(cuda_reduced_M2_hits);
  cudaFree(cuda_reduced_M1_hits);
  cudaFree(cuda_edge_sum);
  cudaFree(cuda_nb_edges);
  cudaFree(cuda_M2_hits);
  cudaFree(cuda_M1_hits);
  cudaFree(cuda_graph_vertices_sum);

  detail::CUDA_edge_data<float> edge_data;
  edge_data.size = nb_graph_edges;
  edge_data.cuda_graph_M1_hits = cuda_graph_M1_hits;
  edge_data.cuda_graph_M2_hits = cuda_graph_M2_hits;
  edge_data.cuda_graph_dR = cuda_graph_dR;
  edge_data.cuda_graph_dz = cuda_graph_dz;
  edge_data.cuda_graph_deta = cuda_graph_deta;
  edge_data.cuda_graph_dphi = cuda_graph_dphi;
  edge_data.cuda_graph_r_phi_slope = cuda_graph_r_phi_slope;
  edge_data.cuda_graph_phi_slope = cuda_graph_phi_slope;

  detail::CUDA_hit_data<float> hit_data;
  hit_data.m_size = nb_graph_hits;
  hit_data.m_cuda_R = cuda_graph_R;
  hit_data.m_cuda_z = cuda_graph_z;
  hit_data.m_cuda_eta = cuda_graph_eta;
  hit_data.m_cuda_phi = cuda_graph_phi;
  hit_data.m_cuda_hit_id = cuda_graph_hitsID;
  CUDA_PRINTV(edge_data.cuda_graph_M1_hits, 10);
  CUDA_PRINTV(edge_data.cuda_graph_M2_hits, 10);
  return {edge_data, hit_data};
}

}  // namespace Acts
