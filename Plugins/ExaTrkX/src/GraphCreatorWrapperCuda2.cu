// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/detail/GraphCreatorWrapper.hpp"

// #include <CUDA_graph_creator_new>
#include <CUDA_graph_creator>
#include <TTree_hits>
#include <algorithm>
#include <graph>
#include <ranges>

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <torch/torch.h>

#include "oldApiHelper.hpp"

namespace {

inline void cudaAssert(cudaError_t code, const char *file, int line) {
  if (code != cudaSuccess) {
    std::stringstream ss;
    ss << "CUDA error: " << cudaGetErrorString(code) << ", " << file << ":"
       << line;
    // throw std::runtime_error(ss.str());
    std::cout << ss.str() << std::endl;
  }
  cudaDeviceSynchronize();
}

#define CUDA_CHECK(ans) \
  { cudaAssert((ans), __FILE__, __LINE__); }

template <typename T>
cudaError_t cudaMallocT(T **ptr, std::size_t size) {
  return cudaMalloc((void **)ptr, size);
}

template <class T>
__global__ void computeXandY(std::size_t nbHits, T *cuda_x, T *cuda_y,
                             const T *cuda_R, const T *cuda_phi) {
  std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i >= nbHits) {
    return;
  }

  double r = cuda_R[i];
  double phi = cuda_phi[i];

  cuda_x[i] = r * std::cos(phi);
  cuda_y[i] = r * std::sin(phi);
}

template <class T>
__global__ void rescaleFeature(std::size_t nbHits, T *data, T scale) {
  std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i >= nbHits) {
    return;
  }

  data[i] *= scale;
}

constexpr float g_pi = 3.141592654f;

template <typename T>
__device__ T resetAngle(T angle) {
  if (angle > g_pi) {
    return angle - 2.f * g_pi;
  }
  if (angle < -g_pi) {
    return angle + 2.f * g_pi;
  }
  return angle;
};

template <typename T>
__global__ void makeEdgeFeatures(std::size_t nEdges, const int *srcEdges,
                                 const int *tgtEdges, std::size_t nNodeFeatures,
                                 const T *nodeFeatures, T *edgeFeatures) {
  enum NodeFeatures { r = 0, phi, z, eta };
  constexpr static int nEdgeFeatures = 6;

  std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i >= nEdges) {
    return;
  }

  const int src = srcEdges[i];
  const int tgt = tgtEdges[i];

  const T *srcNodeFeatures = nodeFeatures + src * nNodeFeatures;
  const T *tgtNodeFeatures = nodeFeatures + tgt * nNodeFeatures;

  T dr = tgtNodeFeatures[r] - srcNodeFeatures[r];
  T dphi =
      resetAngle(g_pi * (tgtNodeFeatures[phi] - srcNodeFeatures[phi])) / g_pi;
  T dz = tgtNodeFeatures[z] - srcNodeFeatures[z];
  T deta = tgtNodeFeatures[eta] - srcNodeFeatures[eta];
  T phislope = 0.0;
  T rphislope = 0.0;

  if (dr != 0.0) {
    phislope = std::clamp(dphi / dr, -100.f, 100.f);
    T avgR = T{0.5} * (tgtNodeFeatures[r] + srcNodeFeatures[r]);
    rphislope = avgR * phislope;
  }

  T *efPtr = edgeFeatures + i * nEdgeFeatures;
  efPtr[0] = dr;
  efPtr[1] = dphi;
  efPtr[2] = dz;
  efPtr[3] = deta;
  efPtr[4] = phislope;
  efPtr[5] = rphislope;
}

__global__ void setHitId(std::size_t nHits, std::uint64_t *hit_ids) {
  std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= nHits) {
    return;
  }
  hit_ids[i] = i;
}

__global__ void remapEdges(std::size_t nEdges, int *srcNodes, int *tgtNodes,
                           const std::uint64_t *hit_ids, std::size_t nAllNodes,
                           std::size_t nCompressedNodes) {
  std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= nEdges) {
    return;
  }

  srcNodes[i] = hit_ids[srcNodes[i]];
  tgtNodes[i] = hit_ids[tgtNodes[i]];
}

template <typename T>
void copyFromDeviceAndPrint(T *data, std::size_t size, std::string_view name) {
  std::vector<T> data_cpu(size);
  CUDA_CHECK(cudaMemcpy(data_cpu.data(), data, size * sizeof(T),
                        cudaMemcpyDeviceToHost));
  std::cout << name << "[" << size << "]: ";
  for (int i = 0; i < size; ++i) {
    std::cout << data_cpu.at(i) << "  ";
  }
  std::cout << std::endl;
}

}  // namespace

template <typename T>
void printPrefixSum(const std::vector<T> &vec) {
  T prev = vec.front();
  std::cout << 0 << ":" << vec[0] << "  ";
  for (auto i = 0ul; i < vec.size(); ++i) {
    if (vec[i] != prev) {
      std::cout << i << ":" << vec[i] << "  ";
      prev = vec[i];
    }
  }
}

namespace Acts::detail {

GraphCreatorWrapperCuda::GraphCreatorWrapperCuda(const std::string &path,
                                                 int device, int blocks) {
  m_graphCreator =
      // std::make_unique<CUDA_graph_creator<float>>(blocks, device, path);
      std::make_unique<CUDA_graph_creator<float>>(blocks, device, path, 0ul,
                                                  std::pair<float, float>{});
}

GraphCreatorWrapperCuda::~GraphCreatorWrapperCuda() {}

std::pair<at::Tensor, at::Tensor> GraphCreatorWrapperCuda::build(
    const std::vector<float> &features,
    const std::vector<std::uint64_t> &moduleIds, const Acts::Logger &logger) {
  auto lastError = cudaGetLastError();
  using GC = CUDA_graph_creator<float>;
  const auto nHits = moduleIds.size();
  const auto nFeatures = features.size() / moduleIds.size();

  const dim3 blockDim = 512;
  const dim3 gridDimHits = (nHits + blockDim.x - 1) / blockDim.x;

  assert(std::is_sorted(moduleIds.begin(), moduleIds.end()));

  // TODO understand this algorithm
  std::vector<int> hit_indice;
  {
    const auto &module_map =
        m_graphCreator->get_module_map_doublet().module_map();
    std::vector<int> nb_hits(module_map.size(), 0);
    std::vector<int> hits_bool_mask(nHits, 0);

    for (auto i = 0ul; i < nHits; ++i) {
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
  }

  assert(!std::all_of(hit_indice.begin(), hit_indice.end(),
                      [](auto v) { return v == 0; }));

  // std::cout << "my prefix sum (" << hit_indice.size() << "): ";
  // printPrefixSum(hit_indice);
  // std::cout << std::endl;

  float *cudaNodeFeatures{};
  CUDA_CHECK(cudaMalloc(&cudaNodeFeatures, features.size() * sizeof(float)));
  CUDA_CHECK(cudaMemcpy(cudaNodeFeatures, features.data(),
                        features.size() * sizeof(float),
                        cudaMemcpyHostToDevice));

#if 0
  CUDA_hit_data<float> inputData;
  inputData.m_size = nHits;

  std::size_t rOffset = 0;
  std::size_t phiOffset = 1;
  std::size_t zOffset = 2;
  std::size_t etaOffset = 3;

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

  CUDA_CHECK(cudaMallocT(&inputData.m_cuda_R, nHits * sizeof(float)));
  CUDA_CHECK(cudaMemcpy(inputData.m_cuda_R, hostR.data(), nHits * sizeof(float),
                        cudaMemcpyHostToDevice));
  // CUDA_CHECK(cudaMemcpy2D(inputData.cuda_R, dstStride,
  //                         features.data() + rOffset, srcStride, width,
  //                         height, cudaMemcpyHostToDevice));

  CUDA_CHECK(cudaMallocT(&inputData.m_cuda_phi, nHits * sizeof(float)));
  CUDA_CHECK(cudaMemcpy(inputData.m_cuda_phi, hostPhi.data(),
                        nHits * sizeof(float), cudaMemcpyHostToDevice));
  // CUDA_CHECK(cudaMemcpy2D(inputData.cuda_phi, dstStride,
  //                         features.data() + phiOffset, srcStride, width,
  //                         height, cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMallocT(&inputData.m_cuda_z, nHits * sizeof(float)))
  CUDA_CHECK(cudaMemcpy(inputData.m_cuda_z, hostZ.data(), nHits * sizeof(float),
                        cudaMemcpyHostToDevice));
  // CUDA_CHECK(cudaMemcpy2D(inputData.cuda_z, dstStride,
  //                         features.data() + zOffset, srcStride, width,
  //                         height, cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMallocT(&inputData.m_cuda_eta, nHits * sizeof(float)));
  CUDA_CHECK(cudaMemcpy(inputData.m_cuda_eta, hostEta.data(),
                        nHits * sizeof(float), cudaMemcpyHostToDevice));
  // CUDA_CHECK(cudaMemcpy2D(inputData.cuda_eta, dstStride,
  //                         features.data() + etaOffset, srcStride, width,
  //                         height, cudaMemcpyHostToDevice));

  CUDA_CHECK(cudaDeviceSynchronize());
  std::cout << "gridDimHits: " << gridDimHits.x << ", blockDim: " << blockDim.x << std::endl;
  rescaleFeature<<<gridDimHits, blockDim>>>(nHits, inputData.m_cuda_z, 1000.f);
  CUDA_CHECK(cudaGetLastError());
  rescaleFeature<<<gridDimHits, blockDim>>>(nHits, inputData.m_cuda_R, 1000.f);
  CUDA_CHECK(cudaGetLastError());
  rescaleFeature<<<gridDimHits, blockDim>>>(nHits, inputData.m_cuda_phi,
                                        3.14159f);
  CUDA_CHECK(cudaGetLastError());

  CUDA_CHECK(cudaMallocT(&inputData.m_cuda_x, nHits * sizeof(float)));
  CUDA_CHECK(cudaMallocT(&inputData.m_cuda_y, nHits * sizeof(float)));

  computeXandY<<<gridDimHits, blockDim>>>(nHits, inputData.m_cuda_x, inputData.m_cuda_y,
                                      inputData.m_cuda_R, inputData.m_cuda_phi);
  CUDA_CHECK(cudaGetLastError());

  CUDA_CHECK(cudaMallocT(&inputData.m_cuda_hit_id, nHits*sizeof(std::uint64_t)));
  setHitId<<<gridDimHits, blockDim>>>(nHits, inputData.m_cuda_hit_id);
  CUDA_CHECK(cudaGetLastError());

  int *cuda_hit_indice = nullptr;
  CUDA_CHECK(cudaMallocT(&cuda_hit_indice, hit_indice.size() * sizeof(int)));
  CUDA_CHECK(cudaMemcpy(cuda_hit_indice, hit_indice.data(),
                        hit_indice.size() * sizeof(int),
                        cudaMemcpyHostToDevice));

#else
  // Ref hit_indice
  auto hitsTree =
      makeTTreeHits(features, moduleIds, 1000.f, 3.141592654f, 1000.f);

  CUDA_TTree_hits<float> input_hits(
      hitsTree, m_graphCreator->get_module_map_doublet().module_map());

  // std::string event = "0";
  // input_hits.add_event(event, hitsTree,
  //                      m_graphCreator->get_module_map_doublet().module_map());
  input_hits.HostToDevice();
  cudaDeviceSynchronize();

  dim3 grid_dim = ((input_hits.size() + blockDim.x - 1) / blockDim.x);
  TTree_hits_constants<<<grid_dim, blockDim>>>(
      input_hits.size(), input_hits.cuda_x(), input_hits.cuda_y(),
      input_hits.cuda_z(), input_hits.cuda_R(), input_hits.cuda_eta(),
      input_hits.cuda_phi());
  CUDA_CHECK(cudaGetLastError());
  cudaDeviceSynchronize();

  CUDA_hit_data<float> input_hit_data;
  input_hit_data.m_size = input_hits.size();
  input_hit_data.m_cuda_R = input_hits.cuda_R();
  input_hit_data.m_cuda_z = input_hits.cuda_z();
  input_hit_data.m_cuda_eta = input_hits.cuda_eta();
  input_hit_data.m_cuda_phi = input_hits.cuda_phi();
  input_hit_data.m_cuda_x = input_hits.cuda_x();
  input_hit_data.m_cuda_y = input_hits.cuda_y();
  input_hit_data.m_cuda_hit_id = input_hits.cuda_hit_id();

  auto &inputData = input_hit_data;
  int *cuda_hit_indice = input_hits.cuda_hit_indice();
#endif
  /*
    std::vector<int> ref_hit_indice(
        m_graphCreator->get_module_map_doublet().module_map().size() + 1);
    CUDA_CHECK(
        cudaMemcpy(ref_hit_indice.data(), cuda_hit_indice,
                   (m_graphCreator->get_module_map_doublet().module_map().size()
    + 1) * sizeof(int), cudaMemcpyDeviceToHost));

  if( hit_indice == ref_hit_indice ) {
    std::cout << "hit_indice equal on cpu!!!" << std::endl;
  } else {
  assert(hit_indice.size() == ref_hit_indice.size());
  for(auto i=0ul; i< hit_indice.size(); ++i) {
    if(hit_indice.at(i) != ref_hit_indice.at(i)) {
        std::cout << "i: " << hit_indice.at(i) << " != ref " <<
  ref_hit_indice.at(i) << std::endl;
    }
  }
  }*/

  cudaDeviceSynchronize();
  CUDA_graph_creator<float>::graph_building_stats stats;
  const auto [edgeData, hitData] =
      m_graphCreator->build_impl2(inputData, cuda_hit_indice, stats, true);
  CUDA_CHECK(cudaGetLastError());
  CUDA_CHECK(cudaDeviceSynchronize());
  copyFromDeviceAndPrint(edgeData.cuda_graph_M1_hits, 10,
                         "M1 directly after build");
  copyFromDeviceAndPrint(edgeData.cuda_graph_M2_hits, 10,
                         "M2 directly after build");
  std::cout << "Made " << edgeData.size << " edges" << std::endl;
  std::cout << "number of hits in output data: " << hitData.m_size
            << ", before: " << nHits << std::endl;
  assert(edgeData.size > 0 && edgeData.size < 100'000'000);

  dim3 gridDimEdges = (edgeData.size + blockDim.x - 1) / blockDim.x;

  remapEdges<<<gridDimEdges, blockDim>>>(
      edgeData.size, edgeData.cuda_graph_M1_hits, edgeData.cuda_graph_M2_hits,
      hitData.m_cuda_hit_id, moduleIds.size(), hitData.m_size);
  CUDA_CHECK(cudaGetLastError());
  CUDA_CHECK(cudaDeviceSynchronize());

  {
    std::vector<int> edges_M1(edgeData.size), edges_M2(edgeData.size);
    std::vector<std::uint64_t> new_hit_indices(edgeData.size);
    CUDA_CHECK(cudaMemcpy(edges_M1.data(), edgeData.cuda_graph_M1_hits,
                          edgeData.size * sizeof(int), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(edges_M2.data(), edgeData.cuda_graph_M2_hits,
                          edgeData.size * sizeof(int), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(new_hit_indices.data(), hitData.m_cuda_hit_id,
                          hitData.m_size * sizeof(std::uint64_t),
                          cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaDeviceSynchronize());

    /*{
        std::ofstream of("edges_new_api.csv");
        of << "src,tgt\n";
        for(auto i=0ul; i<edgeData.size; ++i) {
          of << edges_M1.at(i) << "," << edges_M2.at(i) << std::endl;
        }
    }*/
    std::cout << "M1 is sorted: "
              << std::is_sorted(edges_M1.begin(), edges_M2.end()) << std::endl;
    std::cout << "max edges M1: "
              << *std::max_element(edges_M1.begin(), edges_M1.end())
              << std::endl;
    std::cout << "max edges M2: "
              << *std::max_element(edges_M2.begin(), edges_M2.end())
              << std::endl;
    std::cout << "max element new_hit_indices: "
              << *std::max_element(new_hit_indices.begin(),
                                   new_hit_indices.end())
              << std::endl;

    std::vector<int> idxs(edgeData.size);
    std::iota(idxs.begin(), idxs.end(), 0);
    std::sort(idxs.begin(), idxs.end(),
              [&](auto a, auto b) { return edges_M1.at(a) < edges_M1.at(b); });
    std::vector<int> edges_M1_sorted(edgeData.size),
        edges_M2_sorted(edgeData.size);
    std::transform(idxs.begin(), idxs.end(), edges_M1_sorted.begin(),
                   [&](auto idx) { return edges_M1.at(idx); });
    std::transform(idxs.begin(), idxs.end(), edges_M2_sorted.begin(),
                   [&](auto idx) { return edges_M2.at(idx); });

    std::array edgeVecs = {&edges_M1_sorted, &edges_M2_sorted};
    for (int n = 0; n < 2; ++n) {
      std::cout << "M" << n + 1 << " CPU: ";
      std::array starts = {0ul, edgeVecs[n]->size() - 10ul};
      bool p = true;
      for (auto s : starts) {
        for (int i = s; i < s + 10ul; ++i) {
          std::cout << edgeVecs[n]->at(i) << "\t";
        }
        if (p) {
          std::cout << " ... ";
          p = false;
        }
      }
      std::cout << std::endl;
    }
  }

  using namespace torch::indexing;
#if 0
  int *edgeIndexPtr{};
  CUDA_CHECK(cudaMallocT(&edgeIndexPtr, 2 * edgeData.size * sizeof(int)));
  CUDA_CHECK(cudaMemcpy(edgeIndexPtr, edgeData.cuda_graph_M1_hits,
                        edgeData.size * sizeof(int), cudaMemcpyDeviceToDevice));
  CUDA_CHECK(cudaMemcpy(edgeIndexPtr + edgeData.size,
                        edgeData.cuda_graph_M2_hits,
                        edgeData.size * sizeof(int), cudaMemcpyDeviceToDevice));
  CUDA_CHECK(cudaDeviceSynchronize());

  std::cout << "Creator edge Tensor" << std::endl;
  auto edgeIndex =
      torch::from_blob(edgeIndexPtr, 2 * static_cast<long>(edgeData.size),
                       at::TensorOptions().device(at::kCUDA).dtype(at::kInt));

  edgeIndex =
      edgeIndex.reshape({2, static_cast<long>(edgeData.size)}).to(at::kLong);

#else
  auto M1 =
      torch::from_blob(edgeData.cuda_graph_M1_hits, edgeData.size,
                       at::TensorOptions().device(at::kCUDA).dtype(at::kInt));
  auto M2 =
      torch::from_blob(edgeData.cuda_graph_M2_hits, edgeData.size,
                       at::TensorOptions().device(at::kCUDA).dtype(at::kInt));
  //std::cout << "M1:\n" << M1.index({Slice(0,10)}) << std::endl;
  //std::cout << "M2:\n" << M2.index({Slice(0,10)}) << std::endl;
  auto edgeIndex = torch::stack({M1, M2}, 1).transpose(1, 0).to(torch::kLong);
#endif
  std::cout << "edge index reshaped:\n"
            << edgeIndex.index({Slice(), Slice(0, 10)}) << std::endl;

  float *edgeFeaturePtr{};
  CUDA_CHECK(cudaMallocT(&edgeFeaturePtr, 6 * edgeData.size * sizeof(float)));

  makeEdgeFeatures<<<gridDimEdges, blockDim>>>(
      edgeData.size, edgeData.cuda_graph_M1_hits, edgeData.cuda_graph_M2_hits,
      nFeatures, cudaNodeFeatures, edgeFeaturePtr);
  CUDA_CHECK(cudaDeviceSynchronize());
  CUDA_CHECK(cudaGetLastError());

  copyFromDeviceAndPrint(edgeFeaturePtr, 20, "start of edge feature ptr");
  /*
    CUDA_CHECK(cudaMemcpy2D(edgeFeaturePtr, 6 * sizeof(float),
                            edgeData.cuda_graph_dR, sizeof(float),
    sizeof(float), edgeData.size, cudaMemcpyDeviceToDevice));
    CUDA_CHECK(cudaMemcpy2D(
        edgeFeaturePtr + 1, 6 * sizeof(float), edgeData.cuda_graph_dphi,
        sizeof(float), sizeof(float), edgeData.size, cudaMemcpyDeviceToDevice));
    CUDA_CHECK(cudaMemcpy2D(edgeFeaturePtr + 2, 6 * sizeof(float),
                            edgeData.cuda_graph_dz, sizeof(float),
    sizeof(float), edgeData.size, cudaMemcpyDeviceToDevice));
    CUDA_CHECK(cudaMemcpy2D(
        edgeFeaturePtr + 3, 6 * sizeof(float), edgeData.cuda_graph_deta,
        sizeof(float), sizeof(float), edgeData.size, cudaMemcpyDeviceToDevice));
    CUDA_CHECK(cudaMemcpy2D(
        edgeFeaturePtr + 4, 6 * sizeof(float), edgeData.cuda_graph_phi_slope,
        sizeof(float), sizeof(float), edgeData.size, cudaMemcpyDeviceToDevice));
    CUDA_CHECK(cudaMemcpy2D(
        edgeFeaturePtr + 5, 6 * sizeof(float), edgeData.cuda_graph_r_phi_slope,
        sizeof(float), sizeof(float), edgeData.size, cudaMemcpyDeviceToDevice));
    CUDA_CHECK(cudaDeviceSynchronize());
  */

  std::cout << "Create edge feature tensor" << std::endl;
  auto edgeFeatures =
      torch::from_blob(edgeFeaturePtr, 6 * static_cast<long>(edgeData.size),
                       //[](void *ptr) { CUDA_CHECK(cudaFree(ptr)); },
                       at::TensorOptions().device(at::kCUDA).dtype(at::kFloat));
  std::cout << "Reshape edge feature tensor" << std::endl;
  edgeFeatures = edgeFeatures.reshape({static_cast<long>(edgeData.size), 6});

  std::cout << "edge featuers reshaped:\n"
            << edgeFeatures.index({Slice(None, 5), Slice()}) << std::endl;
  auto edgeFeaturesNew = edgeFeatures.clone();
  // copyFromDeviceAndPrint(edgeData.cuda_graph_dR, edgeData.size,
  // "cuda_graph_dR");
  //std::cout << "edgeIndex:\n" << edgeIndex << std::endl;
  CUDA_CHECK(cudaDeviceSynchronize());

  // std::cout << "dR (0->1): " << hostR.at(1) - hostR.at(0) << std::endl;
  // std::cout << "dR (0->2): " << hostR.at(2) - hostR.at(0) << std::endl;
#if 1
  {
    cudaGetLastError();
    std::cout << "Run old API for reference:\n";
    auto builder = [&](auto &hits, bool print) {
      CUDA_graph_creator<float>::graph_building_stats stats;
      return m_graphCreator->build_impl(hits, stats, print);
    };
    auto [oldApiEdges, oldApiEdgeFeatures] = oldApiBuild(
        features, moduleIds, logger, builder, 1000.f, 3.14159f, 1000.f);
    std::cout << "old API edge shape: " << oldApiEdges.size(0) << ", "
              << oldApiEdges.size(1) << std::endl;
    std::cout << "old API edges:\n"
              << oldApiEdges.index({Slice(), Slice(None, 10)}) << std::endl;
    // return {oldApiEdges, oldApiEdgeFeatures};
  }
#endif
  return {edgeIndex.clone(), edgeFeaturesNew};
}

}  // namespace Acts::detail
