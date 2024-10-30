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

#include <CUDA_graph_creator>
#include <TTree_hits>
#include <graph>

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <torch/torch.h>

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
                             const T *cuda_R, const T *cuda_phi,
                             double phi_scale) {
  std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i >= nbHits) {
    return;
  }

  double r = cuda_R[i];
  double phi = cuda_phi[i] * phi_scale;

  cuda_x[i] = r * std::cos(phi);
  cuda_y[i] = r * std::sin(phi);

  printf("%lu: %f, %f -> %f, %f\n", i, cuda_R[i], cuda_phi[i], cuda_x[i],
         cuda_y[i]);
}

__global__ void setHitId(std::size_t nbHits, std::size_t *hitIds) {
  std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i >= nbHits) {
    return;
  }

  hitIds[i] = i;
}

}  // namespace

namespace Acts::detail {

GraphCreatorWrapperCuda::GraphCreatorWrapperCuda(const std::string &path,
                                                 int device, int blocks) {
  m_graphCreator = std::make_unique<CUDA_graph_creator<float>>(
      blocks, device, path, 10,
      std::pair<float, float>{0.f, std::numeric_limits<float>::max()});
}

GraphCreatorWrapperCuda::~GraphCreatorWrapperCuda() {}

std::pair<at::Tensor, at::Tensor> GraphCreatorWrapperCuda::build(
    const std::vector<float> &features,
    const std::vector<std::uint64_t> &moduleIds, const Acts::Logger &logger) {
  using GC = CUDA_graph_creator<float>;
  const auto nHits = moduleIds.size();
  const auto nFeatures = features.size() / moduleIds.size();

  dim3 blockDim = 512;
  dim3 gridDim = (nHits + blockDim.x - 1) / blockDim.x;

  dim3 block_dim = blockDim;
  dim3 grid_dim = gridDim;
#if 0
  // TODO understand this algorithm
  std::vector<int> hit_indice;
  {
    const auto &module_map =
        m_graphCreator->get_module_map_doublet().module_map();
    std::vector<int> nb_hits(module_map.size(), 0);

    for (std::size_t i = 0; i < nHits; ++i) {
      const auto it = module_map.find(moduleIds[i]);
      if (it != module_map.end()) {
        nb_hits[it->second] += 1;
      }
    }

    hit_indice.push_back(0);
    for (std::size_t i = 0; i < nb_hits.size(); i++) {
      hit_indice.push_back(hit_indice[i] + nb_hits[i]);
    }
  }

  GC::input_hit_data_t inputData;
  inputData.nb_graph_hits = moduleIds.size();

  std::size_t rOffset = 0;
  std::size_t phiOffset = 1;
  std::size_t zOffset = 2;
  std::size_t etaOffset = 3;

  //const auto srcStride = sizeof(float) * nFeatures;
  //const auto dstStride = sizeof(float);  // contiguous in destination
  //const auto width = sizeof(float);      // only copy 1 column
  //const auto height = nHits;

  std::vector<float> hostR(nHits), hostPhi(nHits), hostZ(nHits), hostEta(nHits);
  for(auto i=0ul; i<nHits; ++i) {
    hostR.at(i) = features.at(i*nFeatures + rOffset);
    hostPhi.at(i) = features.at(i*nFeatures + phiOffset);
    hostZ.at(i) = features.at(i*nFeatures + zOffset);
    hostEta.at(i) = features.at(i*nFeatures + etaOffset);
  }


  CUDA_CHECK(cudaMallocT(&inputData.cuda_R, nHits * sizeof(float)));
  CUDA_CHECK(cudaMemcpy(inputData.cuda_R, hostR.data(), nHits * sizeof(float), cudaMemcpyHostToDevice));
  //CUDA_CHECK(cudaMemcpy2D(inputData.cuda_R, dstStride,
  //                        features.data() + rOffset, srcStride, width, height,
  //                        cudaMemcpyHostToDevice));

  CUDA_CHECK(cudaMallocT(&inputData.cuda_phi, nHits * sizeof(float)));
  CUDA_CHECK(cudaMemcpy(inputData.cuda_phi, hostPhi.data(), nHits * sizeof(float), cudaMemcpyHostToDevice));
  //CUDA_CHECK(cudaMemcpy2D(inputData.cuda_phi, dstStride,
  //                        features.data() + phiOffset, srcStride, width, height,
  //                        cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMallocT(&inputData.cuda_z, nHits * sizeof(float)))
  CUDA_CHECK(cudaMemcpy(inputData.cuda_z, hostZ.data(), nHits * sizeof(float), cudaMemcpyHostToDevice));
  //CUDA_CHECK(cudaMemcpy2D(inputData.cuda_z, dstStride,
  //                        features.data() + zOffset, srcStride, width, height,
  //                        cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMallocT(&inputData.cuda_eta, nHits * sizeof(float)));
  CUDA_CHECK(cudaMemcpy(inputData.cuda_eta, hostEta.data(), nHits * sizeof(float), cudaMemcpyHostToDevice));
  //CUDA_CHECK(cudaMemcpy2D(inputData.cuda_eta, dstStride,
  //                        features.data() + etaOffset, srcStride, width, height,
  //                        cudaMemcpyHostToDevice));

  CUDA_CHECK(cudaMallocT(&inputData.cuda_x, nHits * sizeof(float)));
  CUDA_CHECK(cudaMallocT(&inputData.cuda_y, nHits * sizeof(float)));


  computeXandY<<<gridDim, blockDim>>>(nHits, inputData.cuda_x, inputData.cuda_y,
                                      inputData.cuda_R, inputData.cuda_phi, 3.14152654);
  CUDA_CHECK(cudaGetLastError());

  CUDA_CHECK(cudaMallocT(&inputData.cuda_hit_id, nHits * sizeof(std::size_t)));
  setHitId<<<gridDim, blockDim>>>(nHits, inputData.cuda_hit_id);
  CUDA_CHECK(cudaGetLastError());

  CUDA_CHECK(cudaMallocT(&inputData.cuda_hit_indices, nHits * sizeof(int)));
  CUDA_CHECK(cudaMemcpy(inputData.cuda_hit_indices, hit_indice.data(),
                        nHits * sizeof(int), cudaMemcpyHostToDevice));

  CUDA_CHECK(cudaGetLastError());
#else
  hits<float> hitsCollection(false, false);

  for (auto i = 0ul; i < nHits; ++i) {
    // TODO Use std::span when we move to C++20
    const float *hitFeatures = features.data() + i * nFeatures;

    int hitId = static_cast<int>(i);

    // Needs to be rescaled because ModuleMapGraph expects unscaled features
    float r = hitFeatures[0] * 1000.f;         // rScale;
    float phi = hitFeatures[1] * 3.141592654;  // phiScale;
    float z = hitFeatures[2] * 1000.f;         // zScale;

    float x = r * std::cos(phi);
    float y = r * std::sin(phi);

    std::uint64_t particleId = 0;  // We do not know
    std::uint64_t moduleId = moduleIds[i];
    std::string hardware = "";      // now hardware
    int barrelEndcap = 0;           // unclear, is this a flag???
    std::uint64_t particleID1 = 0;  // unclear
    std::uint64_t particleID2 = 0;  // unclear

    hit<float> hit(hitId, x, y, z, particleId, moduleId, hardware, barrelEndcap,
                   particleID1, particleID2);

    hitsCollection += hit;
  }

  TTree_hits<float> hitsTree = hitsCollection;

  CUDA_TTree_hits<float> input_hits;
  std::string event = "0";
  input_hits.add_event(event, hitsTree,
                       m_graphCreator->get_module_map_doublet().module_map());
  input_hits.HostToDevice();

  TTree_hits_constants<<<grid_dim, block_dim>>>(
      input_hits.size(), input_hits.cuda_x(), input_hits.cuda_y(),
      input_hits.cuda_z(), input_hits.cuda_R(), input_hits.cuda_eta(),
      input_hits.cuda_phi());
  cudaDeviceSynchronize();

  using input_hit_data_t = CUDA_graph_creator<float>::input_hit_data_t;
  input_hit_data_t input_hit_data;
  input_hit_data.nb_graph_hits = input_hits.size();
  input_hit_data.cuda_R = input_hits.cuda_R();
  input_hit_data.cuda_z = input_hits.cuda_z();
  input_hit_data.cuda_eta = input_hits.cuda_eta();
  input_hit_data.cuda_phi = input_hits.cuda_phi();
  input_hit_data.cuda_hit_id = input_hits.cuda_hit_id();
  input_hit_data.cuda_hit_indices = input_hits.cuda_hit_indice();
  input_hit_data.cuda_x = input_hits.cuda_x();
  input_hit_data.cuda_y = input_hits.cuda_y();

  auto &inputData = input_hit_data;
#endif

  GC::statistics_t stats;

  auto [edgeData, outHitData] =
      m_graphCreator->build_implementation(inputData, stats);
  CUDA_CHECK(cudaGetLastError());
  cudaDeviceSynchronize();

  std::cout << "Made " << stats.nb_doublet_edges << " doublet edges "
            << std::endl;
  std::cout << "Made " << edgeData.nb_graph_edges << " edges" << std::endl;
  assert(edgeData.nb_graph_edges > 0 && edgeData.nb_graph_edges < 100'000'000);

  int *edgeIndexPtr{};
  CUDA_CHECK(
      cudaMallocT(&edgeIndexPtr, 2 * edgeData.nb_graph_edges * sizeof(int)));
  CUDA_CHECK(cudaMemcpy(edgeIndexPtr, edgeData.cuda_graph_M1_hits,
                        edgeData.nb_graph_edges * sizeof(int),
                        cudaMemcpyDeviceToDevice));
  CUDA_CHECK(cudaMemcpy(
      edgeIndexPtr + edgeData.nb_graph_edges, edgeData.cuda_graph_M2_hits,
      edgeData.nb_graph_edges * sizeof(int), cudaMemcpyDeviceToDevice));
  CUDA_CHECK(cudaDeviceSynchronize());

  // std::cout << "Make edge index from blob" << std::endl;
  auto edgeIndex = torch::from_blob(
      edgeIndexPtr, 2 * static_cast<long>(edgeData.nb_graph_edges),
      at::TensorOptions().device(at::kCUDA).dtype(at::kInt));

  edgeIndex = edgeIndex.reshape({2, static_cast<long>(edgeData.nb_graph_edges)})
                  .to(at::kLong);
  // CUDA_CHECK(cudaFree(edgeIndexPtr));

  float *edgeFeaturePtr{};
  cudaMallocT(&edgeFeaturePtr, 6 * edgeData.nb_graph_edges * sizeof(float));

  CUDA_CHECK(cudaMemcpy2D(edgeFeaturePtr, 6 * sizeof(float),
                          edgeData.cuda_graph_dR, sizeof(float), sizeof(float),
                          edgeData.nb_graph_edges, cudaMemcpyDeviceToDevice));
  CUDA_CHECK(cudaMemcpy2D(edgeFeaturePtr + 1, 6 * sizeof(float),
                          edgeData.cuda_graph_dphi, sizeof(float),
                          sizeof(float), edgeData.nb_graph_edges,
                          cudaMemcpyDeviceToDevice));
  CUDA_CHECK(cudaMemcpy2D(edgeFeaturePtr + 2, 6 * sizeof(float),
                          edgeData.cuda_graph_dz, sizeof(float), sizeof(float),
                          edgeData.nb_graph_edges, cudaMemcpyDeviceToDevice));
  CUDA_CHECK(cudaMemcpy2D(edgeFeaturePtr + 3, 6 * sizeof(float),
                          edgeData.cuda_graph_deta, sizeof(float),
                          sizeof(float), edgeData.nb_graph_edges,
                          cudaMemcpyDeviceToDevice));
  CUDA_CHECK(cudaMemcpy2D(edgeFeaturePtr + 4, 6 * sizeof(float),
                          edgeData.cuda_graph_phi_slope, sizeof(float),
                          sizeof(float), edgeData.nb_graph_edges,
                          cudaMemcpyDeviceToDevice));
  CUDA_CHECK(cudaMemcpy2D(edgeFeaturePtr + 5, 6 * sizeof(float),
                          edgeData.cuda_graph_r_phi_slope, sizeof(float),
                          sizeof(float), edgeData.nb_graph_edges,
                          cudaMemcpyDeviceToDevice));
  CUDA_CHECK(cudaDeviceSynchronize());

  // std::cout << "Make edge features from blob" << std::endl;
  auto edgeFeatures = torch::from_blob(
      edgeFeaturePtr, {static_cast<long>(edgeData.nb_graph_edges), 6},
      //[](void *ptr) { CUDA_CHECK(cudaFree(ptr)); },
      at::TensorOptions().device(at::kCUDA).dtype(at::kFloat));

#if 0
  CUDA_CHECK(cudaFree(edgeData.cuda_graph_dR));
  CUDA_CHECK(cudaFree(edgeData.cuda_graph_dphi));
  CUDA_CHECK(cudaFree(edgeData.cuda_graph_dz));
  CUDA_CHECK(cudaFree(edgeData.cuda_graph_deta));
  CUDA_CHECK(cudaFree(edgeData.cuda_graph_phi_slope));
  CUDA_CHECK(cudaFree(edgeData.cuda_graph_r_phi_slope));
  CUDA_CHECK(cudaFree(edgeData.cuda_graph_M1_hits));
  CUDA_CHECK(cudaFree(edgeData.cuda_graph_M2_hits));
#endif

  CUDA_CHECK(cudaFree(outHitData.cuda_graph_R));
  CUDA_CHECK(cudaFree(outHitData.cuda_graph_phi));
  CUDA_CHECK(cudaFree(outHitData.cuda_graph_z));
  CUDA_CHECK(cudaFree(outHitData.cuda_graph_eta));
  CUDA_CHECK(cudaFree(outHitData.cuda_graph_phi_slope));
  CUDA_CHECK(cudaFree(outHitData.cuda_graph_r_phi_slope));

  auto edgeFeaturesNew = edgeFeatures.clone();

  // std::cout << edgeIndex << std::endl;
  // std::cout << edgeFeaturesNew << std::endl;
  CUDA_CHECK(cudaDeviceSynchronize());

  return {edgeIndex.clone(), edgeFeaturesNew};
}

}  // namespace Acts::detail
