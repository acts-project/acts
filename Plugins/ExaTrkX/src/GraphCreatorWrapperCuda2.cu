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
#include <torch/csrc/autograd/generated/variable_factories.h>

namespace {

inline void cudaAssert(cudaError_t code, const char *file, int line) {
  if (code != cudaSuccess) {
    std::stringstream ss;
    ss << "CUDA error: " << cudaGetErrorString(code) << ", " << file << ":"
       << line;
    throw std::runtime_error(ss.str());
    // std::cout << ss.str() << std::endl;
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

__global__ void setHitId(std::size_t nbHits, int *hitIds) {
  std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i >= nbHits) {
    return;
  }

  hitIds[i] = i;
}

}  // namespace

namespace Acts::detail {

GraphCreatorWrapperCuda2::GraphCreatorWrapperCuda2(const std::string &path,
                                                   int device, int blocks) {
  m_graphCreator = std::make_unique<CUDA_graph_creator<float>>(
      blocks, device, path, 10,
      std::pair<float, float>{0.f, std::numeric_limits<float>::max()});
}

GraphCreatorWrapperCuda2::~GraphCreatorWrapperCuda2() {}

std::pair<at::Tensor, at::Tensor> GraphCreatorWrapperCuda2::build(
    const std::vector<float> &features,
    const std::vector<std::uint64_t> &moduleIds, const Acts::Logger &logger) {
  using GC = CUDA_graph_creator<float>;
  const auto nHits = moduleIds.size();
  const auto nFeatures = features.size() / moduleIds.size();

  GC::input_hit_data_t inputData;
  inputData.nb_graph_hits = moduleIds.size();

  std::size_t rOffset = 0;
  std::size_t phiOffset = 1;
  std::size_t zOffset = 2;
  std::size_t etaOffset = 3;

  const auto srcStride = sizeof(float) * nFeatures;
  const auto dstStride = sizeof(float);  // contiguous in destination
  const auto width = sizeof(float);      // only copy 1 column
  const auto height = nHits;

  CUDA_CHECK(cudaMallocT(&inputData.cuda_R, nHits * sizeof(float)));
  CUDA_CHECK(cudaMemcpy2D(inputData.cuda_R, dstStride,
                          features.data() + rOffset, srcStride, width, height,
                          cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMallocT(&inputData.cuda_phi, nHits * sizeof(float)));
  CUDA_CHECK(cudaMemcpy2D(inputData.cuda_phi, dstStride,
                          features.data() + phiOffset, srcStride, width, height,
                          cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMallocT(&inputData.cuda_z, nHits * sizeof(float)))
  CUDA_CHECK(cudaMemcpy2D(inputData.cuda_z, dstStride,
                          features.data() + zOffset, srcStride, width, height,
                          cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMallocT(&inputData.cuda_eta, nHits * sizeof(float)));
  CUDA_CHECK(cudaMemcpy2D(inputData.cuda_eta, dstStride,
                          features.data() + etaOffset, srcStride, width, height,
                          cudaMemcpyHostToDevice));

  CUDA_CHECK(cudaMallocT(&inputData.cuda_x, nHits * sizeof(float)));
  CUDA_CHECK(cudaMallocT(&inputData.cuda_y, nHits * sizeof(float)));

  dim3 blockDim = 512;
  dim3 gridDim = (nHits + blockDim.x - 1) / blockDim.x;
  computeXandY<<<gridDim, blockDim>>>(nHits, inputData.cuda_x, inputData.cuda_y,
                                      inputData.cuda_R, inputData.cuda_phi);
  CUDA_CHECK(cudaGetLastError());

  CUDA_CHECK(cudaMallocT(&inputData.cuda_hit_id, nHits * sizeof(std::size_t)));
  CUDA_CHECK(cudaMemcpy(inputData.cuda_hit_id, moduleIds.data(),
                        nHits * sizeof(std::size_t), cudaMemcpyHostToDevice));

  CUDA_CHECK(cudaMallocT(&inputData.cuda_hit_indices, nHits * sizeof(int)));
  setHitId<<<gridDim, blockDim>>>(nHits * sizeof(int),
                                  inputData.cuda_hit_indices);
  CUDA_CHECK(cudaGetLastError());

  GC::statistics_t stats;

  auto [edgeData, outHitData] =
      m_graphCreator->build_implementation(inputData, stats);

  int *edgeIndexPtr{};
  CUDA_CHECK(
      cudaMallocT(&edgeIndexPtr, 2 * edgeData.nb_graph_edges * sizeof(int)));
  CUDA_CHECK(cudaMemcpy(edgeIndexPtr, edgeData.cuda_graph_M1_hits,
                        edgeData.nb_graph_edges * sizeof(int),
                        cudaMemcpyDeviceToDevice));
  CUDA_CHECK(cudaMemcpy(
      edgeIndexPtr + edgeData.nb_graph_edges, edgeData.cuda_graph_M2_hits,
      edgeData.nb_graph_edges * sizeof(int), cudaMemcpyDeviceToDevice));

  auto edgeIndex =
      torch::from_blob(edgeIndexPtr,
                       {2, static_cast<long>(edgeData.nb_graph_edges)},
                       at::TensorOptions().device(at::kCUDA).dtype(at::kInt))
          .to(at::kLong);

  CUDA_CHECK(cudaFree(edgeIndexPtr));

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

  auto edgeFeatures = torch::from_blob(
      edgeFeaturePtr, {static_cast<long>(edgeData.nb_graph_edges), 6},
      [](void *ptr) { CUDA_CHECK(cudaFree(ptr)); },
      at::TensorOptions().device(at::kCUDA).dtype(at::kFloat));

  cudaFree(edgeData.cuda_graph_dR);
  cudaFree(edgeData.cuda_graph_dphi);
  cudaFree(edgeData.cuda_graph_dz);
  cudaFree(edgeData.cuda_graph_deta);
  cudaFree(edgeData.cuda_graph_phi_slope);
  cudaFree(edgeData.cuda_graph_r_phi_slope);
  cudaFree(edgeData.cuda_graph_M1_hits);
  cudaFree(edgeData.cuda_graph_M2_hits);

  cudaFree(outHitData.cuda_graph_R);
  cudaFree(outHitData.cuda_graph_phi);
  cudaFree(outHitData.cuda_graph_z);
  cudaFree(outHitData.cuda_graph_eta);
  cudaFree(outHitData.cuda_graph_phi_slope);
  cudaFree(outHitData.cuda_graph_r_phi_slope);

  return {edgeIndex, edgeFeatures};
}

}  // namespace Acts::detail
