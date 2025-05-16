// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/Tensor.hpp"
#include "Acts/Plugins/ExaTrkX/detail/CudaUtils.hpp"

#include <thrust/copy.h>
#include <thrust/count.h>
#include <thrust/execution_policy.h>

namespace {

__global__ void sigmoidImpl(std::size_t size, float *array) {
  std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= size) {
    return;
  }

  array[i] = 1.f / (1.f + __expf(-array[i]));
}

__global__ void applyCut(std::size_t size, float cutoff, const float *array,
                         bool *mask) {
  std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= size) {
    return;
  }

  mask[i] = array[i] > cutoff;
}

}  // namespace

namespace Acts::detail {

void cudaSigmoid(Tensor<float> &tensor, const ExecutionContext &execContext) {
  dim3 blockDim = 1024;
  dim3 gridDim = (tensor.size() + blockDim.x - 1) / blockDim.x;
  sigmoidImpl<<<blockDim, gridDim, 0, *execContext.stream>>>(tensor.size(),
                                                             tensor.data());
  ACTS_CUDA_CHECK(cudaGetLastError());
}

std::pair<Tensor<float>, Tensor<std::int64_t>> cudaApplyScoreCut(
    const Tensor<float> &scores, const Tensor<std::int64_t> &edgeIndex,
    float cut, const ExecutionContext &execContext) {
  dim3 blockDim = 1024;
  dim3 gridDim = (scores.size() + blockDim.x - 1) / blockDim.x;

  bool *mask{};
  ACTS_CUDA_CHECK(cudaMallocAsync(&mask, scores.size() * sizeof(bool),
                                  *execContext.stream));

  applyCut<<<blockDim, gridDim, 0, *execContext.stream>>>(scores.size(), cut,
                                                          scores.data(), mask);
  ACTS_CUDA_CHECK(cudaGetLastError());

  const std::size_t nEdgesAfter = thrust::count(
      thrust::device.on(*execContext.stream), mask, mask + scores.size(), true);

  auto outputScores = Tensor<float>::Create({nEdgesAfter, 1}, execContext);
  auto outputEdgeIndex =
      Tensor<std::int64_t>::Create({2, nEdgesAfter}, execContext);

  auto pred = [] __device__(bool x) { return x; };
  thrust::copy_if(thrust::device.on(*execContext.stream), scores.data(),
                  scores.data() + scores.size(), mask, outputScores.data(),
                  pred);

  const auto edgesBefore = edgeIndex.size() / 2;
  thrust::copy_if(thrust::device.on(*execContext.stream), edgeIndex.data(),
                  edgeIndex.data() + edgesBefore, mask, outputEdgeIndex.data(),
                  pred);
  thrust::copy_if(thrust::device.on(*execContext.stream),
                  edgeIndex.data() + edgesBefore,
                  edgeIndex.data() + 2 * edgesBefore, mask,
                  outputEdgeIndex.data() + nEdgesAfter, pred);

  ACTS_CUDA_CHECK(cudaFreeAsync(mask, *execContext.stream));
  return {std::move(outputScores), std::move(outputEdgeIndex)};
}

}  // namespace Acts::detail
