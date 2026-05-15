// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Gnn/Tensor.hpp"
#include "ActsPlugins/Gnn/detail/CudaUtils.hpp"

#include <thrust/copy.h>
#include <thrust/count.h>
#include <thrust/execution_policy.h>
#include <thrust/iterator/counting_iterator.h>

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

/// Gather rows from src into dst using a device-side index list.
/// One thread per output row; cols are iterated in a short inner loop.
/// Assumes nCols is small (e.g. node/edge feature count) — not suitable for
/// large nCols since the inner loop becomes the bottleneck.
template <typename T>
__global__ void gatherRowsKernel(const std::size_t *indices, std::size_t nOut,
                                 std::size_t nCols, const T *src, T *dst) {
  const std::size_t k = blockIdx.x * blockDim.x + threadIdx.x;
  if (k >= nOut) {
    return;
  }
  const std::size_t srcRow = indices[k];
  for (std::size_t j = 0; j < nCols; ++j) {
    dst[k * nCols + j] = src[srcRow * nCols + j];
  }
}

/// Gather columns from src into dst using a device-side column index list.
/// One thread per output column; rows are iterated in a short inner loop.
/// Assumes nRows is small (e.g. 2 for an edge index tensor) — not suitable for
/// large nRows since the inner loop becomes the bottleneck.
template <typename T>
__global__ void gatherColsKernel(const std::size_t *indices, std::size_t nRows,
                                 std::size_t nColsSrc, std::size_t nColsDst,
                                 const T *src, T *dst) {
  const std::size_t col = blockIdx.x * blockDim.x + threadIdx.x;
  if (col >= nColsDst) {
    return;
  }
  for (std::size_t row = 0; row < nRows; ++row) {
    dst[row * nColsDst + col] = src[row * nColsSrc + indices[col]];
  }
}

}  // namespace

namespace ActsPlugins::detail {

void cudaSigmoid(Tensor<float> &tensor, cudaStream_t stream) {
  dim3 blockDim = 1024;
  dim3 gridDim = (tensor.size() + blockDim.x - 1) / blockDim.x;
  sigmoidImpl<<<gridDim, blockDim, 0, stream>>>(tensor.size(), tensor.data());
  ACTS_CUDA_CHECK(cudaGetLastError());
}

Tensor<bool> cudaScoreMask(const Tensor<float> &scores, float cut,
                           cudaStream_t stream) {
  ExecutionContext execContext{scores.device(), stream};
  auto mask = Tensor<bool>::Create(scores.shape(), execContext);

  dim3 blockDim = 1024;
  dim3 gridDim = (scores.size() + blockDim.x - 1) / blockDim.x;
  applyCut<<<gridDim, blockDim, 0, stream>>>(scores.size(), cut, scores.data(),
                                             mask.data());
  ACTS_CUDA_CHECK(cudaGetLastError());
  return mask;
}

template <typename T>
Tensor<T> cudaSelectRows(const Tensor<T> &tensor, const Tensor<bool> &mask,
                         const ExecutionContext &execContext) {
  const auto nCols = tensor.shape()[1];
  const auto nRows = tensor.shape()[0];
  const auto stream = execContext.stream.value();

  auto pred = [] __device__(bool x) { return x; };
  const std::size_t nSelected = thrust::count(
      thrust::device.on(stream), mask.data(), mask.data() + nRows, true);

  auto result = Tensor<T>::Create({nSelected, nCols}, execContext);

  if (nSelected == 0) {
    return result;
  }

  if (nCols == 1) {
    // Fast path: direct element selection via thrust stencil
    thrust::copy_if(thrust::device.on(stream), tensor.data(),
                    tensor.data() + nRows, mask.data(), result.data(), pred);
  } else {
    // General path: collect surviving row indices, then gather rows
    std::size_t *devRowIndices{};
    ACTS_CUDA_CHECK(cudaMallocAsync(&devRowIndices,
                                    nSelected * sizeof(std::size_t), stream));
    thrust::copy_if(thrust::device.on(stream),
                    thrust::make_counting_iterator<std::size_t>(0),
                    thrust::make_counting_iterator<std::size_t>(nRows),
                    mask.data(), devRowIndices, pred);

    dim3 blockDim = 1024;
    dim3 gridDim = (nSelected + blockDim.x - 1) / blockDim.x;
    gatherRowsKernel<<<gridDim, blockDim, 0, stream>>>(
        devRowIndices, nSelected, nCols, tensor.data(), result.data());
    ACTS_CUDA_CHECK(cudaGetLastError());

    ACTS_CUDA_CHECK(cudaFreeAsync(devRowIndices, stream));
  }

  return result;
}

template <typename T>
Tensor<T> cudaSelectCols(const Tensor<T> &tensor, const Tensor<bool> &mask,
                         const ExecutionContext &execContext) {
  const auto nRows = tensor.shape()[0];
  const auto nColsSrc = tensor.shape()[1];
  const auto stream = execContext.stream.value();

  auto pred = [] __device__(bool x) { return x; };
  const std::size_t nSelected = thrust::count(
      thrust::device.on(stream), mask.data(), mask.data() + nColsSrc, true);

  auto result = Tensor<T>::Create({nRows, nSelected}, execContext);

  if (nSelected == 0) {
    return result;
  }

  // Collect surviving column indices on device, then use gather kernel
  std::size_t *devColIndices{};
  ACTS_CUDA_CHECK(
      cudaMallocAsync(&devColIndices, nSelected * sizeof(std::size_t), stream));
  thrust::copy_if(thrust::device.on(stream),
                  thrust::make_counting_iterator<std::size_t>(0),
                  thrust::make_counting_iterator<std::size_t>(nColsSrc),
                  mask.data(), devColIndices, pred);

  dim3 blockDim = 1024;
  dim3 gridDim = (nSelected + blockDim.x - 1) / blockDim.x;
  gatherColsKernel<<<gridDim, blockDim, 0, stream>>>(
      devColIndices, nRows, nColsSrc, nSelected, tensor.data(), result.data());
  ACTS_CUDA_CHECK(cudaGetLastError());

  ACTS_CUDA_CHECK(cudaFreeAsync(devColIndices, stream));
  return result;
}

template Tensor<float> cudaSelectRows(const Tensor<float> &,
                                      const Tensor<bool> &,
                                      const ExecutionContext &);
template Tensor<std::int64_t> cudaSelectRows(const Tensor<std::int64_t> &,
                                             const Tensor<bool> &,
                                             const ExecutionContext &);
template Tensor<float> cudaSelectCols(const Tensor<float> &,
                                      const Tensor<bool> &,
                                      const ExecutionContext &);
template Tensor<std::int64_t> cudaSelectCols(const Tensor<std::int64_t> &,
                                             const Tensor<bool> &,
                                             const ExecutionContext &);

}  // namespace ActsPlugins::detail
