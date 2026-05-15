// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Gnn/Tensor.hpp"

#ifdef ACTS_GNN_WITH_CUDA
#include "ActsPlugins/Gnn/detail/CudaUtils.hpp"
#endif

#include <cstring>
#include <span>

namespace ActsPlugins {

namespace detail {

TensorPtr createTensorMemory(std::size_t nbytes,
                             const ExecutionContext &execContext) {
  if (execContext.device.type == Device::Type::eCPU) {
    void *ptr = new std::byte[nbytes];
    if (ptr == nullptr) {
      throw std::bad_alloc{};
    }
    return TensorPtr(ptr,
                     [](void *p) { delete[] static_cast<std::byte *>(p); });
  } else {
#ifdef ACTS_GNN_WITH_CUDA
    assert(execContext.stream.has_value());
    auto stream = *execContext.stream;
    void *ptr{};
    ACTS_CUDA_CHECK(cudaMallocAsync(&ptr, nbytes, stream));
    return TensorPtr(
        ptr, [stream](void *p) { ACTS_CUDA_CHECK(cudaFreeAsync(p, stream)); });
#else
    throw std::runtime_error(
        "Cannot create CUDA tensor, library was not compiled with CUDA");
#endif
  }
}

TensorPtr cloneTensorMemory(const TensorPtr &ptr, std::size_t nbytes,
                            Device devFrom, const ExecutionContext &to) {
  auto clone = createTensorMemory(nbytes, to);
  if (devFrom.isCpu() && to.device.isCpu()) {
    std::memcpy(clone.get(), ptr.get(), nbytes);
  } else {
#ifdef ACTS_GNN_WITH_CUDA
    assert(to.stream.has_value());
    if (devFrom.isCuda() && to.device.isCuda()) {
      ACTS_CUDA_CHECK(cudaMemcpyAsync(clone.get(), ptr.get(), nbytes,
                                      cudaMemcpyDeviceToDevice, *to.stream));
    } else if (devFrom.isCpu() && to.device.isCuda()) {
      ACTS_CUDA_CHECK(cudaMemcpyAsync(clone.get(), ptr.get(), nbytes,
                                      cudaMemcpyHostToDevice, *to.stream));
    } else if (devFrom.isCuda() && to.device.isCpu()) {
      ACTS_CUDA_CHECK(cudaMemcpyAsync(clone.get(), ptr.get(), nbytes,
                                      cudaMemcpyDeviceToHost, *to.stream));
    }
#else
    throw std::runtime_error(
        "Cannot clone CUDA tensor, library was not compiled with CUDA");
#endif
  }
  return clone;
}

void cudaSigmoid(Tensor<float> &tensor, cudaStream_t stream);

Tensor<bool> cudaScoreMask(const Tensor<float> &scores, float cut,
                           cudaStream_t stream);

template <typename T>
Tensor<T> cudaSelectRows(const Tensor<T> &tensor, const Tensor<bool> &mask,
                         const ExecutionContext &execContext);

template <typename T>
Tensor<T> cudaSelectCols(const Tensor<T> &tensor, const Tensor<bool> &mask,
                         const ExecutionContext &execContext);

}  // namespace detail

void sigmoid(Tensor<float> &tensor, std::optional<cudaStream_t> stream) {
  if (tensor.device().type == Device::Type::eCUDA) {
#ifdef ACTS_GNN_WITH_CUDA
    return ActsPlugins::detail::cudaSigmoid(tensor, stream.value());
#else
    throw std::runtime_error(
        "Cannot apply sigmoid to CUDA tensor, library was not compiled with "
        "CUDA");
#endif
  }

  for (auto it = tensor.data(); it != tensor.data() + tensor.size(); ++it) {
    *it = 1.f / (1.f + std::exp(-*it));
  }
}

std::pair<Tensor<std::int64_t>, std::optional<Tensor<float>>> applyEdgeLimit(
    const Tensor<std::int64_t> &edgeIndex,
    const std::optional<Tensor<float>> &edgeFeatures, std::size_t maxEdges,
    std::optional<cudaStream_t> stream) {
  if (edgeFeatures.has_value() &&
      edgeIndex.device() != edgeFeatures->device()) {
    throw std::invalid_argument(
        "limitEdges: edgeIndex and edgeFeatures must be on the same device!");
  }
  if (edgeFeatures.has_value() &&
      edgeFeatures->shape().at(0) != edgeIndex.shape().at(1)) {
    throw std::invalid_argument("limitEdges: inconsistent number of edges");
  }

  const auto nEdgeFeatures =
      edgeFeatures.has_value() ? edgeFeatures->shape().at(1) : 0;
  const auto nEdgesOld = edgeIndex.shape().at(1);

  std::optional<Tensor<std::int64_t>> newEdgeIndexTensor;
  std::optional<Tensor<float>> newEdgeFeatureTensor;

  if (nEdgesOld <= maxEdges) {
    // No need to limit edges, just clone the original tensors
    newEdgeIndexTensor = edgeIndex.clone({edgeIndex.device(), stream});
    if (edgeFeatures.has_value()) {
      newEdgeFeatureTensor =
          edgeFeatures->clone({edgeFeatures->device(), stream});
    }
  } else if (edgeIndex.device().isCpu()) {
    ExecutionContext cpuCtx{Device::Cpu(), {}};

    std::span<const std::int64_t> edge0(edgeIndex.data(), maxEdges);
    std::span<const std::int64_t> edge1(edgeIndex.data() + nEdgesOld, maxEdges);

    newEdgeIndexTensor = Tensor<std::int64_t>::Create({2, maxEdges}, cpuCtx);
    std::copy(edge0.begin(), edge0.end(), newEdgeIndexTensor->data());
    std::copy(edge1.begin(), edge1.end(),
              newEdgeIndexTensor->data() + maxEdges);

    if (edgeFeatures.has_value()) {
      std::span<const float> edgeFeaturesResized(edgeFeatures->data(),
                                                 maxEdges * nEdgeFeatures);

      newEdgeFeatureTensor =
          Tensor<float>::Create({maxEdges, nEdgeFeatures}, cpuCtx);
      std::copy(edgeFeaturesResized.begin(), edgeFeaturesResized.end(),
                newEdgeFeatureTensor->data());
    }
  } else {
#ifdef ACTS_GNN_WITH_CUDA
    ExecutionContext gpuCtx{edgeIndex.device(), stream};

    newEdgeIndexTensor = Tensor<std::int64_t>::Create({2, maxEdges}, gpuCtx);
    ACTS_CUDA_CHECK(cudaMemcpyAsync(newEdgeIndexTensor->data(),
                                    edgeIndex.data(),
                                    maxEdges * sizeof(std::int64_t),
                                    cudaMemcpyDeviceToDevice, stream.value()));
    ACTS_CUDA_CHECK(cudaMemcpyAsync(newEdgeIndexTensor->data() + maxEdges,
                                    edgeIndex.data() + nEdgesOld,
                                    maxEdges * sizeof(std::int64_t),
                                    cudaMemcpyDeviceToDevice, stream.value()));

    if (edgeFeatures.has_value()) {
      newEdgeFeatureTensor =
          Tensor<float>::Create({maxEdges, nEdgeFeatures}, gpuCtx);

      ACTS_CUDA_CHECK(
          cudaMemcpyAsync(newEdgeFeatureTensor->data(), edgeFeatures->data(),
                          maxEdges * nEdgeFeatures * sizeof(float),
                          cudaMemcpyDeviceToDevice, stream.value()));
    }
#else
    throw std::runtime_error(
        "Cannot apply edge limit to CUDA tensors, library was not compiled "
        "with CUDA");
#endif
  }

  return {std::move(newEdgeIndexTensor.value()),
          std::move(newEdgeFeatureTensor)};
}

Tensor<bool> scoreMask(const Tensor<float> &scores, float cut,
                       std::optional<cudaStream_t> stream) {
  if (scores.shape()[1] != 1) {
    throw std::invalid_argument("scoreMask: scores must have shape [N, 1]");
  }
  ExecutionContext execContext{scores.device(), stream};

  if (scores.device().type == Device::Type::eCUDA) {
#ifdef ACTS_GNN_WITH_CUDA
    return detail::cudaScoreMask(scores, cut, stream.value());
#else
    throw std::runtime_error(
        "Cannot compute score mask on CUDA tensor, library was not compiled "
        "with CUDA");
#endif
  }

  auto mask = Tensor<bool>::Create(scores.shape(), execContext);
  for (std::size_t i = 0; i < scores.size(); ++i) {
    mask.data()[i] = scores.data()[i] > cut;
  }
  return mask;
}

template <Acts::Concepts::arithmetic T>
Tensor<T> selectRows(const Tensor<T> &tensor, const Tensor<bool> &mask,
                     const ExecutionContext &execContext) {
  detail::checkMaskCompatibility(tensor, mask, 0);
  const auto nCols = tensor.shape()[1];

  if (tensor.device().type == Device::Type::eCUDA) {
#ifdef ACTS_GNN_WITH_CUDA
    return detail::cudaSelectRows(tensor, mask, execContext);
#else
    throw std::runtime_error(
        "Cannot selectRows on CUDA tensor, library was not compiled with CUDA");
#endif
  }

  const std::size_t n =
      std::count(mask.data(), mask.data() + mask.size(), true);
  auto result = Tensor<T>::Create({n, nCols}, execContext);
  std::size_t out = 0;
  for (std::size_t row = 0; row < tensor.shape()[0]; ++row) {
    if (mask.data()[row]) {
      std::copy_n(tensor.data() + row * nCols, nCols,
                  result.data() + out * nCols);
      ++out;
    }
  }
  return result;
}

template <Acts::Concepts::arithmetic T>
Tensor<T> selectCols(const Tensor<T> &tensor, const Tensor<bool> &mask,
                     const ExecutionContext &execContext) {
  detail::checkMaskCompatibility(tensor, mask, 1);
  const auto nRows = tensor.shape()[0];
  const auto nColsSrc = tensor.shape()[1];

  if (tensor.device().type == Device::Type::eCUDA) {
#ifdef ACTS_GNN_WITH_CUDA
    return detail::cudaSelectCols(tensor, mask, execContext);
#else
    throw std::runtime_error(
        "Cannot selectCols on CUDA tensor, library was not compiled with CUDA");
#endif
  }

  const std::size_t n =
      std::count(mask.data(), mask.data() + mask.size(), true);
  auto result = Tensor<T>::Create({nRows, n}, execContext);
  for (std::size_t row = 0; row < nRows; ++row) {
    std::size_t out = 0;
    for (std::size_t col = 0; col < nColsSrc; ++col) {
      if (mask.data()[col]) {
        result.data()[row * n + out] = tensor.data()[row * nColsSrc + col];
        ++out;
      }
    }
  }
  return result;
}

template Tensor<float> selectRows(const Tensor<float> &, const Tensor<bool> &,
                                  const ExecutionContext &);
template Tensor<std::int64_t> selectRows(const Tensor<std::int64_t> &,
                                         const Tensor<bool> &,
                                         const ExecutionContext &);
template Tensor<float> selectCols(const Tensor<float> &, const Tensor<bool> &,
                                  const ExecutionContext &);
template Tensor<std::int64_t> selectCols(const Tensor<std::int64_t> &,
                                         const Tensor<bool> &,
                                         const ExecutionContext &);

}  // namespace ActsPlugins
