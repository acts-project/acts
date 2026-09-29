// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Gnn/Tensor.hpp"

#include <cstring>
#include <span>

#include "DeviceOps.hpp"

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
  }
  return cudaCreateTensorMemory(nbytes, execContext);
}

TensorPtr cloneTensorMemory(const TensorPtr &ptr, std::size_t nbytes,
                            Device devFrom, const ExecutionContext &to) {
  auto clone = createTensorMemory(nbytes, to);
  if (devFrom.isCpu() && to.device.isCpu()) {
    std::memcpy(clone.get(), ptr.get(), nbytes);
  } else {
    cudaCopyTensorMemory(clone.get(), ptr.get(), nbytes, devFrom, to);
  }
  return clone;
}

}  // namespace detail

void sigmoid(Tensor<float> &tensor, std::optional<cudaStream_t> stream) {
  if (tensor.device().type == Device::Type::eCUDA) {
    return detail::cudaSigmoid(tensor, stream.value());
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
    ExecutionContext gpuCtx{edgeIndex.device(), stream};
    const Device devFrom = edgeIndex.device();

    newEdgeIndexTensor = Tensor<std::int64_t>::Create({2, maxEdges}, gpuCtx);
    detail::cudaCopyTensorMemory(newEdgeIndexTensor->data(), edgeIndex.data(),
                                 maxEdges * sizeof(std::int64_t), devFrom,
                                 gpuCtx);
    detail::cudaCopyTensorMemory(
        newEdgeIndexTensor->data() + maxEdges, edgeIndex.data() + nEdgesOld,
        maxEdges * sizeof(std::int64_t), devFrom, gpuCtx);

    if (edgeFeatures.has_value()) {
      newEdgeFeatureTensor =
          Tensor<float>::Create({maxEdges, nEdgeFeatures}, gpuCtx);

      detail::cudaCopyTensorMemory(
          newEdgeFeatureTensor->data(), edgeFeatures->data(),
          maxEdges * nEdgeFeatures * sizeof(float), devFrom, gpuCtx);
    }
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
    return detail::cudaScoreMask(scores, cut, stream.value());
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
    return detail::cudaSelectRows(tensor, mask, execContext);
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
    return detail::cudaSelectCols(tensor, mask, execContext);
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

template <Acts::Concepts::arithmetic T>
Tensor<T> mulPerColumn(const Tensor<T> &src, const std::vector<T> &scales,
                       const ExecutionContext &execContext) {
  const auto rows = src.shape()[0];
  const auto cols = src.shape()[1];
  if (cols != static_cast<std::size_t>(scales.size())) {
    throw std::invalid_argument("mulPerColumn: scales size must match cols");
  }

  if (execContext.device.type == Device::Type::eCUDA) {
    // Move scales to device (as Tensor)
    auto scalesTensor = Tensor<T>::Create({cols, 1ul}, execContext);
    detail::cudaCopyTensorMemory(scalesTensor.data(), scales.data(),
                                 cols * sizeof(T), Device::Cpu(), execContext);
    return detail::cudaMulPerColumn(src, scalesTensor, execContext);
  }

  auto result = Tensor<T>::Create({rows, cols}, execContext);
  for (std::size_t r = 0; r < rows; ++r) {
    const std::size_t base = r * cols;
    for (std::size_t c = 0; c < cols; ++c) {
      result.data()[base + c] = src.data()[base + c] * scales[c];
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
template Tensor<float> mulPerColumn(const Tensor<float> &,
                                    const std::vector<float> &,
                                    const ExecutionContext &);
template Tensor<std::int64_t> mulPerColumn(const Tensor<std::int64_t> &,
                                           const std::vector<std::int64_t> &,
                                           const ExecutionContext &);

}  // namespace ActsPlugins
