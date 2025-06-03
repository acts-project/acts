// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/Tensor.hpp"

#ifdef ACTS_EXATRKX_WITH_CUDA
#include "Acts/Plugins/ExaTrkX/detail/CudaUtils.hpp"
#endif

#include <cstring>
#include <numeric>

namespace Acts {

namespace detail {

TensorPtr createTensorMemory(std::size_t nbytes,
                             const ExecutionContext &execContext) {
  if (execContext.device.type == Acts::Device::Type::eCPU) {
    void *ptr = new std::byte[nbytes];
    if (ptr == nullptr) {
      throw std::bad_alloc{};
    }
    return TensorPtr(ptr,
                     [](void *p) { delete[] static_cast<std::byte *>(p); });
  } else {
#ifdef ACTS_EXATRKX_WITH_CUDA
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
#ifdef ACTS_EXATRKX_WITH_CUDA
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

std::pair<Tensor<float>, Tensor<std::int64_t>> cudaApplyScoreCut(
    const Tensor<float> &scores, const Tensor<std::int64_t> &edgeIndex,
    float cut, cudaStream_t stream);

}  // namespace detail

void sigmoid(Tensor<float> &tensor, std::optional<cudaStream_t> stream) {
  if (tensor.device().type == Acts::Device::Type::eCUDA) {
#ifdef ACTS_EXATRKX_WITH_CUDA
    return Acts::detail::cudaSigmoid(tensor, stream.value());
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

std::pair<Tensor<float>, Tensor<std::int64_t>> applyScoreCut(
    const Tensor<float> &scores, const Tensor<std::int64_t> &edgeIndex,
    float cut, std::optional<cudaStream_t> stream) {
  assert(scores.shape()[1] == 1);
  assert(edgeIndex.shape()[0] == 2);
  assert(edgeIndex.shape()[1] == scores.shape()[0]);
  assert(scores.device() == edgeIndex.device());
  ExecutionContext execContext{scores.device(), stream};

  if (scores.device().type == Acts::Device::Type::eCUDA) {
#ifdef ACTS_EXATRKX_WITH_CUDA
    return detail::cudaApplyScoreCut(scores, edgeIndex, cut, stream.value());
#else
    throw std::runtime_error(
        "Cannot apply score cut to CUDA tensor, library was not compiled with "
        "CUDA");
#endif
  }

  std::vector<std::size_t> indices(scores.size());
  std::iota(indices.begin(), indices.end(), 0);
  indices.erase(
      std::remove_if(indices.begin(), indices.end(),
                     [&](std::size_t i) { return scores.data()[i] < cut; }),
      indices.end());
  auto n = indices.size();
  auto outputScores =
      Tensor<float>::Create({static_cast<std::size_t>(n), 1}, execContext);
  auto outputEdges = Tensor<std::int64_t>::Create(
      {2, static_cast<std::size_t>(n)}, execContext);

  auto scoreIt = outputScores.data();
  auto edgeIt1 = outputEdges.data();
  auto edgeIt2 = outputEdges.data() + n;
  for (auto i : indices) {
    *scoreIt = scores.data()[i];
    *edgeIt1 = edgeIndex.data()[i];
    *edgeIt2 = edgeIndex.data()[i + scores.size()];
    ++scoreIt;
    ++edgeIt1;
    ++edgeIt2;
  }

  return {std::move(outputScores), std::move(outputEdges)};
}

}  // namespace Acts
