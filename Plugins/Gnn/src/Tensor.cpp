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
#include <format>
#include <fstream>
#include <numeric>
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

std::pair<Tensor<float>, Tensor<std::int64_t>> cudaApplyScoreCut(
    const Tensor<float> &scores, const Tensor<std::int64_t> &edgeIndex,
    float cut, cudaStream_t stream);

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

std::pair<Tensor<float>, Tensor<std::int64_t>> applyScoreCut(
    const Tensor<float> &scores, const Tensor<std::int64_t> &edgeIndex,
    float cut, std::optional<cudaStream_t> stream) {
  assert(scores.shape()[1] == 1);
  assert(edgeIndex.shape()[0] == 2);
  assert(edgeIndex.shape()[1] == scores.shape()[0]);
  assert(scores.device() == edgeIndex.device());
  ExecutionContext execContext{scores.device(), stream};

  if (scores.device().type == Device::Type::eCUDA) {
#ifdef ACTS_GNN_WITH_CUDA
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

void detail::dumpNpy(const std::string &filename, const std::string &type,
                     std::span<const std::byte> data,
                     const std::array<std::size_t, 2> &shape) {
  // Simple NPY header for 2D array
  std::ofstream ofs(filename, std::ios::binary);
  if (!ofs.is_open()) {
    throw std::runtime_error("Could not open file for writing: " + filename);
  }

  // NPY header for version 1.0
  const char vMajor = 1;
  const char vMinor = 0;
  const std::array<char, 8> magicString = {'\x93', 'N', 'U',    'M',
                                           'P',    'Y', vMajor, vMinor};
  ofs.write(magicString.data(), magicString.size());

  // Construct the dictionary
  std::string dict = std::format(
      "{{'descr': '{}', 'fortran_order': False, 'shape': ({}, {}), }}", type,
      shape[0], shape[1]);

  // Pad the dictionary to be 16-byte aligned
  std::size_t padding = 16 - (10 + dict.size()) % 16;
  dict.append(padding, ' ');
  dict.push_back('\n');

  // Write the length of the dictionary
  static_assert(std::endian::native == std::endian::little);
  std::uint16_t dictLen = static_cast<std::uint16_t>(dict.size());
  ofs.write(reinterpret_cast<const char *>(&dictLen), sizeof(dictLen));

  // Write the dictionary
  ofs.write(dict.data(), dict.size());

  // Write the data
  ofs.write(reinterpret_cast<const char *>(data.data()), data.size());
}

}  // namespace ActsPlugins
