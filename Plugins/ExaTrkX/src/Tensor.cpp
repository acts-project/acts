// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/Tensor.hpp"

namespace Acts {

namespace detail {
void cudaSigmoid(Tensor<float> &tensor, const ExecutionContext &execContext);
std::pair<Tensor<float>, Tensor<std::int64_t>> cudaApplyScoreCut(
    const Tensor<float> &scores, const Tensor<std::int64_t> &edgeIndex,
    float cut, const ExecutionContext &execContext);
}  // namespace detail

void sigmoid(Tensor<float> &tensor, const ExecutionContext &execContext) {
  if (execContext.device.type == Acts::Device::Type::eCUDA) {
#ifdef ACTS_EXATRKX_WITH_CUDA
    return detail::cudaSigmoid(tensor, execContext);
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
    float cut, const ExecutionContext &execContext) {
  assert(scores.shape()[1] == 1);
  assert(edgeIndex.shape()[0] == 2);
  assert(edgeIndex.shape()[1] == scores.shape()[0]);

  if (execContext.device.type == Acts::Device::Type::eCUDA) {
#ifdef ACTS_EXATRKX_WITH_CUDA
    return detail::cudaApplyScoreCut(scores, edgeIndex, cut, execContext);
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
  auto outputScores = Tensor<float>::Create({static_cast<std::size_t>(n), 1},
                                            {Acts::Device::Cpu(), {}});
  auto outputEdges = Tensor<std::int64_t>::Create(
      {2, static_cast<std::size_t>(n)}, {Acts::Device::Cpu(), {}});

  auto scoreIt = outputScores.data();
  auto edgeIt1 = outputEdges.data();
  auto edgeIt2 = outputEdges.data() + n;
  for (auto i : indices) {
    *scoreIt = scores.data()[i];
    *edgeIt1 = edgeIndex.data()[i];
    *edgeIt2 = edgeIndex.data()[i + n];
    ++scoreIt;
    ++edgeIt1;
    ++edgeIt2;
  }

  return {std::move(outputScores), std::move(outputEdges)};
}

}  // namespace Acts
