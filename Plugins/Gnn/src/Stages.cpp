// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Gnn/Stages.hpp"

#include <algorithm>
#include <stdexcept>

namespace ActsPlugins {

namespace {

PipelineTensors cpuRemoveUnusedNodes(PipelineTensors &&tensors,
                                     std::vector<int> &spacePointIds,
                                     const ExecutionContext & /*execCtx*/) {
  const auto nNodes = tensors.nodeFeatures.shape()[0];
  const auto nEdges = tensors.edgeIndex.shape()[1];
  auto *edgeData = tensors.edgeIndex.data();

  std::vector<std::int64_t> sortedIdx(edgeData, edgeData + 2 * nEdges);
  std::ranges::sort(sortedIdx);
  auto dups = std::ranges::unique(sortedIdx);
  sortedIdx.erase(dups.begin(), dups.end());
  const auto nUsed = sortedIdx.size();

  ExecutionContext cpuCtx{Device::Cpu(), {}};
  auto mask = Tensor<bool>::Create({nNodes, 1}, cpuCtx);
  std::fill(mask.data(), mask.data() + nNodes, false);
  std::vector<std::int64_t> oldToNew(nNodes, 0);
  for (std::size_t newIdx = 0; newIdx < nUsed; ++newIdx) {
    const auto oldIdx = static_cast<std::size_t>(sortedIdx[newIdx]);
    mask.data()[oldIdx] = true;
    oldToNew[oldIdx] = static_cast<std::int64_t>(newIdx);
  }

  auto newNodeFeatures = selectRows(tensors.nodeFeatures, mask, cpuCtx);

  for (std::size_t i = 0; i < 2 * nEdges; ++i) {
    edgeData[i] = oldToNew[static_cast<std::size_t>(edgeData[i])];
  }

  std::vector<int> remapped;
  remapped.reserve(nUsed);
  for (const auto oldIdx : sortedIdx) {
    remapped.push_back(spacePointIds[static_cast<std::size_t>(oldIdx)]);
  }
  spacePointIds = std::move(remapped);

  return {std::move(newNodeFeatures), std::move(tensors.edgeIndex),
          std::move(tensors.edgeFeatures), std::move(tensors.edgeScores)};
}

}  // namespace

#ifdef ACTS_GNN_WITH_CUDA
namespace detail {
PipelineTensors cudaRemoveUnusedNodes(PipelineTensors &&tensors,
                                      std::vector<int> &spacePointIds,
                                      const ExecutionContext &execCtx);
}  // namespace detail
#endif

PipelineTensors removeUnusedNodes(PipelineTensors &&tensors,
                                  std::vector<int> &spacePointIds,
                                  const ExecutionContext &execCtx) {
  if (tensors.edgeIndex.shape()[1] == 0) {
    throw NoEdgesError{};
  }

  if (tensors.nodeFeatures.device().isCuda()) {
#ifdef ACTS_GNN_WITH_CUDA
    return detail::cudaRemoveUnusedNodes(std::move(tensors), spacePointIds,
                                         execCtx);
#else
    throw std::runtime_error(
        "Cannot removeUnusedNodes on CUDA tensor, library not compiled with "
        "CUDA");
#endif
  }

  return cpuRemoveUnusedNodes(std::move(tensors), spacePointIds, execCtx);
}

}  // namespace ActsPlugins
