// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/ModuleMapCpp.hpp"

#include "Acts/Plugins/ExaTrkX/detail/GraphCreatorWrapper.hpp"
#include "Acts/Plugins/ExaTrkX/detail/TensorVectorConversion.hpp"
#include "Acts/Plugins/ExaTrkX/detail/Utils.hpp"

#include <TTree_hits>
#include <chrono>
#include <graph_creator>
#include <map>
#include <module_map_triplet>

using namespace torch::indexing;

namespace Acts {

ModuleMapCpp::ModuleMapCpp(const Config &cfg,
                           std::unique_ptr<const Acts::Logger> logger_)
    : m_cfg(cfg), m_logger(std::move(logger_)) {
  if (!m_cfg.useGpu) {
    m_graphCreator =
        std::make_unique<detail::GraphCreatorWrapperCpu>(m_cfg.moduleMapPath);
  } else {
#ifndef ACTS_EXATRKX_CPUONLY
    m_graphCreator = std::make_unique<detail::GraphCreatorWrapperCuda>(
        m_cfg.moduleMapPath, m_cfg.gpuDevice, m_cfg.gpuBlocks);
#else
    throw std::runtime_error(
        "Cannot use cuda version of GraphModuleMap (CUDA is not enabled in "
        "CMake)");
#endif
  }
}

ModuleMapCpp::~ModuleMapCpp() {}

std::tuple<std::any, std::any, std::any> ModuleMapCpp::operator()(
    std::vector<float> &inputValues, std::size_t numNodes,
    const std::vector<std::uint64_t> &moduleIds, torch::Device /*device*/) {
  if (numNodes != moduleIds.size()) {
    throw std::invalid_argument(
        "Module Ids do not match number of graph nodes");
  }

  if (inputValues.empty()) {
    throw NoEdgesError{};
  }

  const auto numFeatures = inputValues.size() / numNodes;

  assert(inputValues.size() % numFeatures == 0);
  auto nodeFeatures =
      detail::vectorToTensor2D(inputValues, numFeatures).clone();
  ACTS_DEBUG("nodeFeatures: " << detail::TensorDetails{nodeFeatures});

  auto [edgeIndex, edgeFeatures] =
      m_graphCreator->build(inputValues, moduleIds, logger());

  if (edgeIndex.numel() == 0 || (edgeIndex == 0).all().item<bool>()) {
    throw Acts::NoEdgesError{};
  }

  return std::make_tuple(std::move(nodeFeatures), std::move(edgeIndex),
                         std::move(edgeFeatures));
}

}  // namespace Acts
