// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/TorchMetricLearning.hpp"

#include "Acts/Plugins/ExaTrkX/buildEdges.hpp"

#include <torch/script.h>
#include <torch/torch.h>

#include "printCudaMemInfo.hpp"

namespace Acts {

TorchMetricLearning::TorchMetricLearning(const Config &cfg,
                                         std::unique_ptr<const Logger> _logger)
    : m_logger(std::move(_logger)), m_cfg(cfg) {
  c10::InferenceMode guard(true);
  m_deviceType = torch::cuda::is_available() ? torch::kCUDA : torch::kCPU;
  ACTS_DEBUG("Using torch version " << TORCH_VERSION_MAJOR << "."
                                    << TORCH_VERSION_MINOR << "."
                                    << TORCH_VERSION_PATCH);
#ifndef ACTS_EXATRKX_CPUONLY
  if (not torch::cuda::is_available()) {
    ACTS_INFO("CUDA not available, falling back to CPU");
  }
#endif

  try {
    m_model = std::make_unique<torch::jit::Module>();
    *m_model = torch::jit::load(m_cfg.modelPath, m_deviceType);
    m_model->eval();
  } catch (const c10::Error &e) {
    throw std::invalid_argument("Failed to load models: " + e.msg());
  }
}

TorchMetricLearning::~TorchMetricLearning() {}

std::tuple<std::any, std::any> TorchMetricLearning::operator()(
    std::vector<float> &inputValues) {
  ACTS_DEBUG("Start graph construction");
  c10::InferenceMode guard(true);
  const torch::Device device(m_deviceType);

  // printout the r,phi,z of the first spacepoint
  ACTS_VERBOSE("First spacepoint information [r, phi, z]: "
               << inputValues[0] << ", " << inputValues[1] << ", "
               << inputValues[2]);
  ACTS_VERBOSE("Max and min spacepoint: "
               << *std::max_element(inputValues.begin(), inputValues.end())
               << ", "
               << *std::min_element(inputValues.begin(), inputValues.end()))
  printCudaMemInfo(logger());

  // **********
  // Embedding
  // **********

  const int64_t numSpacepoints = inputValues.size() / m_cfg.spacepointFeatures;
  std::vector<torch::jit::IValue> inputTensors;
  auto opts = torch::TensorOptions().dtype(torch::kFloat32);
  torch::Tensor inputTensor =
      torch::from_blob(inputValues.data(),
                       {numSpacepoints, m_cfg.spacepointFeatures}, opts)
          .to(torch::kFloat32)
          .to(device);

  // Clone models (solve memory leak? members can be const...)
  auto model = m_model->clone();
  model.to(device);

  inputTensors.push_back(inputTensor);
  auto output = model.forward(inputTensors).toTensor();
  inputTensors.clear();

  ACTS_VERBOSE("Embedding space of the first SP:\n"
               << output.slice(/*dim=*/0, /*start=*/0, /*end=*/1));
  printCudaMemInfo(logger());

  // ****************
  // Building Edges
  // ****************

  auto edgeList = buildEdges(output, numSpacepoints, m_cfg.embeddingDim,
                             m_cfg.rVal, m_cfg.knnVal);

  ACTS_VERBOSE("Shape of built edges: (" << edgeList.size(0) << ", "
                                         << edgeList.size(1));
  ACTS_VERBOSE("Slice of edgelist:\n" << edgeList.slice(1, 0, 5));
  printCudaMemInfo(logger());

  return {inputTensor, edgeList};
}
}  // namespace Acts
