// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/TorchMetricLearning.hpp"

#include <torch/script.h>

#include "buildEdges.hpp"

namespace Acts {

TorchMetricLearning::TorchMetricLearning(const Config &cfg) : m_cfg(cfg) {
  c10::InferenceMode guard(true);

  try {
    m_model = std::make_unique<torch::jit::Module>();
    *m_model = torch::jit::load(m_cfg.modelPath);
    m_model->eval();
  } catch (const c10::Error &e) {
    throw std::invalid_argument("Failed to load models: " + e.msg());
  }
}

TorchMetricLearning::~TorchMetricLearning() {}

std::tuple<std::any, std::any> TorchMetricLearning::operator()(
    std::vector<float> &inputValues, const Logger &logger) {
  c10::InferenceMode guard(true);
  torch::Device device(torch::kCUDA);

  // Clone models (solve memory leak? members can be const...)
  auto e_model = m_model->clone();
  e_model.to(device);

  // printout the r,phi,z of the first spacepoint
  ACTS_VERBOSE("First spacepoint information [r, phi, z]: "
               << inputValues[0] << ", " << inputValues[1] << ", "
               << inputValues[2]);
  ACTS_VERBOSE("Max and min spacepoint: "
               << *std::max_element(inputValues.begin(), inputValues.end())
               << ", "
               << *std::min_element(inputValues.begin(), inputValues.end()))
  // print_current_cuda_meminfo(logger);

  // ExaTrkXTimer timer(not recordTiming);

  // **********
  // Embedding
  // **********

  // timer.start();
  int64_t numSpacepoints = inputValues.size() / m_cfg.spacepointFeatures;
  std::vector<torch::jit::IValue> eInputTensorJit;
  auto e_opts = torch::TensorOptions().dtype(torch::kFloat32);
  torch::Tensor eLibInputTensor =
      torch::from_blob(inputValues.data(),
                       {numSpacepoints, m_cfg.spacepointFeatures}, e_opts)
          .to(torch::kFloat32);

  eInputTensorJit.push_back(eLibInputTensor.to(device));
  std::optional<at::Tensor> eOutput =
      e_model.forward(eInputTensorJit).toTensor();
  eInputTensorJit.clear();

  ACTS_VERBOSE("Embedding space of the first SP:\n"
               << eOutput->slice(/*dim=*/0, /*start=*/0, /*end=*/1));
  // print_current_cuda_meminfo(logger);

  // timeInfo.embedding = timer.stopAndGetElapsedTime();

  // ****************
  // Building Edges
  // ****************

  // timer.start();

  // At this point, buildEdgesBruteForce could be used instead
  std::optional<torch::Tensor> edgeList = buildEdges(
      *eOutput, numSpacepoints, m_cfg.embeddingDim, m_cfg.rVal, m_cfg.knnVal);
  eOutput.reset();

  ACTS_VERBOSE("Shape of built edges: (" << edgeList->size(0) << ", "
                                         << edgeList->size(1));
  ACTS_VERBOSE("Slice of edgelist:\n" << edgeList->slice(1, 0, 5));
  // print_current_cuda_meminfo(logger);

  // timeInfo.building = timer.stopAndGetElapsedTime();
  return {eLibInputTensor, *edgeList};
}
}  // namespace Acts
