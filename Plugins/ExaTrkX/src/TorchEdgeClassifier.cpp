// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/TorchEdgeClassifier.hpp"

#include <torch/script.h>
#include <torch/torch.h>

#include "printCudaMemInfo.hpp"

using namespace torch::indexing;

namespace Acts {

TorchEdgeClassifier::TorchEdgeClassifier(const Config& cfg) : m_cfg(cfg) {
  c10::InferenceMode guard(true);
  m_deviceType = torch::cuda::is_available() ? torch::kCUDA : torch::kCPU;

  try {
    m_model = std::make_unique<torch::jit::Module>();
    *m_model = torch::jit::load(m_cfg.modelPath.c_str(), m_deviceType);
    m_model->eval();
  } catch (const c10::Error& e) {
    throw std::invalid_argument("Failed to load models: " + e.msg());
  }
}

TorchEdgeClassifier::~TorchEdgeClassifier() {}

std::tuple<std::any, std::any, std::any> TorchEdgeClassifier::operator()(
    std::any inputNodes, std::any inputEdges, const Logger& logger) {
  const torch::Device device(m_deviceType);

  const auto eLibInputTensor = std::any_cast<torch::Tensor>(inputNodes);
  const auto edgeList = std::any_cast<torch::Tensor>(inputEdges);

  c10::InferenceMode guard(true);

  const auto chunks = at::chunk(at::arange(edgeList.size(1)), m_cfg.nChunks);
  std::vector<at::Tensor> results;

  for (const auto& chunk : chunks) {
    std::vector<torch::jit::IValue> fInputTensorJit;
    fInputTensorJit.push_back(eLibInputTensor.to(device));
    fInputTensorJit.push_back(edgeList.index({Slice(), chunk}).to(device));

    results.push_back(m_model->forward(fInputTensorJit).toTensor());
    results.back().squeeze_();
    results.back().sigmoid_();
  }

  auto fOutput = torch::cat(results);
  results.clear();

  ACTS_VERBOSE("Size after filtering network: " << fOutput.size(0));
  ACTS_VERBOSE("Slice of filtered output:\n"
               << fOutput.slice(/*dim=*/0, /*start=*/0, /*end=*/9));
  printCudaMemInfo(logger);

  torch::Tensor filterMask = fOutput > m_cfg.cut;
  torch::Tensor edgesAfterF = edgeList.index({Slice(), filterMask});
  edgesAfterF = edgesAfterF.to(torch::kInt64);
  const int64_t numEdgesAfterF = edgesAfterF.size(1);

  ACTS_VERBOSE("Size after filter cut: " << numEdgesAfterF);
  printCudaMemInfo(logger);

  return {eLibInputTensor, edgesAfterF, fOutput};
}

}  // namespace Acts
