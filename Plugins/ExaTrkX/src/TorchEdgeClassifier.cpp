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

TorchEdgeClassifier::TorchEdgeClassifier(const Config& cfg,
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
    *m_model = torch::jit::load(m_cfg.modelPath.c_str(), m_deviceType);
    m_model->eval();
  } catch (const c10::Error& e) {
    throw std::invalid_argument("Failed to load models: " + e.msg());
  }
}

TorchEdgeClassifier::~TorchEdgeClassifier() {}

std::tuple<std::any, std::any, std::any> TorchEdgeClassifier::operator()(
    std::any inputNodes, std::any inputEdges, int deviceHint) {
  ACTS_DEBUG("Start edge classification");
  c10::InferenceMode guard(true);
  const torch::Device device(m_deviceType, deviceHint);

  auto nodes = std::any_cast<torch::Tensor>(inputNodes).to(device);
  auto edgeList = std::any_cast<torch::Tensor>(inputEdges).to(device);

  if (m_cfg.numFeatures > nodes.size(1)) {
    throw std::runtime_error("requested more features then available");
  }

  std::vector<at::Tensor> results;
  results.reserve(m_cfg.nChunks);

  auto edgeListTmp =
      m_cfg.undirected ? torch::cat({edgeList, edgeList.flip(0)}, 1) : edgeList;

  std::vector<torch::jit::IValue> inputTensors(2);
  inputTensors[0] = m_cfg.numFeatures < nodes.size(1)
                        ? nodes.index({Slice{}, Slice{None, m_cfg.numFeatures}})
                        : nodes;

  const auto chunks = at::chunk(at::arange(edgeListTmp.size(1)), m_cfg.nChunks);
  for (const auto& chunk : chunks) {
    ACTS_VERBOSE("Process chunk");
    inputTensors[1] = edgeListTmp.index({Slice(), chunk});

    results.push_back(m_model->forward(inputTensors).toTensor());
    results.back().squeeze_();
    results.back().sigmoid_();
  }

  auto output = torch::cat(results);

  if (m_cfg.undirected) {
    output = output.index({Slice(None, output.size(0) / 2)});
  }

  ACTS_VERBOSE("Size after classifier: " << output.size(0));
  ACTS_VERBOSE("Slice of classified output:\n"
               << output.slice(/*dim=*/0, /*start=*/0, /*end=*/9));
  printCudaMemInfo(logger());

  torch::Tensor mask = output > m_cfg.cut;
  torch::Tensor edgesAfterCut = edgeList.index({Slice(), mask});
  edgesAfterCut = edgesAfterCut.to(torch::kInt64);

  ACTS_VERBOSE("Size after score cut: " << edgesAfterCut.size(1));
  printCudaMemInfo(logger());

  return {std::move(nodes), std::move(edgesAfterCut),
          output.masked_select(mask)};
}

}  // namespace Acts
