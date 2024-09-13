// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/TorchEdgeClassifierAOT.hpp"

#include "Acts/Plugins/ExaTrkX/detail/Utils.hpp"

#include <chrono>

#ifdef ACTS_EXATRKX_CPUONLY
#error "Cannot build in CPUONLY mode"
#endif

#include <c10/cuda/CUDAGuard.h>
#include <torch/script.h>
#include <torch/torch.h>

#include "printCudaMemInfo.hpp"

using namespace torch::indexing;

namespace Acts {

TorchEdgeClassifierAOT::TorchEdgeClassifierAOT(
    const Config& cfg, std::unique_ptr<const Logger> _logger)
    : m_logger(std::move(_logger)),
      m_cfg(cfg),
      m_device(torch::Device(torch::kCPU)) {
  c10::InferenceMode guard(true);
  if (not torch::cuda::is_available()) {
    throw std::runtime_error("CUDA is required for TorchEdgeClassifierAOT!");
  }

  ACTS_DEBUG("Using torch version " << TORCH_VERSION_MAJOR << "."
                                    << TORCH_VERSION_MINOR << "."
                                    << TORCH_VERSION_PATCH);

  std::make_unique<torch::inductor::AOTIModelContainerRunnerCuda>(
      m_cfg.modelPath.c_str());
}

TorchEdgeClassifierAOT::~TorchEdgeClassifierAOT() {}

std::tuple<std::any, std::any, std::any, std::any>
TorchEdgeClassifierAOT::operator()(std::any inNodeFeatures,
                                   std::any inEdgeIndex,
                                   std::any inEdgeFeatures,
                                   torch::Device device) {
  decltype(std::chrono::high_resolution_clock::now()) t0, t1, t2, t3, t4, t5;
  t0 = std::chrono::high_resolution_clock::now();
  ACTS_DEBUG("Start edge classification, use " << device);
  c10::InferenceMode guard(true);

  // add a protection to avoid calling for kCPU
  std::optional<c10::cuda::CUDAGuard> device_guard;
  if (device.is_cuda()) {
    device_guard.emplace(device.index());
  }

  auto nodeFeatures = std::any_cast<torch::Tensor>(inNodeFeatures).to(device);
  auto edgeIndex = std::any_cast<torch::Tensor>(inEdgeIndex).to(device);

  ACTS_DEBUG("edgeIndex: " << detail::TensorDetails{edgeIndex});

  std::optional<torch::Tensor> edgeFeatures;
  if (inEdgeFeatures.has_value()) {
    edgeFeatures = std::any_cast<torch::Tensor>(inEdgeFeatures).to(device);
    ACTS_DEBUG("edgeFeatures: " << detail::TensorDetails{*edgeFeatures});
  }

  t2 = std::chrono::high_resolution_clock::now();

  torch::Tensor output;

  // Scope this to keep inference objects separate
  {
    std::vector<torch::Tensor> inputs;
    auto selectedFeaturesTensor =
        at::tensor(at::ArrayRef<int>(m_cfg.selectedFeatures));

    inputs.emplace_back(
        !m_cfg.selectedFeatures.empty()
            ? nodeFeatures.index({Slice{}, selectedFeaturesTensor}).clone()
            : nodeFeatures);

    ACTS_DEBUG(
        "selected nodeFeatures: " << detail::TensorDetails{inputs.back()});

    inputs.emplace_back(m_cfg.undirected
                            ? torch::cat({edgeIndex, edgeIndex.flip(0)}, 1)
                            : edgeIndex);

    if (edgeFeatures && m_cfg.useEdgeFeatures) {
      inputs.push_back(*edgeFeatures);
    }

    t3 = std::chrono::high_resolution_clock::now();
    output = m_model->run(inputs).at(0).to(torch::kFloat32);
    t4 = std::chrono::high_resolution_clock::now();
    output.squeeze_();
  }

  output.sigmoid_();

  if (m_cfg.undirected) {
    auto newSize = output.size(0) / 2;
    output = output.index({Slice(None, newSize)});
  }

  ACTS_VERBOSE("Size after classifier: " << output.size(0));
  ACTS_VERBOSE("Slice of classified output:\n"
               << output.slice(/*dim=*/0, /*start=*/0, /*end=*/9));
  printCudaMemInfo(logger());

  torch::Tensor mask = output > m_cfg.cut;
  torch::Tensor edgesAfterCut = edgeIndex.index({Slice(), mask});
  edgesAfterCut = edgesAfterCut.to(torch::kInt64);

  ACTS_VERBOSE("Size after score cut: " << edgesAfterCut.size(1));
  printCudaMemInfo(logger());
  t5 = std::chrono::high_resolution_clock::now();

  auto milliseconds = [](const auto& a, const auto& b) {
    return std::chrono::duration<double, std::milli>(b - a).count();
  };
  ACTS_DEBUG("Time anycast, device guard:  " << milliseconds(t0, t2));
  ACTS_DEBUG("Time jit::IValue creation:   " << milliseconds(t2, t3));
  ACTS_DEBUG("Time model forward:          " << milliseconds(t3, t4));
  ACTS_DEBUG("Time sigmoid and cut:        " << milliseconds(t4, t5));

  return {std::move(nodeFeatures), std::move(edgesAfterCut),
          std::move(inEdgeFeatures), output.masked_select(mask)};
}

}  // namespace Acts
