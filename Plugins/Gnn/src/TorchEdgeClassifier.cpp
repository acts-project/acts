// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Gnn/TorchEdgeClassifier.hpp"

#include "ActsPlugins/Gnn/detail/TensorVectorConversion.hpp"
#include "ActsPlugins/Gnn/detail/Utils.hpp"

#include <chrono>

#ifdef ACTS_GNN_WITH_CUDA
#include <c10/cuda/CUDAGuard.h>
#endif

#include <torch/script.h>
#include <torch/torch.h>

#include "printCudaMemInfo.hpp"

using namespace torch::indexing;

using namespace Acts;

namespace ActsPlugins {

TorchEdgeClassifier::TorchEdgeClassifier(const Config& cfg,
                                         std::unique_ptr<const Logger> _logger)
    : m_logger(std::move(_logger)), m_cfg(cfg) {
  c10::InferenceMode guard(true);
  torch::Device device = torch::kCPU;

  if (cfg.device.isCuda()) {
    if (!torch::cuda::is_available()) {
      throw std::runtime_error(
          "CUDA device requested but CUDA is not available");
    }
    if (cfg.device.index >=
        static_cast<std::size_t>(torch::cuda::device_count())) {
      throw std::runtime_error(
          "CUDA device index " + std::to_string(cfg.device.index) +
          " is out of range (" + std::to_string(torch::cuda::device_count()) +
          " devices available)");
    }
    device = torch::Device(torch::kCUDA, cfg.device.index);
  }

  ACTS_DEBUG("Using torch version " << TORCH_VERSION_MAJOR << "."
                                    << TORCH_VERSION_MINOR << "."
                                    << TORCH_VERSION_PATCH);

  try {
    m_model = std::make_unique<torch::jit::Module>();
    *m_model = torch::jit::load(m_cfg.modelPath, device);
    m_model->eval();
  } catch (const c10::Error& e) {
    throw std::invalid_argument("Failed to load models: " + e.msg());
  }
}

TorchEdgeClassifier::~TorchEdgeClassifier() {}

PipelineTensors TorchEdgeClassifier::operator()(
    PipelineTensors tensors, const ExecutionContext& execContext) {
  const auto device =
      execContext.device.type == Device::Type::eCUDA
          ? torch::Device(torch::kCUDA, execContext.device.index)
          : torch::kCPU;
  decltype(std::chrono::high_resolution_clock::now()) t0, t1, t2, t3, t4;
  t0 = std::chrono::high_resolution_clock::now();
  ACTS_DEBUG("Start edge classification, use " << device);

  if (tensors.edgeIndex.size() == 0) {
    throw NoEdgesError{};
  }

  c10::InferenceMode guard(true);

  // add a protection to avoid calling for kCPU
#ifndef ACTS_GNN_WITH_CUDA
  assert(device == torch::Device(torch::kCPU));
#else
  std::optional<c10::cuda::CUDAGuard> device_guard;
  if (device.is_cuda()) {
    device_guard.emplace(device.index());
  }
#endif

  auto nodeFeatures = detail::actsToNonOwningTorchTensor(tensors.nodeFeatures);
  ACTS_DEBUG("nodeFeatures: " << detail::TensorDetails{nodeFeatures});

  auto edgeIndex = detail::actsToNonOwningTorchTensor(tensors.edgeIndex);
  ACTS_DEBUG("edgeIndex: " << detail::TensorDetails{edgeIndex});

  std::optional<torch::Tensor> edgeFeatures;
  if (tensors.edgeFeatures.has_value()) {
    edgeFeatures = detail::actsToNonOwningTorchTensor(*tensors.edgeFeatures);
    ACTS_DEBUG("edgeFeatures: " << detail::TensorDetails{*edgeFeatures});
  }

  torch::Tensor output;

  // Scope this to keep inference objects separate
  {
    auto edgeIndexTmp = m_cfg.undirected
                            ? torch::cat({edgeIndex, edgeIndex.flip(0)}, 1)
                            : edgeIndex;

    std::vector<torch::jit::IValue> inputTensors(2);
    auto selectedFeaturesTensor =
        at::tensor(at::ArrayRef<int>(m_cfg.selectedFeatures));
    at::Tensor selectedNodeFeatures =
        !m_cfg.selectedFeatures.empty()
            ? nodeFeatures.index({Slice{}, selectedFeaturesTensor}).clone()
            : nodeFeatures;

    ACTS_DEBUG("selected nodeFeatures: "
               << detail::TensorDetails{selectedNodeFeatures});
    inputTensors[0] = selectedNodeFeatures;

    if (edgeFeatures && m_cfg.useEdgeFeatures) {
      inputTensors.push_back(*edgeFeatures);
    }

    t1 = std::chrono::high_resolution_clock::now();

    if (m_cfg.nChunks > 1) {
      std::vector<at::Tensor> results;
      results.reserve(m_cfg.nChunks);

      auto chunks = at::chunk(edgeIndexTmp, m_cfg.nChunks, 1);
      for (auto& chunk : chunks) {
        ACTS_VERBOSE("Process chunk with shape" << chunk.sizes());
        inputTensors[1] = chunk;

        results.push_back(m_model->forward(inputTensors).toTensor());
        results.back().squeeze_();
      }

      output = torch::cat(results);
    } else {
      inputTensors[1] = edgeIndexTmp;
      output = m_model->forward(inputTensors).toTensor().to(torch::kFloat32);
      output.squeeze_();
    }
    t2 = std::chrono::high_resolution_clock::now();
  }

  ACTS_VERBOSE("Slice of classified output before sigmoid:\n"
               << output.slice(/*dim=*/0, /*start=*/0, /*end=*/9));

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

  if (edgesAfterCut.numel() == 0) {
    throw NoEdgesError{};
  }

  ACTS_VERBOSE("Size after score cut: " << edgesAfterCut.size(1));
  printCudaMemInfo(logger());

  std::optional<Tensor<float>> filteredEdgeFeatures;
  if (edgeFeatures.has_value()) {
    auto filtered = edgeFeatures->index({mask, Slice()}).contiguous();
    filteredEdgeFeatures =
        detail::torchToActsTensor<float>(filtered, execContext);
  }

  t3 = std::chrono::high_resolution_clock::now();

  auto milliseconds = [](const auto& a, const auto& b) {
    return std::chrono::duration<double, std::milli>(b - a).count();
  };
  ACTS_DEBUG("Time preparation:    " << milliseconds(t0, t1));
  ACTS_DEBUG("Time inference:      " << milliseconds(t1, t2));
  ACTS_DEBUG("Time postprocessing: " << milliseconds(t2, t3));

  return {std::move(tensors.nodeFeatures),
          detail::torchToActsTensor<std::int64_t>(edgesAfterCut, execContext),
          std::move(filteredEdgeFeatures),
          detail::torchToActsTensor<float>(output.masked_select(mask),
                                           execContext)};
}

}  // namespace ActsPlugins
