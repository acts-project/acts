// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/OnnxEdgeClassifier.hpp"

#include "Acts/Plugins/ExaTrkX/detail/Utils.hpp"

#include <boost/container/static_vector.hpp>
#include <onnxruntime_cxx_api.h>
#include <torch/script.h>

using namespace torch::indexing;
namespace bc = boost::container;

namespace {

Ort::Value torchToOnnx(Ort::MemoryInfo &memoryInfo, at::Tensor &tensor) {
  ONNXTensorElementDataType onnxType = ONNX_TENSOR_ELEMENT_DATA_TYPE_UNDEFINED;

  if (tensor.dtype() == torch::kFloat32) {
    onnxType = ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT;
  } else if (tensor.dtype() == torch::kInt64) {
    onnxType = ONNX_TENSOR_ELEMENT_DATA_TYPE_INT64;
  } else {
    throw std::runtime_error(
        "Cannot convert torch::Tensor to Ort::Value (datatype)");
  }

  bc::static_vector<std::int64_t, 2> shape;
  for (auto size : tensor.sizes()) {
    shape.push_back(size);
  }
  return Ort::Value::CreateTensor(memoryInfo, tensor.data_ptr(),
                                  tensor.nbytes(), shape.data(), shape.size(),
                                  onnxType);
}

}  // namespace

namespace Acts {

OnnxEdgeClassifier::OnnxEdgeClassifier(const Config &cfg,
                                       std::unique_ptr<const Logger> _logger)
    : m_logger(std::move(_logger)), m_cfg(cfg) {
  ACTS_INFO("OnnxEdgeClassifier with ORT API version " << ORT_API_VERSION);

  OrtLoggingLevel onnxLevel = ORT_LOGGING_LEVEL_WARNING;
  switch (m_logger->level()) {
    case Acts::Logging::VERBOSE:
      onnxLevel = ORT_LOGGING_LEVEL_VERBOSE;
      break;
    case Acts::Logging::DEBUG:
      onnxLevel = ORT_LOGGING_LEVEL_INFO;
      break;
    case Acts::Logging::INFO:
      onnxLevel = ORT_LOGGING_LEVEL_WARNING;
      break;
    case Acts::Logging::WARNING:
      onnxLevel = ORT_LOGGING_LEVEL_WARNING;
      break;
    case Acts::Logging::ERROR:
      onnxLevel = ORT_LOGGING_LEVEL_ERROR;
      break;
    case Acts::Logging::FATAL:
      onnxLevel = ORT_LOGGING_LEVEL_FATAL;
      break;
    default:
      throw std::runtime_error("Invalid log level");
  }

  m_env = std::make_unique<Ort::Env>(onnxLevel, "ExaTrkX - edge classifier");

  Ort::SessionOptions sessionOptions;
  sessionOptions.SetIntraOpNumThreads(1);
  sessionOptions.SetGraphOptimizationLevel(
      GraphOptimizationLevel::ORT_ENABLE_EXTENDED);

  if (torch::cuda::is_available()) {
    ACTS_INFO("Try to add ONNX execution provider for CUDA");
    OrtCUDAProviderOptions cuda_options;
    cuda_options.device_id = 0;
    sessionOptions.AppendExecutionProvider_CUDA(cuda_options);
  }

  m_model = std::make_unique<Ort::Session>(*m_env, m_cfg.modelPath.c_str(),
                                           sessionOptions);

  Ort::AllocatorWithDefaultOptions allocator;

  if (m_model->GetInputCount() < 2 || m_model->GetInputCount() > 3) {
    throw std::invalid_argument("ONNX edge classifier needs 2 or 3 inputs!");
  }

  for (std::size_t i = 0; i < m_model->GetInputCount(); ++i) {
    m_inputNames.emplace_back(
        m_model->GetInputNameAllocated(i, allocator).get());
  }

  if (m_model->GetOutputCount() != 1) {
    throw std::invalid_argument(
        "ONNX edge classifier needs exactly one output!");
  }

  m_outputName =
      std::string(m_model->GetOutputNameAllocated(0, allocator).get());
}

OnnxEdgeClassifier::~OnnxEdgeClassifier() {}

std::tuple<std::any, std::any, std::any, std::any>
OnnxEdgeClassifier::operator()(std::any inputNodes, std::any inputEdges,
                               std::any inEdgeFeatures,
                               const ExecutionContext &execContext) {
  const char *deviceStr = "Cpu";
  torch::Device torchDevice(torch::kCPU);
  if (execContext.device.type == Acts::Device::Type::eCUDA) {
    deviceStr = "Cuda";
    torchDevice = torch::Device(torch::kCUDA, execContext.device.index);
  }

  ACTS_DEBUG("Create ORT memory info (" << deviceStr << ")");
  Ort::MemoryInfo memoryInfo(deviceStr, OrtArenaAllocator,
                             execContext.device.index, OrtMemTypeDefault);

  bc::static_vector<Ort::Value, 3> inputTensors;
  bc::static_vector<const char *, 3> inputNames;

  // Node tensor
  auto nodeTensor = std::any_cast<torch::Tensor>(inputNodes).to(torchDevice);
  ACTS_DEBUG("nodes: " << detail::TensorDetails{nodeTensor});
  inputTensors.push_back(torchToOnnx(memoryInfo, nodeTensor));
  inputNames.push_back(m_inputNames.at(0).c_str());

  // Edge tensor
  auto edgeIndex = std::any_cast<torch::Tensor>(inputEdges).to(torchDevice);
  ACTS_DEBUG("edgeIndex: " << detail::TensorDetails{edgeIndex});
  inputTensors.push_back(torchToOnnx(memoryInfo, edgeIndex));
  inputNames.push_back(m_inputNames.at(1).c_str());

  // Edge feature tensor
  std::optional<torch::Tensor> edgeFeatures;
  if (m_inputNames.size() == 3 && inEdgeFeatures.has_value()) {
    edgeFeatures = std::any_cast<torch::Tensor>(inEdgeFeatures).to(torchDevice);
    ACTS_DEBUG("edgeFeatures: " << detail::TensorDetails{*edgeFeatures});
    inputTensors.push_back(torchToOnnx(memoryInfo, *edgeFeatures));
    inputNames.push_back(m_inputNames.at(2).c_str());
  }

  // Output score tensor
  ACTS_DEBUG("Create score tensor");
  auto scores = torch::empty(
      edgeIndex.size(1),
      torch::TensorOptions().device(torchDevice).dtype(torch::kFloat32));
  if (m_model->GetOutputTypeInfo(0)
          .GetTensorTypeAndShapeInfo()
          .GetDimensionsCount() == 2) {
    scores = scores.reshape({scores.numel(), 1});
  }

  std::vector<Ort::Value> outputTensors;
  outputTensors.push_back(torchToOnnx(memoryInfo, scores));
  std::vector<const char *> outputNames{m_outputName.c_str()};

  ACTS_DEBUG("Run model");
  Ort::RunOptions options;
  m_model->Run(options, inputNames.data(), inputTensors.data(),
               inputTensors.size(), outputNames.data(), outputTensors.data(),
               outputNames.size());
  scores = scores.squeeze();

  ACTS_VERBOSE("Slice of classified output before sigmoid:\n"
               << scores.slice(/*dim=*/0, /*start=*/0, /*end=*/9));

  scores.sigmoid_();

  ACTS_DEBUG("scores: " << detail::TensorDetails{scores});
  ACTS_VERBOSE("Slice of classified output:\n"
               << scores.slice(/*dim=*/0, /*start=*/0, /*end=*/9));

  torch::Tensor filterMask = scores > m_cfg.cut;
  torch::Tensor edgesAfterCut = edgeIndex.index({Slice(), filterMask});

  ACTS_DEBUG("Finished edge classification, after cut: "
             << edgesAfterCut.size(1) << " edges.");

  if (edgesAfterCut.size(1) == 0) {
    throw Acts::NoEdgesError{};
  }

  return {std::move(nodeTensor), edgesAfterCut.clone(),
          std::move(inEdgeFeatures), std::move(scores)};
}

}  // namespace Acts
