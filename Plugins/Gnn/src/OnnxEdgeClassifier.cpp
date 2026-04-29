// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Gnn/OnnxEdgeClassifier.hpp"

#include <boost/container/static_vector.hpp>
#include <onnxruntime_cxx_api.h>

namespace bc = boost::container;

namespace {

template <typename T>
Ort::Value toOnnx(Ort::MemoryInfo &memoryInfo, ActsPlugins::Tensor<T> &tensor,
                  std::size_t rank = 2) {
  assert(rank == 1 || rank == 2);

  bc::static_vector<std::int64_t, 2> shape;
  for (auto size : tensor.shape()) {
    // If rank is 1 and we encounter a dimension with size 1, then we skip it
    if (size > 1 || rank == 2) {
      shape.push_back(size);
    }
  }

  assert(shape.size() == rank);
  return Ort::Value::CreateTensor<T>(memoryInfo, tensor.data(), tensor.size(),
                                     shape.data(), shape.size());
}

}  // namespace

using namespace Acts;

namespace ActsPlugins {

OnnxEdgeClassifier::OnnxEdgeClassifier(const Config &cfg,
                                       std::unique_ptr<const Logger> _logger)
    : m_logger(std::move(_logger)), m_cfg(cfg) {
  ACTS_INFO("OnnxEdgeClassifier with ORT API version " << ORT_API_VERSION);

  OrtLoggingLevel onnxLevel = ORT_LOGGING_LEVEL_WARNING;
  switch (m_logger->level()) {
    case Logging::VERBOSE:
      onnxLevel = ORT_LOGGING_LEVEL_VERBOSE;
      break;
    case Logging::DEBUG:
      onnxLevel = ORT_LOGGING_LEVEL_INFO;
      break;
    case Logging::INFO:
      onnxLevel = ORT_LOGGING_LEVEL_WARNING;
      break;
    case Logging::WARNING:
      onnxLevel = ORT_LOGGING_LEVEL_WARNING;
      break;
    case Logging::ERROR:
      onnxLevel = ORT_LOGGING_LEVEL_ERROR;
      break;
    case Logging::FATAL:
      onnxLevel = ORT_LOGGING_LEVEL_FATAL;
      break;
    default:
      throw std::runtime_error("Invalid log level");
  }

  m_env = std::make_unique<Ort::Env>(onnxLevel, "Gnn - edge classifier");

  Ort::SessionOptions sessionOptions;
  sessionOptions.SetIntraOpNumThreads(1);
  sessionOptions.SetGraphOptimizationLevel(
      GraphOptimizationLevel::ORT_ENABLE_EXTENDED);
  sessionOptions.SetExecutionMode(ORT_SEQUENTIAL);

  if (m_cfg.device.isCuda()) {
    ACTS_INFO("Try to add ONNX execution provider for CUDA");
    OrtCUDAProviderOptions cuda_options;
    cuda_options.device_id = m_cfg.device.index;
    sessionOptions.AppendExecutionProvider_CUDA(cuda_options);
  } else {
    ACTS_INFO("Using CPU execution provider for ONNX");
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

PipelineTensors OnnxEdgeClassifier::operator()(
    PipelineTensors tensors, const ExecutionContext &execContext) {
  auto memoryInfo =
      tensors.nodeFeatures.device().isCpu()
          ? Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault)
          : Ort::MemoryInfo("Cuda", OrtArenaAllocator, execContext.device.index,
                            OrtMemTypeDefault);

  bc::static_vector<Ort::Value, 3> inputTensors;
  bc::static_vector<const char *, 3> inputNames;

  // Node tensor
  inputTensors.push_back(toOnnx(memoryInfo, tensors.nodeFeatures));
  inputNames.push_back(m_inputNames.at(0).c_str());
  ACTS_DEBUG("Node features shape: (" << tensors.nodeFeatures.shape()[0] << ", "
                                      << tensors.nodeFeatures.shape()[1]
                                      << ")");

  // Edge tensor
  inputTensors.push_back(toOnnx(memoryInfo, tensors.edgeIndex));
  inputNames.push_back(m_inputNames.at(1).c_str());
  ACTS_DEBUG("Edge index shape: (" << tensors.edgeIndex.shape()[0] << ", "
                                   << tensors.edgeIndex.shape()[1] << ")");

  // If the model has three inputs, we require edge features, otherwise throw
  if (m_inputNames.size() == 3 && !tensors.edgeFeatures.has_value()) {
    throw std::invalid_argument(
        "ONNX edge classifier model has three inputs, but no edge features "
        "provided!");
  }

  // Edge feature tensor
  if (m_inputNames.size() == 3 && tensors.edgeFeatures.has_value()) {
    inputTensors.push_back(toOnnx(memoryInfo, *tensors.edgeFeatures));
    inputNames.push_back(m_inputNames.at(2).c_str());
    ACTS_DEBUG("Edge features shape: ("
               << tensors.edgeFeatures->shape()[0] << ", "
               << tensors.edgeFeatures->shape()[1] << ")");
  }

  // Output score tensor
  ACTS_DEBUG("Create score tensor");
  auto scores =
      Tensor<float>::Create({tensors.edgeIndex.shape()[1], 1ul}, execContext);

  std::vector<Ort::Value> outputTensors;
  auto outputRank = m_model->GetOutputTypeInfo(0)
                        .GetTensorTypeAndShapeInfo()
                        .GetDimensionsCount();
  outputTensors.push_back(toOnnx(memoryInfo, scores, outputRank));
  std::vector<const char *> outputNames{m_outputName.c_str()};

  ACTS_DEBUG("Run model");
  Ort::RunOptions options;
  m_model->Run(options, inputNames.data(), inputTensors.data(),
               inputTensors.size(), outputNames.data(), outputTensors.data(),
               outputNames.size());

  sigmoid(scores, execContext.stream);
  auto mask = scoreMask(scores, m_cfg.cut, execContext.stream);
  auto newScores = selectRows(scores, mask, execContext);
  auto newEdgeIndex = selectCols(tensors.edgeIndex, mask, execContext);
  std::optional<Tensor<float>> newEdgeFeatures;
  if (tensors.edgeFeatures.has_value()) {
    newEdgeFeatures = selectRows(*tensors.edgeFeatures, mask, execContext);
  }

  ACTS_DEBUG("Finished edge classification, after cut: "
             << newEdgeIndex.shape()[1] << " edges.");

  if (newEdgeIndex.shape()[1] == 0) {
    throw NoEdgesError{};
  }

  return {std::move(tensors.nodeFeatures), std::move(newEdgeIndex),
          std::move(newEdgeFeatures), std::move(newScores)};
}

}  // namespace ActsPlugins
