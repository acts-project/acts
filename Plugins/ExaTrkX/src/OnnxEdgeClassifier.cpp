// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/OnnxEdgeClassifier.hpp"

#include <boost/container/static_vector.hpp>
#include <onnxruntime_cxx_api.h>

namespace bc = boost::container;

namespace {

template <typename T>
Ort::Value toOnnx(Ort::MemoryInfo &memoryInfo, Acts::Tensor<T> &tensor,
                  std::size_t rank = 2) {
  assert(rank == 1 || rank == 2);
  ONNXTensorElementDataType onnxType = ONNX_TENSOR_ELEMENT_DATA_TYPE_UNDEFINED;

  if constexpr (std::is_same_v<T, float>) {
    onnxType = ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT;
  } else if constexpr (std::is_same_v<T, std::int64_t>) {
    onnxType = ONNX_TENSOR_ELEMENT_DATA_TYPE_INT64;
  } else {
    throw std::runtime_error(
        "Cannot convert Acts::Tensor to Ort::Value (datatype)");
  }

  bc::static_vector<std::int64_t, 2> shape;
  for (auto size : tensor.shape()) {
    // If rank is 1 and we encounter a dimension with size 1, then we skip it
    if (size > 1 || rank == 2) {
      shape.push_back(size);
    }
  }

  assert(shape.size() == rank);
  return Ort::Value::CreateTensor(memoryInfo, tensor.data(), tensor.nbytes(),
                                  shape.data(), shape.size(), onnxType);
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

#ifndef ACTS_EXATRKX_CPUONLY
  ACTS_INFO("Try to add ONNX execution provider for CUDA");
  OrtCUDAProviderOptions cuda_options;
  cuda_options.device_id = 0;
  sessionOptions.AppendExecutionProvider_CUDA(cuda_options);
#endif

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
  const char *deviceStr = "Cpu";
  if (execContext.device.type == Acts::Device::Type::eCUDA) {
    deviceStr = "Cuda";
  }

  ACTS_DEBUG("Create ORT memory info (" << deviceStr << ")");
  Ort::MemoryInfo memoryInfo(deviceStr, OrtArenaAllocator,
                             execContext.device.index, OrtMemTypeDefault);

  bc::static_vector<Ort::Value, 3> inputTensors;
  bc::static_vector<const char *, 3> inputNames;

  // Node tensor
  inputTensors.push_back(toOnnx(memoryInfo, tensors.nodeFeatures));
  inputNames.push_back(m_inputNames.at(0).c_str());

  // Edge tensor
  inputTensors.push_back(toOnnx(memoryInfo, tensors.edgeIndex));
  inputNames.push_back(m_inputNames.at(1).c_str());

  // Edge feature tensor
  std::optional<Acts::Tensor<float>> edgeFeatures;
  if (m_inputNames.size() == 3 && tensors.edgeFeatures.has_value()) {
    inputTensors.push_back(toOnnx(memoryInfo, *tensors.edgeFeatures));
    inputNames.push_back(m_inputNames.at(2).c_str());
  }

  // Output score tensor
  ACTS_DEBUG("Create score tensor");
  auto scores = Acts::Tensor<float>::Create({tensors.edgeIndex.shape()[1], 1ul},
                                            execContext);

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
  auto [newScores, newEdgeIndex] =
      applyScoreCut(scores, tensors.edgeIndex, m_cfg.cut, execContext.stream);

  ACTS_DEBUG("Finished edge classification, after cut: "
             << newEdgeIndex.shape()[1] << " edges.");

  if (newEdgeIndex.shape()[1] == 0) {
    throw Acts::NoEdgesError{};
  }

  return {std::move(tensors.nodeFeatures),
          std::move(newEdgeIndex),
          {},
          std::move(newScores)};
}

}  // namespace Acts
