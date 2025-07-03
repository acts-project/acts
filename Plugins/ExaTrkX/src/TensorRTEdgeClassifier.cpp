// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/TensorRTEdgeClassifier.hpp"

#include "Acts/Plugins/ExaTrkX/detail/CudaUtils.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>

#include <NvInfer.h>
#include <NvInferPlugin.h>
#include <NvInferRuntimeBase.h>
#include <cuda_runtime.h>

#include "printCudaMemInfo.hpp"

namespace {

class TensorRTLogger : public nvinfer1::ILogger {
  std::unique_ptr<const Acts::Logger> m_logger;

 public:
  TensorRTLogger(Acts::Logging::Level lvl)
      : m_logger(Acts::getDefaultLogger("TensorRT", lvl)) {}

  void log(Severity severity, const char *msg) noexcept override {
    const auto &logger = *m_logger;
    switch (severity) {
      case Severity::kVERBOSE:
        ACTS_DEBUG(msg);
        break;
      case Severity::kINFO:
        ACTS_INFO(msg);
        break;
      case Severity::kWARNING:
        ACTS_WARNING(msg);
        break;
      case Severity::kERROR:
        ACTS_ERROR(msg);
        break;
      case Severity::kINTERNAL_ERROR:
        ACTS_FATAL(msg);
        break;
    }
  }
};

}  // namespace

namespace Acts {

TensorRTEdgeClassifier::TensorRTEdgeClassifier(
    const Config &cfg, std::unique_ptr<const Logger> _logger)
    : m_logger(std::move(_logger)),
      m_cfg(cfg),
      m_trtLogger(std::make_unique<TensorRTLogger>(m_logger->level())) {
  auto status = initLibNvInferPlugins(m_trtLogger.get(), "");
  if (!status) {
    throw std::runtime_error("Failed to initialize TensorRT plugins");
  }

  std::size_t fsize =
      std::filesystem::file_size(std::filesystem::path(m_cfg.modelPath));
  std::vector<char> engineData(fsize);

  ACTS_DEBUG("Load '" << m_cfg.modelPath << "' with size " << fsize);

  std::ifstream engineFile(m_cfg.modelPath);
  if (!engineFile) {
    throw std::runtime_error("Failed to open engine file");
  } else if (!engineFile.read(engineData.data(), fsize)) {
    throw std::runtime_error("Failed to read engine file");
  }

  m_runtime.reset(nvinfer1::createInferRuntime(*m_trtLogger));
  if (!m_runtime) {
    throw std::runtime_error("Failed to create TensorRT runtime");
  }

  m_engine.reset(m_runtime->deserializeCudaEngine(engineData.data(), fsize));
  if (!m_engine) {
    throw std::runtime_error("Failed to deserialize CUDA engine");
  }

  for (auto i = 0ul; i < m_cfg.numExecutionContexts; ++i) {
    ACTS_DEBUG("Create execution context " << i);
    m_contexts.emplace_back(m_engine->createExecutionContext());
    if (!m_contexts.back()) {
      throw std::runtime_error("Failed to create execution context");
    }
  }

  std::size_t freeMem{}, totalMem{};
  ACTS_CUDA_CHECK(cudaMemGetInfo(&freeMem, &totalMem));
  ACTS_DEBUG("Used CUDA memory after TensorRT initialization: "
             << (totalMem - freeMem) * 1e-9 << " / " << totalMem * 1e-9
             << " GB");
}

TensorRTEdgeClassifier::~TensorRTEdgeClassifier() {}

PipelineTensors TensorRTEdgeClassifier::operator()(
    PipelineTensors tensors, const ExecutionContext &execContext) {
  assert(execContext.device.type == Acts::Device::Type::eCUDA);

  decltype(std::chrono::high_resolution_clock::now()) t0, t1, t2, t3, t4;
  t0 = std::chrono::high_resolution_clock::now();

  // get a context from the list of contexts
  std::unique_ptr<nvinfer1::IExecutionContext> context;
  while (context == nullptr) {
    std::lock_guard<std::mutex> lock(m_contextMutex);
    if (!m_contexts.empty()) {
      context = std::move(m_contexts.back());
      m_contexts.pop_back();
    }
  }
  assert(context != nullptr);

  context->setInputShape(
      "x", nvinfer1::Dims2{static_cast<long>(tensors.nodeFeatures.shape()[0]),
                           static_cast<long>(tensors.nodeFeatures.shape()[1])});
  context->setTensorAddress("x", tensors.nodeFeatures.data());

  context->setInputShape(
      "edge_index",
      nvinfer1::Dims2{static_cast<long>(tensors.edgeIndex.shape()[0]),
                      static_cast<long>(tensors.edgeIndex.shape()[1])});
  context->setTensorAddress("edge_index", tensors.edgeIndex.data());

  if (tensors.edgeFeatures.has_value()) {
    context->setInputShape(
        "edge_attr",
        nvinfer1::Dims2{static_cast<long>(tensors.edgeFeatures->shape()[0]),
                        static_cast<long>(tensors.edgeFeatures->shape()[1])});
    context->setTensorAddress("edge_attr", tensors.edgeFeatures->data());
  }

  auto scores =
      Tensor<float>::Create({tensors.edgeIndex.shape()[1], 1ul}, execContext);
  context->setTensorAddress("output", scores.data());

  t2 = std::chrono::high_resolution_clock::now();

  auto stream = execContext.stream.value();
  auto status = context->enqueueV3(stream);
  if (!status) {
    throw std::runtime_error("Failed to execute TensorRT model");
  }
  ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));

  t3 = std::chrono::high_resolution_clock::now();

  {
    std::lock_guard<std::mutex> lock(m_contextMutex);
    m_contexts.push_back(std::move(context));
  }

  sigmoid(scores, execContext.stream);

  ACTS_VERBOSE("Size after classifier: " << scores.shape()[0]);
  printCudaMemInfo(logger());

  auto [newScores, newEdgeIndex] =
      applyScoreCut(scores, tensors.edgeIndex, m_cfg.cut, execContext.stream);
  ACTS_VERBOSE("Size after score cut: " << newEdgeIndex.shape()[1]);
  printCudaMemInfo(logger());

  t4 = std::chrono::high_resolution_clock::now();

  auto milliseconds = [](const auto &a, const auto &b) {
    return std::chrono::duration<double, std::milli>(b - a).count();
  };
  ACTS_DEBUG("Time anycast:  " << milliseconds(t0, t1));
  ACTS_DEBUG("Time alloc, set shape " << milliseconds(t1, t2));
  ACTS_DEBUG("Time inference:       " << milliseconds(t2, t3));
  ACTS_DEBUG("Time sigmoid and cut: " << milliseconds(t3, t4));

  return {std::move(tensors.nodeFeatures), std::move(newEdgeIndex),
          std::nullopt, std::move(newScores)};
}

}  // namespace Acts
