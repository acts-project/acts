// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Gnn/TensorRTEdgeClassifier.hpp"

#include "ActsPlugins/Gnn/detail/CudaUtils.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>

#include <NvInfer.h>
#include <NvInferPlugin.h>
#include <NvInferRuntimeBase.h>
#include <cuda_runtime.h>

#include "printCudaMemInfo.hpp"

using namespace Acts;

namespace {

class TensorRTLogger : public nvinfer1::ILogger {
  std::unique_ptr<const Logger> m_logger;

 public:
  TensorRTLogger(Logging::Level lvl)
      : m_logger(getDefaultLogger("TensorRT", lvl)) {}

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

namespace ActsPlugins {

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

  ACTS_INFO("Device memory required by TRT context: "
            << m_engine->getDeviceMemorySizeV2() * 1e-9 << " GB");

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

  m_count.emplace(m_contexts.size());

  if (m_engine->getNbOptimizationProfiles() > 1) {
    ACTS_WARNING("Cannot handle more then one optimization profile for now");
  }

  m_maxNodes =
      m_engine->getProfileShape("x", 0, nvinfer1::OptProfileSelector::kMAX)
          .d[0];
  ACTS_INFO("Maximum number of nodes: " << m_maxNodes);

  auto maxEdgesA =
      m_engine
          ->getProfileShape("edge_index", 0, nvinfer1::OptProfileSelector::kMAX)
          .d[1];
  auto maxEdgesB =
      m_engine
          ->getProfileShape("edge_attr", 0, nvinfer1::OptProfileSelector::kMAX)
          .d[0];

  if (maxEdgesA != maxEdgesB) {
    throw std::invalid_argument(
        "Inconsistent max edges definition in engine for 'edge_index' and "
        "'edge_attr'");
  }

  m_maxEdges = maxEdgesA;
  ACTS_INFO("Maximum number of edges: " << m_maxEdges);
}

TensorRTEdgeClassifier::~TensorRTEdgeClassifier() {}

PipelineTensors TensorRTEdgeClassifier::operator()(
    PipelineTensors tensors, const ExecutionContext &execContext) {
  assert(execContext.device.type == Device::Type::eCUDA);

  decltype(std::chrono::high_resolution_clock::now()) t0, t1, t2, t3, t4;
  t0 = std::chrono::high_resolution_clock::now();

  // Curing this would require more complicated handling, and should happen
  // almost never
  if (auto nNodes = tensors.nodeFeatures.shape().at(0); nNodes > m_maxNodes) {
    ACTS_WARNING("Number of nodes ("
                 << nNodes << ") exceeds configured maximum, return 0 edges");
    throw NoEdgesError{};
  }

  if (auto nEdges = tensors.edgeIndex.shape().at(1); nEdges > m_maxEdges) {
    ACTS_WARNING("Number of edges ("
                 << nEdges << ") exceeds maximum, shrink edge tensor to "
                 << m_maxEdges);

    auto [newEdgeIndex, newEdgeFeatures] =
        applyEdgeLimit(tensors.edgeIndex, tensors.edgeFeatures, m_maxEdges,
                       execContext.stream);
    tensors.edgeIndex = std::move(newEdgeIndex);
    tensors.edgeFeatures = std::move(newEdgeFeatures);
  }

  // get a context from the list of contexts
  std::unique_ptr<nvinfer1::IExecutionContext> context;

  // Try to decrement the count before getting a context. By this it should be
  // garantueed that a context is available
  m_count->acquire();
  assert(!m_contexts.empty());

  // Protect access to the context by a mutex
  {
    std::lock_guard<std::mutex> lock(m_contextMutex);
    context = std::move(m_contexts.back());
    m_contexts.pop_back();
  }

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

  // Increment the count after the context have been put back in the vector
  m_count->release();

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

}  // namespace ActsPlugins
