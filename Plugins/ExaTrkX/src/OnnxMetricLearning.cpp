// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/OnnxMetricLearning.hpp"

#include <core/session/onnxruntime_cxx_api.h>
#include <torch/script.h>

#include "buildEdges.hpp"
#include "runSessionWithIoBinding.hpp"

namespace Acts {

OnnxMetricLearning::OnnxMetricLearning(Config) {
  m_env = std::make_unique<Ort::Env>(ORT_LOGGING_LEVEL_WARNING,
                                     "ExaTrkX - metric learning");

  Ort::SessionOptions session_options;
  session_options.SetIntraOpNumThreads(1);
  // OrtStatus* status =
  // OrtSessionOptionsAppendExecutionProvider_CUDA(session_options, 0);
  session_options.SetGraphOptimizationLevel(
      GraphOptimizationLevel::ORT_ENABLE_EXTENDED);

  m_model = std::make_unique<Ort::Session>(*m_env, m_cfg.modelPath.c_str(),
                                           session_options);
}

OnnxMetricLearning::~OnnxMetricLearning() {}

void OnnxMetricLearning::buildEdgesWrapper(std::vector<float>& embedFeatures,
                                           std::vector<int64_t>& edgeList,
                                           int64_t numSpacepoints) const {
  torch::Device device(torch::kCUDA);
  auto options =
      torch::TensorOptions().dtype(torch::kFloat32).device(torch::kCUDA);

  torch::Tensor embedTensor =
      torch::tensor(embedFeatures, options)
          .reshape({1, numSpacepoints, m_cfg.embeddingDim});

  auto stackedEdges = buildEdges(embedTensor, numSpacepoints,
                                 m_cfg.embeddingDim, m_cfg.rVal, m_cfg.knnVal);

  stackedEdges = stackedEdges.toType(torch::kInt64).to(torch::kCPU);

  ACTS_INFO("copy edges to std::vector");
  std::copy(stackedEdges.data_ptr<int64_t>(),
            stackedEdges.data_ptr<int64_t>() + stackedEdges.numel(),
            std::back_inserter(edgeList));
}

std::tuple<std::any, std::any> OnnxMetricLearning::operator()(
    std::vector<float>& inputValues) {
  Ort::AllocatorWithDefaultOptions allocator;
  auto memoryInfo = Ort::MemoryInfo::CreateCpu(
      OrtAllocatorType::OrtArenaAllocator, OrtMemType::OrtMemTypeDefault);

  // printout the r,phi,z of the first spacepoint
  ACTS_INFO("First spacepoint information: " << inputValues.size() << "\n\t");
  for (size_t i = 0; i < 3; i++) {
    ACTS_INFO("\t" << inputValues[i]);
  }

  // ************
  // Embedding
  // ************

  int64_t numSpacepoints = inputValues.size() / m_cfg.spacepointFeatures;
  std::vector<int64_t> eInputShape{numSpacepoints, m_cfg.spacepointFeatures};

  std::vector<const char*> eInputNames{"sp_features"};
  auto eInputTensor = std::make_shared<std::vector<Ort::Value>>();
  eInputTensor->push_back(Ort::Value::CreateTensor<float>(
      memoryInfo, inputValues.data(), inputValues.size(), eInputShape.data(),
      eInputShape.size()));

  std::vector<float> eOutputData(numSpacepoints * m_cfg.embeddingDim);
  std::vector<const char*> eOutputNames{"embedding_output"};
  std::vector<int64_t> eOutputShape{numSpacepoints, m_cfg.embeddingDim};
  std::vector<Ort::Value> eOutputTensor;
  eOutputTensor.push_back(Ort::Value::CreateTensor<float>(
      memoryInfo, eOutputData.data(), eOutputData.size(), eOutputShape.data(),
      eOutputShape.size()));
  runSessionWithIoBinding(*m_model, eInputNames, *eInputTensor, eOutputNames,
                          eOutputTensor);

  ACTS_INFO("Embedding space of the first SP: ");
  for (size_t i = 0; i < 3; i++) {
    ACTS_INFO("\t" << eOutputData[i]);
  }
  // if (debug) {
  //   std::fstream out(embedding_outname, out.out);
  //   if (!out.is_open()) {
  //     ACTS_ERROR("failed to open " << embedding_outname);
  //   } else {
  //     std::copy(eOutputData.begin(), eOutputData.end(),
  //               std::ostream_iterator<float>(out, " "));
  //   }
  // }

  // ************
  // Building Edges
  // ************
  std::vector<int64_t> edgeList;
  buildEdgesWrapper(eOutputData, edgeList, numSpacepoints);
  int64_t numEdges = edgeList.size() / 2;
  ACTS_INFO("Built " << numEdges << " edges.");

  for (size_t i = 0; i < 10; i++) {
    ACTS_INFO(edgeList[i]);
  }
  for (size_t i = 0; i < 10; i++) {
    ACTS_INFO(edgeList[numEdges + i]);
  }

  // if (debug) {
  //   std::fstream out(edgelist_outname, out.out);
  //   if (!out.is_open()) {
  //     ACTS_ERROR("failed to open " << edgelist_outname);
  //   } else {
  //     std::copy(edgeList.begin(), edgeList.end(),
  //               std::ostream_iterator<int64_t>(out, " "));
  //   }
  // }

  return {eInputTensor, edgeList};
}

}  // namespace Acts
