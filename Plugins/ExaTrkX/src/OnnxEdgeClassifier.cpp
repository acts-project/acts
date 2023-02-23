// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/OnnxEdgeClassifier.hpp"

#include <core/session/onnxruntime_cxx_api.h>

namespace Acts {
    
OnnxEdgeClassifier::OnnxEdgeClassifier(Config cfg, const Logger &logger) : EdgeClassificationBase(logger), m_cfg(cfg) {
  m_env = std::make_unique<Ort::Env>(ORT_LOGGING_LEVEL_WARNING,
                                     "ExaTrkX - edge classifier");

  Ort::SessionOptions session_options;
  session_options.SetIntraOpNumThreads(1);
  // OrtStatus* status =
  // OrtSessionOptionsAppendExecutionProvider_CUDA(session_options, 0);
  session_options.SetGraphOptimizationLevel(
      GraphOptimizationLevel::ORT_ENABLE_EXTENDED);

  m_model = std::make_unique<Ort::Session>(*m_env, m_cfg.modelPath.c_str(),
                                           session_options);
}

OnnxEdgeClassifier::~OnnxEdgeClassifier() {}

std::tuple<std::any, std::any, std::any> OnnxEdgeClassifier::operator()(std::any inputNodes,
                                            std::any inputEdges) {
    
  // ************
  // Filtering
  // ************
  std::vector<const char*> fInputNames{"f_nodes", "f_edges"};
  std::vector<Ort::Value> fInputTensor;
  fInputTensor.push_back(std::move(eInputTensor[0]));
  std::vector<int64_t> fEdgeShape{2, numEdges};
  fInputTensor.push_back(Ort::Value::CreateTensor<int64_t>(
      memoryInfo, edgeList.data(), edgeList.size(), fEdgeShape.data(),
      fEdgeShape.size()));

  // filtering outputs
  std::vector<const char*> fOutputNames{"f_edge_score"};
  std::vector<float> fOutputData(numEdges);
  std::vector<int64_t> fOutputShape{numEdges, 1};
  std::vector<Ort::Value> fOutputTensor;
  fOutputTensor.push_back(Ort::Value::CreateTensor<float>(
      memoryInfo, fOutputData.data(), fOutputData.size(), fOutputShape.data(),
      fOutputShape.size()));
  runSessionWithIoBinding(*m_filterSession, fInputNames, fInputTensor,
                          fOutputNames, fOutputTensor);

  ACTS_INFO("Get scores for " << numEdges << " edges.");
  // However, I have to convert those numbers to a score by applying sigmoid!
  // Use torch::tensor
  torch::Tensor edgeListCTen = torch::tensor(edgeList, {torch::kInt64});
  edgeListCTen = edgeListCTen.reshape({2, numEdges});

  torch::Tensor fOutputCTen = torch::tensor(fOutputData, {torch::kFloat32});
  fOutputCTen = fOutputCTen.sigmoid();

  if (debug) {
    std::fstream out(filtering_outname, out.out);
    if (!out.is_open()) {
      ACTS_ERROR("failed to open " << filtering_outname);
    } else {
      std::copy(fOutputCTen.data_ptr<float>(),
                fOutputCTen.data_ptr<float>() + fOutputCTen.numel(),
                std::ostream_iterator<float>(out, " "));
    }
  }
  // std::cout << fOutputCTen.slice(0, 0, 3) << std::endl;
  torch::Tensor filterMask = fOutputCTen > m_cfg.filterCut;
  torch::Tensor edgesAfterFCTen = edgeListCTen.index({Slice(), filterMask});

  std::vector<int64_t> edgesAfterFiltering;
  std::copy(edgesAfterFCTen.data_ptr<int64_t>(),
            edgesAfterFCTen.data_ptr<int64_t>() + edgesAfterFCTen.numel(),
            std::back_inserter(edgesAfterFiltering));

  int64_t numEdgesAfterF = edgesAfterFiltering.size() / 2;
  ACTS_INFO("After filtering: " << numEdgesAfterF << " edges.");
}

}
