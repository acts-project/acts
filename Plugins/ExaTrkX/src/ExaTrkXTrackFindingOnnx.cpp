// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/ExaTrkXTrackFindingOnnx.hpp"

#include <filesystem>

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <onnxruntime_cxx_api.h>
#include <torch/script.h>
#include <torch/torch.h>

#include "buildEdges.hpp"
#include "weaklyConnectedComponentsCugraph.hpp"

using namespace torch::indexing;

Acts::ExaTrkXTrackFindingOnnx::ExaTrkXTrackFindingOnnx(
    const Config& config, std::unique_ptr<const Logger> _logger)
    : Acts::ExaTrkXTrackFindingBase("ExaTrkXOnnx"),
      m_cfg(config),
      m_logger{std::move(_logger)} {
  ACTS_INFO("Model input directory: " << m_cfg.modelDir);
  ACTS_INFO("Spacepoint features: " << m_cfg.spacepointFeatures);
  ACTS_INFO("Embedding Dimension: " << m_cfg.embeddingDim);
  ACTS_INFO("radius value       : " << m_cfg.rVal);
  ACTS_INFO("k-nearest neigbour : " << m_cfg.knnVal);
  ACTS_INFO("filtering cut      : " << m_cfg.filterCut);

  m_env = std::make_unique<Ort::Env>(ORT_LOGGING_LEVEL_WARNING, "ExaTrkX");

  using Path = std::filesystem::path;
  const Path embedModelPath = Path(m_cfg.modelDir) / "embedding.onnx";
  const Path filterModelPath = Path(m_cfg.modelDir) / "filtering.onnx";
  const Path gnnModelPath = Path(m_cfg.modelDir) / "gnn.onnx";

  Ort::SessionOptions session_options;
  session_options.SetIntraOpNumThreads(1);
  // OrtStatus* status =
  // OrtSessionOptionsAppendExecutionProvider_CUDA(session_options, 0);
  session_options.SetGraphOptimizationLevel(
      GraphOptimizationLevel::ORT_ENABLE_EXTENDED);

  m_embeddingSession = std::make_unique<Ort::Session>(
      *m_env, embedModelPath.c_str(), session_options);
  m_filterSession = std::make_unique<Ort::Session>(
      *m_env, filterModelPath.c_str(), session_options);
  m_gnnSession = std::make_unique<Ort::Session>(*m_env, gnnModelPath.c_str(),
                                                session_options);
}

Acts::ExaTrkXTrackFindingOnnx::~ExaTrkXTrackFindingOnnx() {}

void Acts::ExaTrkXTrackFindingOnnx::runSessionWithIoBinding(
    Ort::Session& sess, std::vector<const char*>& inputNames,
    std::vector<Ort::Value>& inputData, std::vector<const char*>& outputNames,
    std::vector<Ort::Value>& outputData) const {
  if (inputNames.size() < 1) {
    throw std::runtime_error("Onnxruntime input data maping cannot be empty");
  }
  if (inputNames.size() != inputData.size()) {
    throw std::runtime_error("inputData size mismatch");
  }

  Ort::IoBinding iobinding(sess);
  for (size_t idx = 0; idx < inputNames.size(); ++idx) {
    iobinding.BindInput(inputNames[idx], inputData[idx]);
  }

  for (size_t idx = 0; idx < outputNames.size(); ++idx) {
    iobinding.BindOutput(outputNames[idx], outputData[idx]);
  }

  sess.Run(Ort::RunOptions{nullptr}, iobinding);
}

void Acts::ExaTrkXTrackFindingOnnx::buildEdges(
    std::vector<float>& embedFeatures, std::vector<int64_t>& edgeList,
    int64_t numSpacepoints) const {
  torch::Device device(torch::kCUDA);
  auto options =
      torch::TensorOptions().dtype(torch::kFloat32).device(torch::kCUDA);

  torch::Tensor embedTensor =
      torch::tensor(embedFeatures, options)
          .reshape({1, numSpacepoints, m_cfg.embeddingDim});

  auto stackedEdges =
      Acts::buildEdges(embedTensor, numSpacepoints, m_cfg.embeddingDim,
                       m_cfg.rVal, m_cfg.knnVal);

  stackedEdges = stackedEdges.toType(torch::kInt64).to(torch::kCPU);

  ACTS_INFO("copy edges to std::vector");
  std::copy(stackedEdges.data_ptr<int64_t>(),
            stackedEdges.data_ptr<int64_t>() + stackedEdges.numel(),
            std::back_inserter(edgeList));
}

std::optional<Acts::ExaTrkXTime> Acts::ExaTrkXTrackFindingOnnx::getTracks(
    std::vector<float>& inputValues, std::vector<int>& spacepointIDs,
    std::vector<std::vector<int> >& trackCandidates, const Logger& logger,
    bool recordTiming) const {
  ExaTrkXTime timeInfo;
  // hardcoded debugging information
  bool debug = true;
  const std::string embedding_outname = "debug_embedding_outputs.txt";
  const std::string edgelist_outname = "debug_edgelist_outputs.txt";
  const std::string filtering_outname = "debug_filtering_scores.txt";

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
  std::vector<Ort::Value> eInputTensor;
  eInputTensor.push_back(Ort::Value::CreateTensor<float>(
      memoryInfo, inputValues.data(), inputValues.size(), eInputShape.data(),
      eInputShape.size()));

  std::vector<float> eOutputData(numSpacepoints * m_cfg.embeddingDim);
  std::vector<const char*> eOutputNames{"embedding_output"};
  std::vector<int64_t> eOutputShape{numSpacepoints, m_cfg.embeddingDim};
  std::vector<Ort::Value> eOutputTensor;
  eOutputTensor.push_back(Ort::Value::CreateTensor<float>(
      memoryInfo, eOutputData.data(), eOutputData.size(), eOutputShape.data(),
      eOutputShape.size()));
  runSessionWithIoBinding(*m_embeddingSession, eInputNames, eInputTensor,
                          eOutputNames, eOutputTensor);

  ACTS_INFO("Embedding space of the first SP: ");
  for (size_t i = 0; i < 3; i++) {
    ACTS_INFO("\t" << eOutputData[i]);
  }
  if (debug) {
    std::fstream out(embedding_outname, out.out);
    if (!out.is_open()) {
      ACTS_ERROR("failed to open " << embedding_outname);
    } else {
      std::copy(eOutputData.begin(), eOutputData.end(),
                std::ostream_iterator<float>(out, " "));
    }
  }

  // ************
  // Building Edges
  // ************
  std::vector<int64_t> edgeList;
  buildEdges(eOutputData, edgeList, numSpacepoints);
  int64_t numEdges = edgeList.size() / 2;
  ACTS_INFO("Built " << numEdges << " edges.");

  for (size_t i = 0; i < 10; i++) {
    ACTS_INFO(edgeList[i]);
  }
  for (size_t i = 0; i < 10; i++) {
    ACTS_INFO(edgeList[numEdges + i]);
  }

  if (debug) {
    std::fstream out(edgelist_outname, out.out);
    if (!out.is_open()) {
      ACTS_ERROR("failed to open " << edgelist_outname);
    } else {
      std::copy(edgeList.begin(), edgeList.end(),
                std::ostream_iterator<int64_t>(out, " "));
    }
  }

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

  // ************
  // GNN
  // ************
  std::vector<const char*> gInputNames{"g_nodes", "g_edges"};
  std::vector<Ort::Value> gInputTensor;
  gInputTensor.push_back(std::move(fInputTensor[0]));
  std::vector<int64_t> gEdgeShape{2, numEdgesAfterF};
  gInputTensor.push_back(Ort::Value::CreateTensor<int64_t>(
      memoryInfo, edgesAfterFiltering.data(), edgesAfterFiltering.size(),
      gEdgeShape.data(), gEdgeShape.size()));
  // gnn outputs
  std::vector<const char*> gOutputNames{"gnn_edge_score"};
  std::vector<float> gOutputData(numEdgesAfterF);
  std::vector<int64_t> gOutputShape{numEdgesAfterF};
  std::vector<Ort::Value> gOutputTensor;
  gOutputTensor.push_back(Ort::Value::CreateTensor<float>(
      memoryInfo, gOutputData.data(), gOutputData.size(), gOutputShape.data(),
      gOutputShape.size()));

  ACTS_INFO("run ONNX session");
  runSessionWithIoBinding(*m_gnnSession, gInputNames, gInputTensor,
                          gOutputNames, gOutputTensor);
  ACTS_INFO("done with ONNX session");

  torch::Tensor gOutputCTen = torch::tensor(gOutputData, {torch::kFloat32});
  gOutputCTen = gOutputCTen.sigmoid();
  ACTS_INFO(gOutputCTen.slice(0, 0, 3));

  // ************
  // Track Labeling with cugraph::connected_components
  // ************
  std::vector<int32_t> rowIndices;
  std::vector<int32_t> colIndices;
  std::vector<float> edgeWeights;
  std::vector<int32_t> trackLabels(numSpacepoints);
  std::copy(edgesAfterFiltering.begin(),
            edgesAfterFiltering.begin() + numEdgesAfterF,
            std::back_insert_iterator(rowIndices));
  std::copy(edgesAfterFiltering.begin() + numEdgesAfterF,
            edgesAfterFiltering.end(), std::back_insert_iterator(colIndices));
  std::copy(gOutputCTen.data_ptr<float>(),
            gOutputCTen.data_ptr<float>() + numEdgesAfterF,
            std::back_insert_iterator(edgeWeights));

  ACTS_INFO("run weaklyConnectedComponents");
  weaklyConnectedComponents<int32_t, int32_t, float>(rowIndices, colIndices,
                                                     edgeWeights, trackLabels);

  ACTS_INFO("size of components: " << trackLabels.size());
  if (trackLabels.size() == 0) {
    return timeInfo;
  }

  trackCandidates.clear();

  int existTrkIdx = 0;
  // map labeling from MCC to customized track id.
  std::map<int, int> trackLableToIds;

  for (int idx = 0; idx < numSpacepoints; ++idx) {
    int trackLabel = trackLabels[idx];
    int spacepointID = spacepointIDs[idx];

    int trkId;
    if (trackLableToIds.find(trackLabel) != trackLableToIds.end()) {
      trkId = trackLableToIds[trackLabel];
      trackCandidates[trkId].push_back(spacepointID);
    } else {
      // a new track, assign the track id
      // and create a vector
      trkId = existTrkIdx;
      trackCandidates.push_back(std::vector<int>{trkId});
      trackLableToIds[trackLabel] = trkId;
      existTrkIdx++;
    }
  }

  return timeInfo;
}
