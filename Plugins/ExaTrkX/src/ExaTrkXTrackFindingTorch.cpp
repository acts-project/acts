// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/ExaTrkXTrackFindingTorch.hpp"

#include <filesystem>

#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <grid/counting_sort.h>
#include <grid/find_nbrs.h>
#include <grid/grid.h>
#include <grid/insert_points.h>
#include <grid/prefix_sum.h>
#include <torch/script.h>
#include <torch/torch.h>

#include "buildEdges.hpp"
#include "weaklyConnectedComponentsBoost.hpp"

using namespace torch::indexing;

namespace {
void print_current_cuda_meminfo(Acts::LoggerWrapper& logger) {
  constexpr int kb = 1024;
  constexpr int mb = kb * kb;

  int device;
  std::size_t free, total;
  cudaMemGetInfo(&free, &total);
  cudaGetDevice(&device);

  ACTS_VERBOSE("Current CUDA device: " << device);
  ACTS_VERBOSE("Memory (used / total) [in MB]: " << (total - free) / mb << " / "
                                                 << total / mb);
}
}  // namespace

namespace Acts {

ExaTrkXTrackFindingTorch::ExaTrkXTrackFindingTorch(
    const ExaTrkXTrackFindingTorch::Config& config)
    : ExaTrkXTrackFindingBase("ExaTrkXTorch"), m_cfg(config) {
  using Path = std::filesystem::path;

  const Path embedModelPath = Path(m_cfg.modelDir) / "embed.pt";
  const Path filterModelPath = Path(m_cfg.modelDir) / "filter.pt";
  const Path gnnModelPath = Path(m_cfg.modelDir) / "gnn.pt";
  c10::InferenceMode guard(true);

  try {
    m_embeddingModel = std::make_unique<torch::jit::Module>();
    *m_embeddingModel = torch::jit::load(embedModelPath.c_str());
    m_embeddingModel->eval();

    m_filterModel = std::make_unique<torch::jit::Module>();
    *m_filterModel = torch::jit::load(filterModelPath.c_str());
    m_filterModel->eval();

    m_gnnModel = std::make_unique<torch::jit::Module>();
    *m_gnnModel = torch::jit::load(gnnModelPath.c_str());
    m_gnnModel->eval();
  } catch (const c10::Error& e) {
    throw std::invalid_argument("Failed to load models: " + e.msg());
  }
}

ExaTrkXTrackFindingTorch::~ExaTrkXTrackFindingTorch() {}

std::optional<ExaTrkXTime> ExaTrkXTrackFindingTorch::getTracks(
    std::vector<float>& inputValues, std::vector<int>& spacepointIDs,
    std::vector<std::vector<int> >& trackCandidates, LoggerWrapper logger,
    bool recordTiming) const {
  ExaTrkXTime timeInfo;

  ExaTrkXTimer totalTimer(not recordTiming);
  totalTimer.start();

  c10::InferenceMode guard(true);
  torch::Device device(torch::kCUDA);

  // Clone models (solve memory leak? members can be const...)
  auto e_model = m_embeddingModel->clone();
  auto f_model = m_filterModel->clone();
  auto g_model = m_gnnModel->clone();

  // printout the r,phi,z of the first spacepoint
  ACTS_VERBOSE("First spacepoint information [r, phi, z]: "
               << inputValues[0] << ", " << inputValues[1] << ", "
               << inputValues[2]);
  ACTS_VERBOSE("Max and min spacepoint: "
               << *std::max_element(inputValues.begin(), inputValues.end())
               << ", "
               << *std::min_element(inputValues.begin(), inputValues.end()))
  print_current_cuda_meminfo(logger);

  ExaTrkXTimer timer(not recordTiming);

  // **********
  // Embedding
  // **********

  timer.start();
  int64_t numSpacepoints = inputValues.size() / m_cfg.spacepointFeatures;
  std::vector<torch::jit::IValue> eInputTensorJit;
  auto e_opts = torch::TensorOptions().dtype(torch::kFloat32);
  torch::Tensor eLibInputTensor =
      torch::from_blob(inputValues.data(),
                       {numSpacepoints, m_cfg.spacepointFeatures}, e_opts)
          .to(torch::kFloat32);

  eInputTensorJit.push_back(eLibInputTensor.to(device));
  std::optional<at::Tensor> eOutput =
      e_model.forward(eInputTensorJit).toTensor();
  eInputTensorJit.clear();

  ACTS_VERBOSE("Embedding space of the first SP:\n"
               << eOutput->slice(/*dim=*/0, /*start=*/0, /*end=*/1));
  print_current_cuda_meminfo(logger);

  timeInfo.embedding = timer.stopAndGetElapsedTime();

  // ****************
  // Building Edges
  // ****************

  timer.start();

  // At this point, buildEdgesBruteForce could be used instead
  std::optional<torch::Tensor> edgeList = buildEdges(
      *eOutput, numSpacepoints, m_cfg.embeddingDim, m_cfg.rVal, m_cfg.knnVal);
  eOutput.reset();

  ACTS_VERBOSE("Shape of built edges: (" << edgeList->size(0) << ", "
                                         << edgeList->size(1));
  ACTS_VERBOSE("Slice of edgelist:\n" << edgeList->slice(1, 0, 5));
  print_current_cuda_meminfo(logger);

  timeInfo.building = timer.stopAndGetElapsedTime();

  // **********
  // Filtering
  // **********

  timer.start();

  const auto chunks = at::chunk(at::arange(edgeList->size(1)), m_cfg.n_chunks);
  std::vector<at::Tensor> results;

  for (const auto& chunk : chunks) {
    std::vector<torch::jit::IValue> fInputTensorJit;
    fInputTensorJit.push_back(eLibInputTensor.to(device));
    fInputTensorJit.push_back(edgeList->index({Slice(), chunk}).to(device));

    results.push_back(f_model.forward(fInputTensorJit).toTensor());
    results.back().squeeze_();
    results.back().sigmoid_();
  }

  auto fOutput = torch::cat(results);
  results.clear();

  ACTS_VERBOSE("Size after filtering network: " << fOutput.size(0));
  ACTS_VERBOSE("Slice of filtered output:\n"
               << fOutput.slice(/*dim=*/0, /*start=*/0, /*end=*/9));
  print_current_cuda_meminfo(logger);

  torch::Tensor filterMask = fOutput > m_cfg.filterCut;
  torch::Tensor edgesAfterF = edgeList->index({Slice(), filterMask});
  edgeList.reset();
  edgesAfterF = edgesAfterF.to(torch::kInt64);
  const int64_t numEdgesAfterF = edgesAfterF.size(1);

  ACTS_VERBOSE("Size after filter cut: " << numEdgesAfterF)
  print_current_cuda_meminfo(logger);

  timeInfo.filtering = timer.stopAndGetElapsedTime();

  // ****
  // GNN
  // ****

  timer.start();

  auto bidirEdgesAfterF = torch::cat({edgesAfterF, edgesAfterF.flip(0)}, 1);

  ACTS_VERBOSE("Bi-directional edges shape: ("
               << bidirEdgesAfterF.size(0) << ", " << bidirEdgesAfterF.size(1)
               << ")")
  print_current_cuda_meminfo(logger);

  std::vector<torch::jit::IValue> gInputTensorJit;
  gInputTensorJit.push_back(eLibInputTensor.to(device));
  gInputTensorJit.push_back(bidirEdgesAfterF.to(device));

  auto gOutputBidir = g_model.forward(gInputTensorJit).toTensor();
  gInputTensorJit.clear();
  gOutputBidir.sigmoid_();
  gOutputBidir = gOutputBidir.cpu();

  auto gOutput = gOutputBidir.index({Slice(None, gOutputBidir.size(0) / 2)});

  timeInfo.gnn = timer.stopAndGetElapsedTime();

  ACTS_VERBOSE("GNN scores size: " << gOutput.size(0) << " (bidir: "
                                   << gOutputBidir.size(0) << ")");
  ACTS_VERBOSE("Score output slice:\n" << gOutput.slice(0, 0, 5));
  print_current_cuda_meminfo(logger);

  // ***************
  // Track Labeling
  // ***************

  timer.start();

  using vertex_t = int32_t;
  std::vector<vertex_t> rowIndices;
  std::vector<vertex_t> colIndices;
  std::vector<float> edgeWeights;
  std::vector<vertex_t> trackLabels(numSpacepoints);
  std::copy(edgesAfterF.data_ptr<int64_t>(),
            edgesAfterF.data_ptr<int64_t>() + numEdgesAfterF,
            std::back_insert_iterator(rowIndices));
  std::copy(edgesAfterF.data_ptr<int64_t>() + numEdgesAfterF,
            edgesAfterF.data_ptr<int64_t>() + numEdgesAfterF + numEdgesAfterF,
            std::back_insert_iterator(colIndices));
  std::copy(gOutput.data_ptr<float>(),
            gOutput.data_ptr<float>() + numEdgesAfterF,
            std::back_insert_iterator(edgeWeights));

  weaklyConnectedComponents<int32_t, int32_t, float>(
      numSpacepoints, rowIndices, colIndices, edgeWeights, trackLabels,
      m_cfg.edgeCut);

  ACTS_VERBOSE("Number of track labels: " << trackLabels.size());
  ACTS_VERBOSE("NUmber of unique track labels: " << [&]() {
    std::vector<vertex_t> sorted(trackLabels);
    std::sort(sorted.begin(), sorted.end());
    sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());
    return sorted.size();
  }());
  print_current_cuda_meminfo(logger);

  if (trackLabels.size() == 0) {
    if (recordTiming) {
      return timeInfo;
    } else {
      return std::nullopt;
    }
  }

  trackCandidates.clear();

  int existTrkIdx = 0;
  // map labeling from MCC to customized track id.
  std::map<int32_t, int32_t> trackLableToIds;

  for (int32_t idx = 0; idx < numSpacepoints; ++idx) {
    int32_t trackLabel = trackLabels[idx];
    int spacepointID = spacepointIDs[idx];

    int trkId;
    if (trackLableToIds.find(trackLabel) != trackLableToIds.end()) {
      trkId = trackLableToIds[trackLabel];
      trackCandidates[trkId].push_back(spacepointID);
    } else {
      // a new track, assign the track id
      // and create a vector
      trkId = existTrkIdx;
      trackCandidates.push_back(std::vector<int>{spacepointID});
      trackLableToIds[trackLabel] = trkId;
      existTrkIdx++;
    }
  }

  timeInfo.labeling = timer.stopAndGetElapsedTime();
  timeInfo.total = totalTimer.stopAndGetElapsedTime();
  c10::cuda::CUDACachingAllocator::emptyCache();

  if (recordTiming) {
    return timeInfo;
  } else {
    return std::nullopt;
  }
}

}  // namespace Acts
