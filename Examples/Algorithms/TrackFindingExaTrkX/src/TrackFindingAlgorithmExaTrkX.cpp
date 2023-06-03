// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFindingExaTrkX/TrackFindingAlgorithmExaTrkX.hpp"

#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

ActsExamples::TrackFindingAlgorithmExaTrkX::TrackFindingAlgorithmExaTrkX(
    Config config, Acts::Logging::Level level)
    : ActsExamples::IAlgorithm("TrackFindingMLBasedAlgorithm", level),
      m_cfg(std::move(config)) {
  if (m_cfg.inputSpacePoints.empty()) {
    throw std::invalid_argument("Missing spacepoint input collection");
  }
  if (m_cfg.outputProtoTracks.empty()) {
    throw std::invalid_argument("Missing protoTrack output collection");
  }
  if (!m_cfg.graphConstructor) {
    throw std::invalid_argument("Missing graph construction module");
  }
  if (!m_cfg.trackBuilder) {
    throw std::invalid_argument("Missing track building module");
  }
  if (m_cfg.edgeClassifiers.empty() or
      not std::all_of(m_cfg.edgeClassifiers.begin(),
                      m_cfg.edgeClassifiers.end(),
                      [](const auto& a) { return static_cast<bool>(a); })) {
    throw std::invalid_argument("Missing graph construction module");
  }

  // Sanitizer run with dummy input to detect configuration issues
  // TODO This would be quite helpful I think, but currently it does not work in
  // general because the stages do not expose the number of node features.
  // However, this must be addressed anyways when we also want to allow to
  // configure this more flexible with e.g. cluster information as input. So for
  // now, we disable this.
#if 0
  if( m_cfg.sanitize ) {
  Eigen::VectorXf dummyInput = Eigen::VectorXf::Random(3 * 15);
  std::vector<float> dummyInputVec(dummyInput.data(),
                                   dummyInput.data() + dummyInput.size());
  std::vector<int> spacepointIDs;
  std::iota(spacepointIDs.begin(), spacepointIDs.end(), 0);
  
  runPipeline(dummyInputVec, spacepointIDs);
  }
#endif

  m_inputSpacePoints.initialize(m_cfg.inputSpacePoints);
  m_outputProtoTracks.initialize(m_cfg.outputProtoTracks);
}

std::vector<std::vector<int>>
ActsExamples::TrackFindingAlgorithmExaTrkX::runPipeline(
    std::vector<float>& inputValues, std::vector<int>& spacepointIDs) const {
  auto [nodes, edges] = (*m_cfg.graphConstructor)(inputValues, logger());
  std::any edge_weights;

  for (auto edgeClassifier : m_cfg.edgeClassifiers) {
    auto [newNodes, newEdges, newWeights] =
        (*edgeClassifier)(nodes, edges, logger());
    nodes = newNodes;
    edges = newEdges;
    edge_weights = newWeights;
  }

  return (*m_cfg.trackBuilder)(nodes, edges, edge_weights, spacepointIDs,
                               logger());
}

ActsExamples::ProcessCode ActsExamples::TrackFindingAlgorithmExaTrkX::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // Read input data
  const auto& spacepoints = m_inputSpacePoints(ctx);

  // Convert Input data to a list of size [num_measurements x
  // measurement_features]
  size_t num_spacepoints = spacepoints.size();
  ACTS_INFO("Received " << num_spacepoints << " spacepoints");

  std::vector<float> inputValues;
  std::vector<int> spacepointIDs;
  inputValues.reserve(spacepoints.size() * 3);
  spacepointIDs.reserve(spacepoints.size());
  for (const auto& sp : spacepoints) {
    float x = sp.x();
    float y = sp.y();
    float z = sp.z();
    float r = sp.r();
    float phi = std::atan2(y, x);

    inputValues.push_back(r / m_cfg.rScale);
    inputValues.push_back(phi / m_cfg.phiScale);
    inputValues.push_back(z / m_cfg.zScale);

    // For now just take the first index since does require one single index per
    // spacepoint
    const auto& islink = sp.sourceLinks()[0].template get<IndexSourceLink>();
    spacepointIDs.push_back(islink.index());
  }

  // Run the pipeline
  const auto trackCandidates = runPipeline(inputValues, spacepointIDs);

  // Make the prototracks
  std::vector<ProtoTrack> protoTracks;
  protoTracks.reserve(trackCandidates.size());
  for (auto& x : trackCandidates) {
    ProtoTrack onetrack;
    std::copy(x.begin(), x.end(), std::back_inserter(onetrack));
    protoTracks.push_back(std::move(onetrack));
  }

  ACTS_INFO("Created " << protoTracks.size() << " proto tracks");
  m_outputProtoTracks(ctx, std::move(protoTracks));

  return ActsExamples::ProcessCode::SUCCESS;
}
