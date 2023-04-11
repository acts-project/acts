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
  if (!m_cfg.trackFinderML) {
    throw std::invalid_argument("Missing track finder");
  }

  m_inputSpacePoints.initialize(m_cfg.inputSpacePoints);
  m_outputProtoTracks.initialize(m_cfg.outputProtoTracks);
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

  // ProtoTrackContainer protoTracks;
  std::vector<std::vector<int> > trackCandidates;
  m_cfg.trackFinderML->getTracks(inputValues, spacepointIDs, trackCandidates,
                                 logger());

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
