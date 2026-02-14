// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFindingML/AmbiguityResolutionMLAlgorithm.hpp"

#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <iterator>
#include <map>

using namespace Acts;
using namespace ActsPlugins;

namespace ActsExamples {

namespace {

static std::size_t sourceLinkHash(const SourceLink& a) {
  return static_cast<std::size_t>(a.get<IndexSourceLink>().index());
}

static bool sourceLinkEquality(const SourceLink& a, const SourceLink& b) {
  return a.get<IndexSourceLink>().index() == b.get<IndexSourceLink>().index();
}

}  // namespace

AmbiguityResolutionMLAlgorithm::AmbiguityResolutionMLAlgorithm(
    const Config& cfg, Logging::Level lvl)
    : IAlgorithm("AmbiguityResolutionMLAlgorithm", lvl),
      m_cfg(cfg),
      m_ambiML(m_cfg.toAmbiguityResolutionMLConfig(), logger().clone()) {
  if (m_cfg.inputTracks.empty()) {
    throw std::invalid_argument("Missing trajectories input collection");
  }
  if (m_cfg.outputTracks.empty()) {
    throw std::invalid_argument("Missing trajectories output collection");
  }
  m_inputTracks.initialize(m_cfg.inputTracks);
  m_outputTracks.initialize(m_cfg.outputTracks);
}

ProcessCode AmbiguityResolutionMLAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  // Read input data
  const auto& tracks = m_inputTracks(ctx);
  // Associate measurement to their respective tracks to prepare the track
  // shared hits based clustering
  std::multimap<int, std::pair<std::size_t, std::vector<std::size_t>>>
      trackMap =
          m_ambiML.mapTrackHits(tracks, &sourceLinkHash, &sourceLinkEquality);
  // Cluster the tracks based on the shared hits
  auto cluster = Acts::detail::clusterDuplicateTracks(trackMap);
  // Select the ID of the track we want to keep
  std::vector<std::size_t> goodTracks =
      m_ambiML.solveAmbiguity(cluster, tracks);
  // Prepare the output track collection from the IDs
  TrackContainer solvedTracks{std::make_shared<VectorTrackContainer>(),
                              std::make_shared<VectorMultiTrajectory>()};
  solvedTracks.ensureDynamicColumns(tracks);
  for (auto iTrack : goodTracks) {
    auto destProxy = solvedTracks.makeTrack();
    auto srcProxy = tracks.getTrack(iTrack);
    destProxy.copyFromWithoutStates(srcProxy);
    destProxy.tipIndex() = srcProxy.tipIndex();
  }

  ConstTrackContainer outputTracks{std::make_shared<ConstVectorTrackContainer>(
                                       std::move(solvedTracks.container())),
                                   tracks.trackStateContainerHolder()};

  m_outputTracks(ctx, std::move(outputTracks));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
