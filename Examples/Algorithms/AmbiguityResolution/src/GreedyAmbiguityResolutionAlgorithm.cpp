// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/AmbiguityResolution/GreedyAmbiguityResolutionAlgorithm.hpp"

#include "Acts/AmbiguityResolution/GreedyAmbiguityResolution.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <cstddef>
#include <stdexcept>

#include <boost/container/flat_map.hpp>
#include <boost/container/flat_set.hpp>

namespace ActsExamples {

namespace {

Acts::GreedyAmbiguityResolution::Config transformConfig(
    const GreedyAmbiguityResolutionAlgorithm::Config& cfg) {
  Acts::GreedyAmbiguityResolution::Config result;
  result.maximumSharedHits = cfg.maximumSharedHits;
  result.maximumIterations = cfg.maximumIterations;
  result.nMeasurementsMin = cfg.nMeasurementsMin;
  return result;
}

std::size_t sourceLinkHash(const Acts::SourceLink& a) {
  return static_cast<std::size_t>(a.get<IndexSourceLink>().index());
}

bool sourceLinkEquality(const Acts::SourceLink& a, const Acts::SourceLink& b) {
  return a.get<IndexSourceLink>().index() == b.get<IndexSourceLink>().index();
}

}  // namespace

GreedyAmbiguityResolutionAlgorithm::GreedyAmbiguityResolutionAlgorithm(
    const Config& cfg, Acts::Logging::Level lvl)
    : IAlgorithm("GreedyAmbiguityResolutionAlgorithm", lvl),
      m_cfg(cfg),
      m_core(transformConfig(cfg), logger().clone()) {
  if (m_cfg.inputTracks.empty()) {
    throw std::invalid_argument("Missing trajectories input collection");
  }
  if (m_cfg.outputTracks.empty()) {
    throw std::invalid_argument("Missing trajectories output collection");
  }
  m_inputTracks.initialize(m_cfg.inputTracks);
  m_outputTracks.initialize(m_cfg.outputTracks);
}

ProcessCode GreedyAmbiguityResolutionAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  const auto& tracks = m_inputTracks(ctx);

  ACTS_VERBOSE("Number of input tracks: " << tracks.size());

  Acts::GreedyAmbiguityResolution::State state;
  m_core.computeInitialState(tracks, state, &sourceLinkHash,
                             &sourceLinkEquality);

  ACTS_VERBOSE("State initialized");

  m_core.resolve(state);

  ACTS_DEBUG("Resolved to " << state.selectedTracks.size() << " tracks from "
                            << tracks.size());

  TrackContainer solvedTracks{std::make_shared<Acts::VectorTrackContainer>(),
                              std::make_shared<Acts::VectorMultiTrajectory>()};
  solvedTracks.ensureDynamicColumns(tracks);

  for (auto iTrack : state.selectedTracks) {
    auto destProxy = solvedTracks.makeTrack();
    auto srcProxy = tracks.getTrack(state.trackTips.at(iTrack));
    destProxy.copyFromWithoutStates(srcProxy);
    destProxy.tipIndex() = srcProxy.tipIndex();
  }

  ConstTrackContainer outputTracks{
      std::make_shared<Acts::ConstVectorTrackContainer>(
          std::move(solvedTracks.container())),
      tracks.trackStateContainerHolder()};

  m_outputTracks(ctx, std::move(outputTracks));
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
