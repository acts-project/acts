// This file is part of the Acts project.
//
// Copyright (C) 2019-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TruthTracking/TrackSelector.hpp"

#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Utilities/ThrowAssert.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

ActsExamples::TrackSelector::TrackSelector(const Config& config,
                                           Acts::Logging::Level level)
    : IAlgorithm("TrackSelector", level), m_cfg(config) {
  if (m_cfg.inputTracks.empty()) {
    throw std::invalid_argument("Input track collection is empty");
  }

  if (m_cfg.outputTracks.empty()) {
    throw std::invalid_argument("Output track collection is empty");
  }

  m_inputTrackContainer.initialize(m_cfg.inputTracks);
  m_outputTrackContainer.initialize(m_cfg.outputTracks);
}

ActsExamples::ProcessCode ActsExamples::TrackSelector::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // helper functions to select tracks
  auto within = [](double x, double min, double max) {
    return (min <= x) and (x < max);
  };
  auto isValidTrack = [&](const auto& trk) {
    const auto theta = trk.theta();
    const auto eta = -std::log(std::tan(theta / 2));
    return within(trk.transverseMomentum(), m_cfg.ptMin, m_cfg.ptMax) and
           within(std::abs(eta), m_cfg.absEtaMin, m_cfg.absEtaMax) and
           within(eta, m_cfg.etaMin, m_cfg.etaMax) and
           within(trk.phi(), m_cfg.phiMin, m_cfg.phiMax) and
           within(trk.loc0(), m_cfg.loc0Min, m_cfg.loc0Max) and
           within(trk.loc1(), m_cfg.loc1Min, m_cfg.loc1Max) and
           within(trk.time(), m_cfg.timeMin, m_cfg.timeMax);
  };

  ACTS_VERBOSE("Reading tracks from: " << m_cfg.inputTracks);

  const auto& inputTracks = m_inputTrackContainer(ctx);

  std::shared_ptr<Acts::ConstVectorMultiTrajectory> trackStateContainer =
      inputTracks.trackStateContainerHolder();

  auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();

  // temporary empty track state container: we don't change the original one,
  // but we need one for filtering
  auto tempTrackStateContainer =
      std::make_shared<Acts::VectorMultiTrajectory>();

  TrackContainer filteredTracks{trackContainer, tempTrackStateContainer};
  filteredTracks.ensureDynamicColumns(inputTracks);

  trackContainer->reserve(inputTracks.size());

  ACTS_VERBOSE("Track container size before filtering: " << inputTracks.size());

  for (auto track : inputTracks) {
    if (!isValidTrack(track)) {
      continue;
    }
    auto destProxy = filteredTracks.getTrack(filteredTracks.addTrack());
    destProxy.copyFrom(track);
  }

  ACTS_VERBOSE(
      "Track container size after filtering: " << filteredTracks.size());

  ConstTrackContainer outputTracks{
      std::make_shared<Acts::ConstVectorTrackContainer>(
          std::move(*trackContainer)),
      trackStateContainer};

  m_outputTrackContainer(ctx, std::move(outputTracks));

  return ProcessCode::SUCCESS;
}
