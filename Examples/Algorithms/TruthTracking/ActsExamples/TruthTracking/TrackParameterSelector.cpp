// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TruthTracking/TrackParameterSelector.hpp"

#include "Acts/Utilities/ThrowAssert.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

ActsExamples::TrackParameterSelector::TrackParameterSelector(
    const Config& config, Acts::Logging::Level level)
    : IAlgorithm("TrackParameterSelector", level), m_cfg(config) {
  if (m_cfg.inputTrackParameters.empty()) {
    throw std::invalid_argument("Missing input track parameters");
  }
  if (m_cfg.outputTrackParameters.empty()) {
    throw std::invalid_argument("Missing output track parameters");
  }

  m_inputTrackParameters.initialize(m_cfg.inputTrackParameters);
  m_outputTrackParameters.initialize(m_cfg.outputTrackParameters);
}

ActsExamples::ProcessCode ActsExamples::TrackParameterSelector::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // helper functions to select tracks
  auto within = [](double x, double min, double max) {
    return (min <= x) and (x < max);
  };
  auto isValidTrack = [&](const auto& trk) {
    const auto theta = trk.template get<Acts::eBoundTheta>();
    const auto eta = -std::log(std::tan(theta / 2));
    // define charge selection
    return within(trk.transverseMomentum(), m_cfg.ptMin, m_cfg.ptMax) and
           within(std::abs(eta), m_cfg.absEtaMin, m_cfg.absEtaMax) and
           within(eta, m_cfg.etaMin, m_cfg.etaMax) and
           within(trk.template get<Acts::eBoundPhi>(), m_cfg.phiMin,
                  m_cfg.phiMax) and
           within(trk.template get<Acts::eBoundLoc0>(), m_cfg.loc0Min,
                  m_cfg.loc0Max) and
           within(trk.template get<Acts::eBoundLoc1>(), m_cfg.loc1Min,
                  m_cfg.loc1Max) and
           within(trk.template get<Acts::eBoundTime>(), m_cfg.timeMin,
                  m_cfg.timeMax);
  };

  const auto& inputTrackParameters = m_inputTrackParameters(ctx);
  TrackParametersContainer outputTrackParameters;
  outputTrackParameters.reserve(inputTrackParameters.size());

  // copy selected tracks and record initial track index
  for (uint32_t i = 0; i < inputTrackParameters.size(); ++i) {
    const auto& trk = inputTrackParameters[i];
    if (isValidTrack(trk)) {
      outputTrackParameters.push_back(trk);
    }
  }
  outputTrackParameters.shrink_to_fit();

  ACTS_DEBUG("event " << ctx.eventNumber << " selected "
                      << outputTrackParameters.size() << " from "
                      << inputTrackParameters.size()
                      << " tracks in track parameters");

  m_outputTrackParameters(ctx, std::move(outputTrackParameters));

  return ProcessCode::SUCCESS;
}
