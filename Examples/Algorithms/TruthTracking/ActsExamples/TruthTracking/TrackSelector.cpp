// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TruthTracking/TrackSelector.hpp"

#include "Acts/Utilities/ThrowAssert.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

ActsExamples::TrackSelector::TrackSelector(const Config& config,
                                           Acts::Logging::Level level)
    : BareAlgorithm("TrackSelector", level), m_cfg(config) {
  if (m_cfg.inputTrackParameters.empty()) {
    throw std::invalid_argument("Missing input track parameters collection");
  }
  if (m_cfg.outputTrackParameters.empty()) {
    throw std::invalid_argument("Missing output track parameters collection");
  }
}

ActsExamples::ProcessCode ActsExamples::TrackSelector::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // helper functions to select tracks
  auto within = [](double x, double min, double max) {
    return (min <= x) and (x < max);
  };
  auto isValidTrack = [&](const auto& trk) {
    const auto theta = trk.template get<Acts::eBoundTheta>();
    const auto eta = -std::log(std::tan(theta / 2));
    // define charge selection
    const bool validNeutral = (trk.charge() == 0) and not m_cfg.removeNeutral;
    const bool validCharged = (trk.charge() != 0) and not m_cfg.removeCharged;
    const bool validCharge = validNeutral or validCharged;
    return validCharge and
           within(trk.transverseMomentum(), m_cfg.ptMin, m_cfg.ptMax) and
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

  // prepare input and output containers
  const auto& inputTrackParameters =
      ctx.eventStore.get<TrackParametersContainer>(m_cfg.inputTrackParameters);
  std::optional<std::vector<std::pair<size_t, size_t>>>
      inputTrackParametersTips;
  if (!m_cfg.inputTrackParametersTips.empty()) {
    inputTrackParametersTips =
        ctx.eventStore.get<std::vector<std::pair<size_t, size_t>>>(
            m_cfg.inputTrackParametersTips);
    throw_assert(
        inputTrackParameters.size() == inputTrackParametersTips->size(),
        "Track parameters tips count does not match track parameters count");
  }

  TrackParametersContainer outputTrackParameters;
  std::vector<std::pair<size_t, size_t>> outputTrackParametersTips;
  outputTrackParameters.reserve(inputTrackParameters.size());
  if (inputTrackParametersTips) {
    outputTrackParametersTips.reserve(inputTrackParametersTips->size());
  }

  // copy selected tracks and record initial track index
  for (uint32_t i = 0; i < inputTrackParameters.size(); ++i) {
    const auto& trk = inputTrackParameters[i];
    if (isValidTrack(trk)) {
      outputTrackParameters.push_back(trk);
      if (inputTrackParametersTips) {
        outputTrackParametersTips.push_back((*inputTrackParametersTips)[i]);
      }
    }
  }
  outputTrackParameters.shrink_to_fit();
  outputTrackParametersTips.shrink_to_fit();

  ACTS_DEBUG("event " << ctx.eventNumber << " selected "
                      << outputTrackParameters.size() << " from "
                      << inputTrackParameters.size() << " tracks");

  ctx.eventStore.add(m_cfg.outputTrackParameters,
                     std::move(outputTrackParameters));
  if (inputTrackParametersTips) {
    ctx.eventStore.add(m_cfg.outputTrackParametersTips,
                       std::move(outputTrackParametersTips));
  }
  return ProcessCode::SUCCESS;
}
