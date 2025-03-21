// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/EDM4hep/EDM4hepMultiTrajectoryOutputConverter.hpp"

#include "ActsExamples/Io/EDM4hep/EDM4hepUtil.hpp"
#include "ActsExamples/Validation/TrackClassification.hpp"

#include <stdexcept>

#include <edm4hep/TrackCollection.h>
#include <podio/Frame.h>

namespace ActsExamples {

EDM4hepMultiTrajectoryOutputConverter::EDM4hepMultiTrajectoryOutputConverter(
    const EDM4hepMultiTrajectoryOutputConverter::Config& config,
    Acts::Logging::Level level)
    : EDM4hepOutputConverter("EDM4hepMultiTrajectoryOutputConverter", level),
      m_cfg(config) {
  if (m_cfg.inputTrajectories.empty()) {
    throw std::invalid_argument("Missing input trajectories collection");
  }

  if (m_cfg.inputMeasurementParticlesMap.empty()) {
    throw std::invalid_argument{"Missing input hit to particle map"};
  }

  if (m_cfg.outputTracks.empty()) {
    throw std::invalid_argument("Missing output tracks collection");
  }

  m_inputMeasurementParticlesMap.initialize(m_cfg.inputMeasurementParticlesMap);
  m_outputTracks.initialize(m_cfg.outputTracks);
  m_inputTrajectories.initialize(m_cfg.inputTrajectories);
}

ProcessCode EDM4hepMultiTrajectoryOutputConverter::execute(
    const AlgorithmContext& context) const {
  podio::Frame frame{};

  const auto& hitParticlesMap = m_inputMeasurementParticlesMap(context);
  const auto& trajectories = m_inputTrajectories(context);

  edm4hep::TrackCollection trackCollection;

  for (const auto& from : trajectories) {
    for (const auto& trackTip : from.tips()) {
      auto to = trackCollection.create();
      EDM4hepUtil::writeTrajectory(context.geoContext, m_cfg.Bz, from, to,
                                   trackTip, m_cfg.particleHypothesis,
                                   hitParticlesMap);
    }
  }

  m_outputTracks(context, std::move(trackCollection));

  return ProcessCode::SUCCESS;
}

std::vector<std::string> EDM4hepMultiTrajectoryOutputConverter::collections()
    const {
  return {m_cfg.outputTracks};
}

}  // namespace ActsExamples
