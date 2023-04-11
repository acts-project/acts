// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TruthTracking/TrackModifier.hpp"

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/Utilities/ThrowAssert.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

ActsExamples::TrackModifier::TrackModifier(const Config& config,
                                           Acts::Logging::Level level)
    : IAlgorithm("TrackModifier", level), m_cfg(config) {
  if (m_cfg.inputTrajectories.empty() == m_cfg.inputTrackParameters.empty()) {
    throw std::invalid_argument(
        "Exactly one of trajectories or track parameters input must be set");
  }
  if (m_cfg.outputTrajectories.empty() == m_cfg.outputTrackParameters.empty()) {
    throw std::invalid_argument(
        "Exactly one of trajectories or track parameters output must be set");
  }
  if (m_cfg.inputTrajectories.empty() != m_cfg.outputTrajectories.empty()) {
    throw std::invalid_argument(
        "Input and output for trajectories and track parameters have to be "
        "used consistently");
  }

  m_inputTrackParameters.maybeInitialize(m_cfg.inputTrackParameters);
  m_inputTrajectories.maybeInitialize(m_cfg.inputTrackParameters);
  m_outputTrackParameters.maybeInitialize(m_cfg.outputTrackParameters);
  m_outputTrajectories.maybeInitialize(m_cfg.outputTrajectories);
}

ActsExamples::ProcessCode ActsExamples::TrackModifier::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  auto modifyTrack = [this](auto trk) {
    {
      auto& params = trk.parameters();

      if (m_cfg.killTime) {
        params[Acts::eBoundTime] = 0;
      }
    }

    {
      auto& optCov = trk.covariance();

      if (optCov) {
        auto& cov = *optCov;

        if (m_cfg.dropCovariance) {
          cov = Acts::BoundSymMatrix(cov.diagonal().asDiagonal());
        }
        if (m_cfg.covScale != 1) {
          cov *= m_cfg.covScale;
        }
        if (m_cfg.killTime) {
          cov.row(Acts::eBoundTime).setZero();
          cov.col(Acts::eBoundTime).setZero();
          cov(Acts::eBoundTime, Acts::eBoundTime) = 1;
        }
      }
    }

    return trk;
  };

  if (!m_cfg.inputTrackParameters.empty()) {
    const auto& inputTrackParameters = m_inputTrackParameters(ctx);
    TrackParametersContainer outputTrackParameters;
    outputTrackParameters.reserve(inputTrackParameters.size());

    for (uint32_t i = 0; i < inputTrackParameters.size(); ++i) {
      const auto& trk = inputTrackParameters[i];
      outputTrackParameters.push_back(modifyTrack(trk));
    }

    m_outputTrackParameters(ctx, std::move(outputTrackParameters));
  } else if (!m_cfg.inputTrajectories.empty()) {
    const auto& inputTrajectories = m_inputTrajectories(ctx);
    TrajectoriesContainer outputTrajectories;
    outputTrajectories.reserve(inputTrajectories.size());

    for (const auto& trajectories : inputTrajectories) {
      std::vector<Acts::MultiTrajectoryTraits::IndexType> tips;
      tips.reserve(trajectories.tips().size());
      Trajectories::IndexedParameters parameters;

      for (auto tip : trajectories.tips()) {
        if (!trajectories.hasTrackParameters(tip)) {
          continue;
        }
        const auto& trk = trajectories.trackParameters(tip);
        tips.push_back(tip);
        parameters.emplace(tip, modifyTrack(trk));
      }

      outputTrajectories.emplace_back(trajectories.multiTrajectory(), tips,
                                      parameters);
    }

    m_outputTrajectories(ctx, std::move(outputTrajectories));
  }

  return ProcessCode::SUCCESS;
}
