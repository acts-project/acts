// This file is part of the Acts project.
//
// Copyright (C) 2019-2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TruthTracking/TrackModifier.hpp"

#include "Acts/Utilities/ThrowAssert.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

ActsExamples::TrackModifier::TrackModifier(const Config& config,
                                           Acts::Logging::Level level)
    : BareAlgorithm("TrackModifier", level), m_cfg(config) {
  if (m_cfg.inputTrackParameters.empty()) {
    throw std::invalid_argument("Missing input track parameters collection");
  }
  if (m_cfg.outputTrackParameters.empty()) {
    throw std::invalid_argument("Missing output track parameters collection");
  }
}

ActsExamples::ProcessCode ActsExamples::TrackModifier::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // prepare input and output containers
  const auto& inputTrackParameters =
      ctx.eventStore.get<TrackParametersContainer>(m_cfg.inputTrackParameters);
  TrackParametersContainer outputTrackParameters;
  outputTrackParameters.reserve(inputTrackParameters.size());

  for (const auto& trk : inputTrackParameters) {
    auto newTrk = trk;

    std::optional<Acts::BoundSymMatrix> optCov = trk.covariance();

    if (optCov) {
      auto& cov = *optCov;

      if (m_cfg.dropCovariance) {
        cov = cov.diagonal().asDiagonal();
      }
      cov *= m_cfg.covScale;

      optCov = cov;
    }

    newTrk.setCovariance(optCov);

    outputTrackParameters.emplace_back(newTrk);
  }

  ctx.eventStore.add(m_cfg.outputTrackParameters,
                     std::move(outputTrackParameters));
  return ProcessCode::SUCCESS;
}
