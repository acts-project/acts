// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "TrackParametersPrinter.hpp"

#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

ActsExamples::TrackParametersPrinter::TrackParametersPrinter(
    const Config& cfg, Acts::Logging::Level level)
    : IAlgorithm("TrackParametersPrinter", level), m_cfg(cfg) {
  if (m_cfg.inputTrackParameters.empty()) {
    throw std::invalid_argument(
        "Input track parameters collection is not configured");
  }

  m_inputTrackParameters.initialize(m_cfg.inputTrackParameters);
}

ActsExamples::ProcessCode ActsExamples::TrackParametersPrinter::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  const auto& trackParameters = m_inputTrackParameters(ctx);

  ACTS_INFO("event " << ctx.eventNumber << " collection '"
                     << m_cfg.inputTrackParameters << "' contains "
                     << trackParameters.size() << " track parameters");
  std::size_t i = 0;
  for (const auto& params : trackParameters) {
    ACTS_INFO("  track " << i++ << "\n" << params);
  }
  return ProcessCode::SUCCESS;
}
