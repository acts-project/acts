// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Json/JsonTrackWriter.hpp"

#include <fstream>
#include <sstream>

#include <nlohmann/json.hpp>

namespace ActsExamples {

JsonTrackWriter::JsonTrackWriter(const JsonTrackWriter::Config &config,
                                 Acts::Logging::Level level)
    : m_cfg(config),
      m_logger(Acts::getDefaultLogger("JsonTrackWriter", level)) {
  m_inputTracks.initialize(m_cfg.inputTracks);
  m_inputMeasurementParticlesMap.initialize(m_cfg.inputMeasurementParticlesMap);
}

std::string JsonTrackWriter::name() const {
  return "JsonTrackWriter";
}

ProcessCode JsonTrackWriter::write(const AlgorithmContext &ctx) {
  const auto &tracks = m_inputTracks(ctx);
  const auto &measPartMap = m_inputMeasurementParticlesMap(ctx);

  nlohmann::json json = nlohmann::json::array();

  for (const auto &track : tracks) {
    nlohmann::json jTrack;
    jTrack["reference_surface"] = track.referenceSurface().geometryId().value();
    jTrack["track_states"] = nlohmann::json::array();

    for (const auto &state : track.trackStates()) {
      if (!state.hasUncalibratedSourceLink()) {
        continue;
      }
      nlohmann::json jState;
      jState["reference_surface"] =
          state.referenceSurface().geometryId().value();

      auto idx = state.getUncalibratedSourceLink()
                     .template get<IndexSourceLink>()
                     .index();
      jState["measurement_idx"] = idx;

      nlohmann::json jParticles = nlohmann::json::array();
      auto [begin, end] = measPartMap.equal_range(idx);
      for (auto it = begin; it != end; ++it) {
        jParticles.push_back(it->second.value());
      }
      jState["measurement_particles"] = jParticles;
      jTrack["track_states"].push_back(jState);
    }

    json.push_back(jTrack);
  }

  std::stringstream fname;
  fname << "event" << ctx.eventNumber << "-tracks.json";
  std::ofstream outfile(fname.str());
  outfile << json.dump(4);

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
