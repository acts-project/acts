// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Json/JsonTrackFinderPerformanceWriter.hpp"

#include "ActsPlugins/Json/HistogramJsonConverter.hpp"

#include <fstream>
#include <stdexcept>

#include <nlohmann/json.hpp>

using ActsPlugins::toJson;

namespace ActsExamples {

namespace {

nlohmann::json serializeTrackSummaryTool(const TrackSummaryPlotTool& tool) {
  nlohmann::json arr = nlohmann::json::array();
  for (const auto& [name, prof] : tool.profiles()) {
    arr.push_back(toJson(prof));
  }
  return arr;
}

}  // namespace

JsonTrackFinderPerformanceWriter::JsonTrackFinderPerformanceWriter(
    JsonTrackFinderPerformanceWriter::Config cfg, Acts::Logging::Level lvl)
    : WriterT(cfg.inputTracks, "JsonTrackFinderPerformanceWriter", lvl),
      m_cfg(std::move(cfg)),
      m_collector(
          TrackFinderPerformanceCollector::Config{
              m_cfg.effPlotToolConfig, m_cfg.fakePlotToolConfig,
              m_cfg.duplicationPlotToolConfig, m_cfg.trackSummaryPlotToolConfig,
              m_cfg.trackQualityPlotToolConfig,
              m_cfg.subDetectorTrackSummaryVolumes},
          logger().clone()) {
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing particles input collection");
  }
  if (m_cfg.inputTrackParticleMatching.empty()) {
    throw std::invalid_argument("Missing input track particles matching");
  }
  if (m_cfg.inputParticleTrackMatching.empty()) {
    throw std::invalid_argument("Missing input particle track matching");
  }
  if (m_cfg.inputParticleMeasurementsMap.empty()) {
    throw std::invalid_argument("Missing input measurement particles map");
  }
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing output filename");
  }

  m_inputParticles.initialize(m_cfg.inputParticles);
  m_inputTrackParticleMatching.initialize(m_cfg.inputTrackParticleMatching);
  m_inputParticleTrackMatching.initialize(m_cfg.inputParticleTrackMatching);
  m_inputParticleMeasurementsMap.initialize(m_cfg.inputParticleMeasurementsMap);
}

JsonTrackFinderPerformanceWriter::~JsonTrackFinderPerformanceWriter() = default;

ProcessCode JsonTrackFinderPerformanceWriter::finalize() {
  m_collector.logSummary();

  nlohmann::json root;

  // Efficiency histograms
  {
    nlohmann::json arr = nlohmann::json::array();
    for (const auto& [name, eff] : m_collector.effPlotTool().efficiencies1D()) {
      arr.push_back(toJson(eff));
    }
    for (const auto& [name, eff] : m_collector.effPlotTool().efficiencies2D()) {
      arr.push_back(toJson(eff));
    }
    for (const auto& eff :
         m_collector.effPlotTool().trackEffVsEtaInPtRanges()) {
      arr.push_back(toJson(eff));
    }
    for (const auto& eff :
         m_collector.effPlotTool().trackEffVsPtInAbsEtaRanges()) {
      arr.push_back(toJson(eff));
    }
    root["efficiency"] = std::move(arr);
  }

  // Fake rate histograms
  {
    nlohmann::json arr = nlohmann::json::array();
    for (const auto& [name, hist] : m_collector.fakePlotTool().histograms()) {
      arr.push_back(toJson(hist));
    }
    for (const auto& [name, eff] : m_collector.fakePlotTool().efficiencies()) {
      arr.push_back(toJson(eff));
    }
    root["fake"] = std::move(arr);
  }

  // Duplication rate histograms
  {
    nlohmann::json arr = nlohmann::json::array();
    for (const auto& [name, prof] :
         m_collector.duplicationPlotTool().profiles()) {
      arr.push_back(toJson(prof));
    }
    for (const auto& [name, eff] :
         m_collector.duplicationPlotTool().efficiencies()) {
      arr.push_back(toJson(eff));
    }
    root["duplication"] = std::move(arr);
  }

  // Track summary histograms (global + per-subdetector)
  root["trackSummary"] =
      serializeTrackSummaryTool(m_collector.trackSummaryPlotTool());
  {
    nlohmann::json subDetJson;
    for (const auto& [key, tool] : m_collector.subDetectorSummaryTools()) {
      subDetJson[key] = serializeTrackSummaryTool(tool);
    }
    root["subDetectorSummary"] = std::move(subDetJson);
  }

  // Track quality histograms
  {
    nlohmann::json arr = nlohmann::json::array();
    for (const auto& [name, prof] :
         m_collector.trackQualityPlotTool().profiles()) {
      arr.push_back(toJson(prof));
    }
    root["trackQuality"] = std::move(arr);
  }

  // Summary scalars
  {
    const auto& s = m_collector.stats();
    auto safeRatio = [](std::size_t num, std::size_t den) -> double {
      return den > 0 ? static_cast<double>(num) / static_cast<double>(den)
                     : 0.0;
    };
    root["scalars"] = {
        {"nTotalTracks", s.nTotalTracks},
        {"nTotalMatchedTracks", s.nTotalMatchedTracks},
        {"nTotalFakeTracks", s.nTotalFakeTracks},
        {"nTotalDuplicateTracks", s.nTotalDuplicateTracks},
        {"nTotalParticles", s.nTotalParticles},
        {"nTotalMatchedParticles", s.nTotalMatchedParticles},
        {"nTotalDuplicateParticles", s.nTotalDuplicateParticles},
        {"nTotalFakeParticles", s.nTotalFakeParticles},
        {"eff_tracks", safeRatio(s.nTotalMatchedTracks, s.nTotalTracks)},
        {"fakeratio_tracks", safeRatio(s.nTotalFakeTracks, s.nTotalTracks)},
        {"duplicateratio_tracks",
         safeRatio(s.nTotalDuplicateTracks, s.nTotalTracks)},
        {"eff_particles",
         safeRatio(s.nTotalMatchedParticles, s.nTotalParticles)},
        {"fakeratio_particles",
         safeRatio(s.nTotalFakeParticles, s.nTotalParticles)},
        {"duplicateratio_particles",
         safeRatio(s.nTotalDuplicateParticles, s.nTotalParticles)},
    };
  }

  // Optional matching details
  if (m_cfg.writeMatchingDetails) {
    nlohmann::json arr = nlohmann::json::array();
    for (const auto& r : m_matchingDetails) {
      arr.push_back({{"event_nr", r.eventNr},
                     {"particle_id_vertex_primary", r.vertexPrimary},
                     {"particle_id_vertex_secondary", r.vertexSecondary},
                     {"particle_id_particle", r.particle},
                     {"particle_id_generation", r.generation},
                     {"particle_id_sub_particle", r.subParticle},
                     {"matched", r.isMatched}});
    }
    root["matchingDetails"] = std::move(arr);
  }

  std::ofstream ofstream(m_cfg.filePath);
  if (!ofstream.is_open()) {
    throw std::ios_base::failure("Could not open '" + m_cfg.filePath + "'");
  }
  ofstream << root.dump(2);
  ACTS_INFO("Wrote performance plots to '" << m_cfg.filePath << "'");

  return ProcessCode::SUCCESS;
}

ProcessCode JsonTrackFinderPerformanceWriter::writeT(
    const AlgorithmContext& ctx, const ConstTrackContainer& tracks) {
  const auto& particles = m_inputParticles(ctx);
  const auto& trackParticleMatching = m_inputTrackParticleMatching(ctx);
  const auto& particleTrackMatching = m_inputParticleTrackMatching(ctx);
  const auto& particleMeasurementsMap = m_inputParticleMeasurementsMap(ctx);

  std::lock_guard<std::mutex> lock(m_writeMutex);

  m_collector.fill(ctx.geoContext, tracks, particles, trackParticleMatching,
                   particleTrackMatching, particleMeasurementsMap);

  if (m_cfg.writeMatchingDetails) {
    for (const auto& particle : particles) {
      auto particleId = particle.particleId();
      bool isMatched = particleTrackMatching.contains(particleId) &&
                       particleTrackMatching.at(particleId).track.has_value();
      m_matchingDetails.push_back(
          {static_cast<std::uint32_t>(ctx.eventNumber),
           particleId.vertexPrimary(), particleId.vertexSecondary(),
           particleId.particle(), particleId.generation(),
           particleId.subParticle(), isMatched});
    }
  }

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
