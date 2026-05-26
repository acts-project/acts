// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Json/JsonTrackFitterPerformanceWriter.hpp"

#include "Acts/Utilities/Helpers.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsFatras/EventData/Barcode.hpp"
#include "ActsPlugins/Json/HistogramJsonConverter.hpp"

#include <cstddef>
#include <fstream>
#include <stdexcept>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::phi;
using ActsPlugins::toJson;

namespace ActsExamples {

JsonTrackFitterPerformanceWriter::JsonTrackFitterPerformanceWriter(
    JsonTrackFitterPerformanceWriter::Config config, Acts::Logging::Level level)
    : WriterT(config.inputTracks, "JsonTrackFitterPerformanceWriter", level),
      m_cfg(std::move(config)),
      m_resPlotTool(m_cfg.resPlotToolConfig, level),
      m_effPlotTool(m_cfg.effPlotToolConfig, level),
      m_trackSummaryPlotTool(m_cfg.trackSummaryPlotToolConfig, level) {
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing particles input collection");
  }
  if (m_cfg.inputTrackParticleMatching.empty()) {
    throw std::invalid_argument("Missing input track particles matching");
  }
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing output filename");
  }

  m_inputParticles.initialize(m_cfg.inputParticles);
  m_inputTrackParticleMatching.initialize(m_cfg.inputTrackParticleMatching);
}

JsonTrackFitterPerformanceWriter::~JsonTrackFitterPerformanceWriter() = default;

ProcessCode JsonTrackFitterPerformanceWriter::finalize() {
  nlohmann::json root;

  // Residual histograms (1D, 2D vs eta/pt, 3D vs eta-phi and eta-pt)
  {
    nlohmann::json arr = nlohmann::json::array();
    for (const auto& [name, hist] : m_resPlotTool.res()) {
      arr.push_back(toJson(hist));
    }
    root["residuals"] = std::move(arr);
  }
  {
    nlohmann::json arr = nlohmann::json::array();
    for (const auto& [name, hist] : m_resPlotTool.resVsEta()) {
      arr.push_back(toJson(hist));
    }
    root["residualsVsEta"] = std::move(arr);
  }
  {
    nlohmann::json arr = nlohmann::json::array();
    for (const auto& [name, hist] : m_resPlotTool.resVsPt()) {
      arr.push_back(toJson(hist));
    }
    root["residualsVsPt"] = std::move(arr);
  }
  {
    nlohmann::json arr = nlohmann::json::array();
    for (const auto& [name, hist] : m_resPlotTool.resVsEtaPhi()) {
      arr.push_back(toJson(hist));
    }
    root["residualsVsEtaPhi"] = std::move(arr);
  }
  {
    nlohmann::json arr = nlohmann::json::array();
    for (const auto& [name, hist] : m_resPlotTool.resVsEtaPt()) {
      arr.push_back(toJson(hist));
    }
    root["residualsVsEtaPt"] = std::move(arr);
  }

  // Pull histograms
  {
    nlohmann::json arr = nlohmann::json::array();
    for (const auto& [name, hist] : m_resPlotTool.pull()) {
      arr.push_back(toJson(hist));
    }
    root["pulls"] = std::move(arr);
  }
  {
    nlohmann::json arr = nlohmann::json::array();
    for (const auto& [name, hist] : m_resPlotTool.pullVsEta()) {
      arr.push_back(toJson(hist));
    }
    root["pullsVsEta"] = std::move(arr);
  }
  {
    nlohmann::json arr = nlohmann::json::array();
    for (const auto& [name, hist] : m_resPlotTool.pullVsPt()) {
      arr.push_back(toJson(hist));
    }
    root["pullsVsPt"] = std::move(arr);
  }
  {
    nlohmann::json arr = nlohmann::json::array();
    for (const auto& [name, hist] : m_resPlotTool.pullVsEtaPhi()) {
      arr.push_back(toJson(hist));
    }
    root["pullsVsEtaPhi"] = std::move(arr);
  }
  {
    nlohmann::json arr = nlohmann::json::array();
    for (const auto& [name, hist] : m_resPlotTool.pullVsEtaPt()) {
      arr.push_back(toJson(hist));
    }
    root["pullsVsEtaPt"] = std::move(arr);
  }

  // Efficiency histograms
  {
    nlohmann::json arr = nlohmann::json::array();
    for (const auto& [name, eff] : m_effPlotTool.efficiencies1D()) {
      arr.push_back(toJson(eff));
    }
    for (const auto& [name, eff] : m_effPlotTool.efficiencies2D()) {
      arr.push_back(toJson(eff));
    }
    for (const auto& eff : m_effPlotTool.trackEffVsEtaInPtRanges()) {
      arr.push_back(toJson(eff));
    }
    for (const auto& eff : m_effPlotTool.trackEffVsPtInAbsEtaRanges()) {
      arr.push_back(toJson(eff));
    }
    root["efficiency"] = std::move(arr);
  }

  // Track summary histograms
  {
    nlohmann::json arr = nlohmann::json::array();
    for (const auto& [name, prof] : m_trackSummaryPlotTool.profiles()) {
      arr.push_back(toJson(prof));
    }
    root["trackSummary"] = std::move(arr);
  }

  std::ofstream ofstream(m_cfg.filePath);
  if (!ofstream.is_open()) {
    throw std::ios_base::failure("Could not open '" + m_cfg.filePath + "'");
  }
  ofstream << root.dump(2);
  ACTS_INFO("Wrote performance plots to '" << m_cfg.filePath << "'");

  return ProcessCode::SUCCESS;
}

ProcessCode JsonTrackFitterPerformanceWriter::writeT(
    const AlgorithmContext& ctx, const ConstTrackContainer& tracks) {
  const auto& particles = m_inputParticles(ctx);
  const auto& trackParticleMatching = m_inputTrackParticleMatching(ctx);

  std::vector<SimBarcode> reconParticleIds;
  reconParticleIds.reserve(particles.size());

  std::lock_guard<std::mutex> lock(m_writeMutex);

  for (const auto& track : tracks) {
    if (!track.hasReferenceSurface()) {
      ACTS_WARNING("No fitted track parameters.");
      continue;
    }
    Acts::BoundTrackParameters fittedParameters =
        track.createParametersAtReference();

    auto imatched = trackParticleMatching.find(track.index());
    if (imatched == trackParticleMatching.end()) {
      ACTS_DEBUG("No truth particle associated with this track, index = "
                 << track.index() << " tip index = " << track.tipIndex());
      continue;
    }
    const auto& particleMatch = imatched->second;
    if (!particleMatch.particle.has_value()) {
      ACTS_DEBUG("No truth particle associated with this track.");
      continue;
    }

    SimBarcode majorityParticleId = particleMatch.particle.value();
    auto ip = particles.find(majorityParticleId);
    if (ip == particles.end()) {
      ACTS_DEBUG("Majority particle not found in the particles collection.");
      continue;
    }

    reconParticleIds.push_back(ip->particleId());
    m_resPlotTool.fill(ctx.geoContext, ip->initialState(), fittedParameters);
    m_trackSummaryPlotTool.fill(fittedParameters, track.nTrackStates(),
                                track.nMeasurements(), track.nOutliers(),
                                track.nHoles(), track.nSharedHits());
  }

  for (const auto& particle : particles) {
    bool isReconstructed =
        Acts::rangeContainsValue(reconParticleIds, particle.particleId());

    double minDeltaR = -1;
    for (const auto& closeParticle : particles) {
      if (closeParticle.particleId() == particle.particleId()) {
        continue;
      }
      double p_phi = phi(particle.direction());
      double p_eta = eta(particle.direction());
      double c_phi = phi(closeParticle.direction());
      double c_eta = eta(closeParticle.direction());
      double distance = std::sqrt(std::pow(p_phi - c_phi, 2) +
                                  std::pow(p_eta - c_eta, 2));
      if (minDeltaR == -1 || distance < minDeltaR) {
        minDeltaR = distance;
      }
    }
    m_effPlotTool.fill(ctx.geoContext, particle.initialState(), minDeltaR,
                       isReconstructed);
  }

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
