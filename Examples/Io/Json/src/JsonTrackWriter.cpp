// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Json/JsonTrackWriter.hpp"

#include "Acts/Surfaces/Surface.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <fstream>

const static std::map<ActsExamples::TrackMatchClassification, const char *>
    tmcToStr{{ActsExamples::TrackMatchClassification::Unknown, "unknown"},
             {ActsExamples::TrackMatchClassification::Matched, "matched"},
             {ActsExamples::TrackMatchClassification::Duplicate, "duplicate"},
             {ActsExamples::TrackMatchClassification::Fake, "fake"}};

namespace ActsExamples {

JsonTrackWriter::JsonTrackWriter(const JsonTrackWriter::Config &config,
                                 Acts::Logging::Level level)
    : m_cfg(config),
      m_logger(Acts::getDefaultLogger("JsonTrackWriter", level)) {
  m_inputTracks.initialize(m_cfg.inputTracks);
  m_inputMeasurementParticlesMap.initialize(m_cfg.inputMeasurementParticlesMap);
  m_inputClusters.maybeInitialize(m_cfg.inputClusters);
  m_inputTrackParticleMatching.maybeInitialize(
      m_cfg.inputTrackParticleMatching);
}

JsonTrackWriter::~JsonTrackWriter() {}

std::string JsonTrackWriter::name() const {
  return "JsonTrackWriter";
}

ProcessCode JsonTrackWriter::write(const AlgorithmContext &ctx) {
  auto boundToGlobPos = [&](const Acts::Surface &surface,
                            const Acts::BoundVector &parameters) {
    using namespace Acts;
    return surface.localToGlobal(
        ctx.geoContext, parameters.segment<2>(eBoundLoc0),
        makeDirectionFromPhiTheta(parameters[eBoundPhi],
                                  parameters[eBoundTheta]));
  };

  const auto &tracks = m_inputTracks(ctx);
  const auto &measPartMap = m_inputMeasurementParticlesMap(ctx);

  TrackParticleMatching emptyMatchingInfo;
  const auto &matchingInfo = m_inputTrackParticleMatching.isInitialized()
                                 ? m_inputTrackParticleMatching(ctx)
                                 : emptyMatchingInfo;

  ClusterContainer emptyClusters;
  const auto &clusters =
      m_inputClusters.isInitialized() ? m_inputClusters(ctx) : emptyClusters;

  nlohmann::json jTracks = nlohmann::json::array();

  for (const auto &track : tracks) {
    nlohmann::json jTrack;
    jTrack["index"] = track.index();
    jTrack["reference_surface"] =
        track.hasReferenceSurface()
            ? track.referenceSurface().geometryId().value()
            : -1;
    jTrack["track_states"] = nlohmann::json::array();
    if (track.hasReferenceSurface()) {
      jTrack["track_parameters_global"] =
          boundToGlobPos(track.referenceSurface(), track.parameters());
    }

    for (const auto &state : track.trackStatesReversed()) {
      if (!state.hasUncalibratedSourceLink()) {
        continue;
      }
      nlohmann::json jState;
      jState["reference_surface"] =
          state.hasReferenceSurface()
              ? state.referenceSurface().geometryId().value()
              : -1;

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

      if (!clusters.empty()) {
        jState["cluster_position"] = clusters.at(idx).globalPosition;
      }

      if (state.hasReferenceSurface()) {
        const auto &srf = state.referenceSurface();
        if (state.hasPredicted()) {
          jState["predicted_position"] = boundToGlobPos(srf, state.predicted());
        }
        if (state.hasFiltered()) {
          jState["filtered_position"] = boundToGlobPos(srf, state.filtered());
        }
        if (state.hasSmoothed()) {
          jState["smoothed_position"] = boundToGlobPos(srf, state.smoothed());
        }
      }

      jTrack["track_states"].push_back(jState);
    }

    if (!matchingInfo.empty()) {
      if (matchingInfo.contains(track.index())) {
        const auto &trackMatchInfo = matchingInfo.at(track.index());
        jTrack["matching_info"]["matched"] = true;
        jTrack["matching_info"]["matched_particle"] =
            trackMatchInfo.particle ? trackMatchInfo.particle->value() : -1;
        jTrack["matching_info"]["majority_particle"] =
            trackMatchInfo.contributingParticles.front().particleId.value();
        jTrack["matching_info"]["matching_probability"] =
            trackMatchInfo.majorityMatchingProbability;
        jTrack["matching_info"]["classification"] =
            tmcToStr.at(trackMatchInfo.classification);
      } else {
        jTrack["matching_info"]["matched"] = false;
      }
    }

    jTracks.push_back(jTrack);
  }

  std::ofstream os(perEventFilepath(m_cfg.outputDir, m_cfg.outputStem + ".json",
                                    ctx.eventNumber));
  os << jTracks;

  return ProcessCode::SUCCESS;
}

ProcessCode JsonTrackWriter::finalize() {
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
