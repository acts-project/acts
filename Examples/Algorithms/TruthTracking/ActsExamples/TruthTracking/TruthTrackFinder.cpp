// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TruthTracking/TruthTrackFinder.hpp"

#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/Range.hpp"

#include <algorithm>
#include <stdexcept>
#include <vector>

using namespace ActsExamples;

TruthTrackFinder::TruthTrackFinder(const Config& config,
                                   Acts::Logging::Level level)
    : IAlgorithm("TruthTrackFinder", level), m_cfg(config) {
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input truth particles collection");
  }
  if (m_cfg.inputMeasurementParticlesMap.empty()) {
    throw std::invalid_argument("Missing input hit-particles map collection");
  }
  if (m_cfg.outputProtoTracks.empty()) {
    throw std::invalid_argument("Missing output proto tracks collection");
  }

  m_inputParticles.initialize(m_cfg.inputParticles);
  m_inputMeasurementParticlesMap.initialize(m_cfg.inputMeasurementParticlesMap);
  m_outputProtoTracks.initialize(m_cfg.outputProtoTracks);
}

ProcessCode TruthTrackFinder::execute(const AlgorithmContext& ctx) const {
  // prepare input collections
  const auto& particles = m_inputParticles(ctx);
  const auto& hitParticlesMap = m_inputMeasurementParticlesMap(ctx);
  // compute particle_id -> {hit_id...} map from the
  // hit_id -> {particle_id...} map on the fly.
  const auto& particleHitsMap = invertIndexMultimap(hitParticlesMap);

  // prepare output collection
  ProtoTrackContainer tracks;
  tracks.reserve(particles.size());

  ACTS_VERBOSE("Create prototracks for " << particles.size() << " particles");
  for (const auto& particle : particles) {
    // find the corresponding hits for this particle
    const auto& hits =
        makeRange(particleHitsMap.equal_range(particle.particleId()));
    ACTS_VERBOSE(" - Prototrack from " << hits.size() << " hits");
    // fill hit indices to create the proto track
    ProtoTrack track;
    track.reserve(hits.size());
    for (const auto& hit : hits) {
      track.emplace_back(hit.second);
    }
    // add proto track to the output collection
    tracks.emplace_back(std::move(track));
  }

  m_outputProtoTracks(ctx, std::move(tracks));
  return ProcessCode::SUCCESS;
}
