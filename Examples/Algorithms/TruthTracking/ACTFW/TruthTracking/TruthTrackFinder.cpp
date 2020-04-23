// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/TruthTracking/TruthTrackFinder.hpp"

#include <algorithm>
#include <stdexcept>
#include <vector>

#include "ACTFW/EventData/IndexContainers.hpp"
#include "ACTFW/EventData/ProtoTrack.hpp"
#include "ACTFW/EventData/SimParticle.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "ACTFW/Utilities/Range.hpp"

using namespace FW;

TruthTrackFinder::TruthTrackFinder(const Config& cfg, Acts::Logging::Level lvl)
    : BareAlgorithm("TruthTrackFinder", lvl), m_cfg(cfg) {
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input truth particles collection");
  }
  if (m_cfg.inputHitParticlesMap.empty()) {
    throw std::invalid_argument("Missing input hit-particles map collection");
  }
  if (m_cfg.outputProtoTracks.empty()) {
    throw std::invalid_argument("Missing output proto tracks collection");
  }
}

ProcessCode TruthTrackFinder::execute(const AlgorithmContext& ctx) const {
  using HitParticlesMap = IndexMultimap<ActsFatras::Barcode>;

  // prepare input collections
  const auto& particles =
      ctx.eventStore.get<SimParticleContainer>(m_cfg.inputParticles);
  const auto& hitParticlesMap =
      ctx.eventStore.get<HitParticlesMap>(m_cfg.inputHitParticlesMap);
  // compute particle_id -> {hit_id...} map from the
  // hit_id -> {particle_id...} map on the fly.
  const auto& particleHitsMap = invertIndexMultimap(hitParticlesMap);

  // prepare output collection
  ProtoTrackContainer tracks;
  tracks.reserve(particles.size());

  // create prototracks for all input particles
  for (const auto& particle : particles) {
    // find the corresponding hits for this particle
    const auto& hits =
        makeRange(particleHitsMap.equal_range(particle.particleId()));
    // fill hit indices to create the proto track
    ProtoTrack track;
    track.reserve(hits.size());
    for (const auto& hit : hits) {
      track.emplace_back(hit.second);
    }
    // add proto track to the output collection
    tracks.emplace_back(std::move(track));
  }

  ctx.eventStore.add(m_cfg.outputProtoTracks, std::move(tracks));
  return ProcessCode::SUCCESS;
}
