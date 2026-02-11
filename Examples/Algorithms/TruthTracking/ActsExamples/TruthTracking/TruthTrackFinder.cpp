// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TruthTracking/TruthTrackFinder.hpp"

#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Utilities/Range.hpp"

#include <ostream>
#include <stdexcept>
#include <utility>

namespace ActsExamples {

TruthTrackFinder::TruthTrackFinder(const Config& config,
                                   Acts::Logging::Level level)
    : IAlgorithm("TruthTrackFinder", level), m_cfg(config) {
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input truth particles collection");
  }
  if (m_cfg.inputParticleMeasurementsMap.empty()) {
    throw std::invalid_argument("Missing input hit-particles map collection");
  }
  if (m_cfg.outputProtoTracks.empty()) {
    throw std::invalid_argument("Missing output proto tracks collection");
  }
  if (m_cfg.inputMeasurements.empty()) {
    throw std::invalid_argument("Missing input measurements collection");
  }
  if (m_cfg.inputSimHits.empty()) {
    throw std::invalid_argument("Missing input simulated hits collection");
  }
  if (m_cfg.inputMeasurementSimHitsMap.empty()) {
    throw std::invalid_argument(
        "Missing input simulated hits measurements map");
  }

  m_inputParticles.initialize(m_cfg.inputParticles);
  m_inputParticleMeasurementsMap.initialize(m_cfg.inputParticleMeasurementsMap);
  m_outputProtoTracks.initialize(m_cfg.outputProtoTracks);
  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_inputSimHits.initialize(m_cfg.inputSimHits);
  m_inputMeasurementSimHitsMap.initialize(m_cfg.inputMeasurementSimHitsMap);
}

ProcessCode TruthTrackFinder::execute(const AlgorithmContext& ctx) const {
  // prepare input collections
  const auto& particles = m_inputParticles(ctx);
  const auto& particleMeasurementsMap = m_inputParticleMeasurementsMap(ctx);
  const auto& measurementsIn = m_inputMeasurements(ctx);
  const auto& simHits = m_inputSimHits(ctx);
  const auto& measurementSimHitsMap = m_inputMeasurementSimHitsMap(ctx);

  // prepare output collection
  ProtoTrackContainer tracks;
  tracks.reserve(particles.size());

  ACTS_VERBOSE("Create prototracks for " << particles.size() << " particles");
  for (const auto& [i, particle] : Acts::enumerate(particles)) {
    // find the corresponding hits for this particle
    const auto& measurements =
        makeRange(particleMeasurementsMap.equal_range(particle.particleId()));
    ACTS_VERBOSE(" - Prototrack from " << measurements.size()
                                       << " measurements for particle "
                                       << particle);
    // fill hit indices to create the proto track
    ProtoTrack track;
    std::vector<const SimHit*> hits;
    track.reserve(measurements.size());
    hits.reserve(measurements.size());
    for (const auto& [barcode, index] : measurements) {
      ConstVariableBoundMeasurementProxy measurement =
          measurementsIn.getMeasurement(index);
      ACTS_VERBOSE("   - Measurement " << index << " with barcode " << barcode
                                       << " at " << measurement.geometryId());
      const auto simHitMapIt = measurementSimHitsMap.find(index);
      if (simHitMapIt == measurementSimHitsMap.end()) {
        ACTS_WARNING("No sim hit found for measurement index " << index);
        continue;
      }

      const auto simHitIdxIt = measurementSimHitsMap.nth(simHitMapIt->second);
      if (simHitIdxIt == measurementSimHitsMap.end()) {
        ACTS_WARNING("No sim hit found for index " << simHitMapIt->second);
        continue;
      }

      const auto simHitIt = simHits.nth(simHitIdxIt->second);
      if (simHitIt == simHits.end()) {
        ACTS_WARNING("No sim hit found for index " << simHitIdxIt->second);
        continue;
      }

      const auto& simHit = *simHitIt;

      track.emplace_back(index);
      hits.emplace_back(&simHit);
    }

    std::vector<std::size_t> indices;
    indices.resize(hits.size());
    std::iota(indices.begin(), indices.end(), 0);

    // std::ranges::sort(indices, [&hits](std::size_t a, std::size_t b) {
    //   return hits[a]->time() < hits[b]->time();
    // });
    ProtoTrack sortedTrack;
    for (const auto& idx : indices) {
      sortedTrack.emplace_back(track[idx]);
    }

    // add proto track to the output collection
    tracks.emplace_back(std::move(sortedTrack));
  }

  m_outputProtoTracks(ctx, std::move(tracks));
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
